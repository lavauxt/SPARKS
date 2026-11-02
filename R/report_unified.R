#' Build the unified, dynamic SPARKS HTML report
#'
#' Produces **one self-contained HTML file** that gathers every result of a
#' pipeline run into tabs: an overview, QC, an interactive embedding explorer
#' (lasso selection, colour by metadata or gene), gene violin and dot plots,
#' markers / differential expression with a live volcano plot, enrichment,
#' composition, every result table, the pipeline's figures, the file index and
#' run information. The file needs no server and no network: plotly.js is
#' inlined once, data are embedded as JSON, and each tab is drawn lazily when
#' opened.
#'
#' Figures are drawn from compact "panels" extracted from the merged Seurat
#' objects ([sparks_panel()]): a down-sampled embedding, a handful of metadata
#' columns and the expression of a limited gene set. Anything larger -- full
#' differential expression on an arbitrary selection, re-clustering, any gene --
#' is what the Shiny explorer ([sparks_explorer()]) is for.
#'
#' @param results_dir Pipeline results directory.
#' @param output_file Output HTML path.
#' @param title Report title (default derived from the folder name).
#' @param max_cells Cells embedded per dataset (random down-sampling above).
#' @param max_genes Genes embedded per dataset for the violin / dot plots.
#' @param max_table_rows Row cap per embedded table (the files on disk are
#'   never truncated).
#' @param embed_images Embed the pipeline's PNG figures as base64.
#' @param image_budget_mb Total size of embedded images.
#' @param max_image_mb Largest single image that is embedded.
#' @param bundle Optional pre-read [sparks_read_results()] bundle.
#' @param quiet Suppress progress messages.
#' @return \code{output_file}, invisibly.
#' @export
#' @examples
#' \dontrun{
#' sparks_report("results")
#'
#' # Try it without running the pipeline:
#' demo <- tempfile("demo_"); sparks_demo_results(demo)
#' sparks_report(demo)
#' }
sparks_report <- function(results_dir,
                          output_file = file.path(results_dir, "SPARKS_all_results.html"),
                          title = NULL, max_cells = 20000L, max_genes = 60L,
                          max_table_rows = 3000L, embed_images = TRUE,
                          image_budget_mb = 25, max_image_mb = 1.5,
                          bundle = NULL, quiet = FALSE) {
  missing_pkgs <- .report_missing_packages()
  if (length(missing_pkgs)) {
    stop("sparks_report() requires package(s): ", paste(missing_pkgs, collapse = ", "),
         ". Install with install.packages(c(",
         paste0("'", missing_pkgs, "'", collapse = ", "), ")).", call. = FALSE)
  }
  say <- function(...) if (!isTRUE(quiet)) message("[sparks_report] ", ...)
  bundle <- bundle %||% sparks_read_results(results_dir)
  say("indexed ", nrow(bundle$catalog), " files")

  datasets <- .report_datasets(bundle, max_cells, max_genes, say)
  tables <- .report_tables(bundle, max_table_rows, say)
  figures <- .report_figures(bundle, datasets$sample_counts, say)
  images <- if (isTRUE(embed_images)) .report_images(bundle, image_budget_mb, max_image_mb, say)
            else .report_images(bundle, 0, 0, say)
  data <- list(
    meta = list(
      title = title %||% paste("SPARKS report -", basename(bundle$dir)),
      generated = format(Sys.time(), "%Y-%m-%d %H:%M"),
      version = as.character(tryCatch(utils::packageVersion("SPARKS"), error = function(e) "dev")),
      results_dir_name = basename(bundle$dir),
      species = bundle$config$pipeline$species_target %||% NA_character_,
      demo = file.exists(file.path(bundle$dir, "DEMO_DATA.txt"))
    ),
    theme_map = as.list(.sparks_theme_map()),
    kind_labels = as.list(.KIND_LABELS),
    summary = .report_summary(bundle, datasets, tables),
    datasets = datasets$datasets,
    tables = lapply(tables, `[[`, "json"),
    figures = lapply(seq_along(figures), function(i) {
      f <- figures[[i]]
      f$spec <- sprintf("@@SPARKS_FIG_%d@@", i)
      f$spec_json <- NULL
      f
    }),
    images = images,
    files = .report_files(bundle),
    run = .report_run(bundle)
  )
  json <- .report_json(data, vapply(figures, `[[`, character(1), "spec_json"))
  html <- .report_html(data$meta, json)
  make_dir(dirname(output_file))
  con <- file(output_file, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw(enc2utf8(html)), con)
  say("wrote ", output_file, " (", format(structure(file.size(output_file), class = "object_size"),
                                          units = "auto"), ")")
  invisible(output_file)
}

.report_missing_packages <- function() {
  need <- c("htmltools", "jsonlite", "plotly", "base64enc")
  need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
}

.KIND_LABELS <- c(
  markers_all = "Cluster markers (one vs rest)", deg_sig = "Significant DEGs",
  deg_counts = "DEG counts", cluster_markers = "Positive cluster markers",
  top_expressed = "Top expressed genes", pseudobulk_deg = "Pseudobulk DEG (edgeR)",
  pseudobulk_counts = "Pseudobulk counts", avg_expr = "Average expression",
  proportions = "Cell proportions", prop_test = "Proportion tests",
  gsea = "GSEA", gsea_status = "GSEA status", escape = "Escape signature scores",
  correlation = "Correlations", other_table = "Other tables",
  explorer = "Explorer analyses"
)

# ---- Datasets --------------------------------------------------------------------

#' Find the Seurat objects saved by the pipeline
#' @param bundle A \code{sparks_bundle}.
#' @return Data frame: label, kind (merged / subset), comparison, path.
#' @keywords internal
.discover_result_objects <- function(bundle) {
  cat <- bundle$catalog
  obj <- cat[cat$ext == "rds" & cat$category %in% "RData" &
               grepl("^(Merged|Subset)_", cat$name), , drop = FALSE]
  if (!nrow(obj)) {
    return(data.frame(label = character(), kind = character(), comparison = character(),
                      path = character(), stringsAsFactors = FALSE))
  }
  data.frame(
    label = sub("\\.rds$", "", sub("^(Merged|Subset)_", "", obj$name)),
    kind = ifelse(startsWith(obj$name, "Merged_"), "merged", "subset"),
    comparison = obj$comparison, path = obj$path, stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.arr <- function(x) I(as.vector(x))

.report_gene_pick <- function(bundle, obj, comparison, max_genes) {
  cfg_genes <- unlist(bundle$config$genes$genes_to_plot, use.names = FALSE)
  markers <- sparks_table(bundle, c("markers_all", "cluster_markers"), comparison = comparison)
  top <- character()
  if (nrow(markers) && all(c("gene", "cluster", "avg_log2FC") %in% names(markers))) {
    first <- markers[markers$.file == markers$.file[[1L]], , drop = FALSE]
    top <- unlist(lapply(split(first, first$cluster), function(m) {
      utils::head(m$gene[order(-m$avg_log2FC)], 3L)
    }), use.names = FALSE)
  }
  vf <- tryCatch(SeuratObject::VariableFeatures(obj), error = function(e) character())
  genes <- unique(c(cfg_genes, top, utils::head(vf, max_genes)))
  genes <- genes[genes %in% rownames(obj)]
  list(genes = utils::head(genes, max_genes),
       default = utils::head(unique(c(intersect(cfg_genes, genes), top)), 12L))
}

.report_datasets <- function(bundle, max_cells, max_genes, say) {
  objs <- .discover_result_objects(bundle)
  out <- list(); sample_counts <- NULL; n_genes <- NA_integer_; clusters <- NA_integer_
  if (nrow(objs)) {
    for (i in seq_len(nrow(objs))) {
      path <- objs$path[[i]]
      say("reading ", objs$label[[i]])
      obj <- tryCatch(readRDS(path), error = function(e) {
        say("  skipped (", conditionMessage(e), ")"); NULL })
      if (is.null(obj) || !inherits(obj, "Seurat") || !length(obj@reductions)) next
      pick <- .report_gene_pick(bundle, obj, objs$comparison[[i]], max_genes)
      panel <- tryCatch(
        sparks_panel(obj, reduction = if ("umap" %in% names(obj@reductions)) "umap" else
          names(obj@reductions)[[1L]], max_cells = max_cells, genes = pick$genes,
          max_genes = max_genes),
        error = function(e) { say("  no panel (", conditionMessage(e), ")"); NULL })
      if (is.null(panel)) next
      if (is.null(sample_counts) && "orig.ident" %in% colnames(obj[[]])) {
        sample_counts <- as.data.frame(table(sample = as.character(obj$orig.ident)),
                                       responseName = "cells", stringsAsFactors = FALSE)
      }
      if (is.na(n_genes)) n_genes <- nrow(obj)
      if (is.na(clusters) && "seurat_clusters" %in% colnames(obj[[]])) {
        clusters <- length(unique(obj$seurat_clusters))
      }
      out[[length(out) + 1L]] <- list(
        id = paste0("ds", length(out) + 1L), label = objs$label[[i]],
        kind = objs$kind[[i]], comparison = objs$comparison[[i]],
        n_total = panel$n_total, reduction = panel$reduction,
        cells = .arr(panel$cells), x = .arr(panel$x), y = .arr(panel$y),
        meta = lapply(panel$meta, function(m) {
          if (identical(m$type, "cat")) list(type = "cat", levels = .arr(m$levels), codes = .arr(m$codes))
          else list(type = "num", values = .arr(m$values))
        }),
        genes = lapply(panel$genes, function(g) list(i = .arr(g$i), v = .arr(g$v))),
        default_genes = .arr(pick$default)
      )
      rm(obj, panel); invisible(gc(FALSE))
    }
  }
  list(datasets = out, sample_counts = sample_counts, n_genes = n_genes, clusters = clusters)
}

# ---- Tables --------------------------------------------------------------------------

.TABLE_KINDS <- c("markers_all", "deg_sig", "cluster_markers", "deg_counts",
                  "top_expressed", "pseudobulk_deg", "pseudobulk_counts", "avg_expr",
                  "proportions", "prop_test", "gsea", "gsea_status", "escape",
                  "correlation", "other_table")

.cap_rows <- function(df, max_rows) {
  n <- nrow(df)
  if (n <= max_rows) return(df)
  lfc <- intersect(c("avg_log2FC", "logFC"), names(df))[1L]
  pcol <- intersect(c("p_val_adj", "padj", "FDR"), names(df))[1L]
  grp <- intersect(c("cluster", "group"), names(df))[1L]
  num <- vapply(df, is.numeric, logical(1))
  if (!is.na(lfc) && !is.na(pcol)) {
    score <- order(df[[pcol]], -abs(df[[lfc]]))
    df <- df[score, , drop = FALSE]
    if (!is.na(grp)) {
      quota <- ceiling(max_rows / max(1L, length(unique(df[[grp]]))))
      keep <- unlist(lapply(split(seq_len(nrow(df)), df[[grp]]), utils::head, quota),
                     use.names = FALSE)
      df <- df[sort(keep), , drop = FALSE]
    }
    return(utils::head(df, max_rows))
  }
  if (sum(num) >= 2L && !(!is.na(pcol))) {
    v <- apply(as.matrix(df[, num, drop = FALSE]), 1L, stats::var, na.rm = TRUE)
    return(df[order(-v)[seq_len(max_rows)], , drop = FALSE])
  }
  if (!is.na(pcol)) return(df[order(df[[pcol]])[seq_len(max_rows)], , drop = FALSE])
  utils::head(df, max_rows)
}

.table_json <- function(df) {
  df <- df[, !startsWith(names(df), "."), drop = FALSE]
  data <- lapply(df, function(col) {
    if (is.numeric(col)) .arr(signif(col, 5L)) else .arr(as.character(col))
  })
  list(cols = .arr(names(df)), data = data, n = nrow(df))
}

.report_tables <- function(bundle, max_rows, say) {
  cat <- bundle$catalog
  cat <- cat[cat$kind %in% .TABLE_KINDS & cat$ext %in% c("txt", "tsv", "csv"), , drop = FALSE]
  out <- list()
  # deg_sig is written one file per group: stack those; keep every other kind per file.
  key <- ifelse(cat$kind == "deg_sig",
                paste(cat$kind, cat$comparison, cat$unit, cat$grouping, sep = "|"), cat$rel)
  for (k in if (nrow(cat)) unique(key) else character()) {
    files <- cat[key == k, , drop = FALSE]
    frames <- lapply(files$path, .read_result_table)
    frames <- Filter(function(f) nrow(f) > 0L, frames)
    if (!length(frames)) next
    df <- if (length(frames) == 1L) frames[[1L]] else .rbind_fill(frames)
    n_total <- nrow(df)
    df <- .cap_rows(df, max_rows)
    kind <- files$kind[[1L]]
    stem <- if (kind == "deg_sig") "" else
      sub("\\.(txt|tsv|csv)$", "", sub("^[A-Za-z_]+?_(?=[A-Za-z0-9])", "", files$name[[1L]], perl = TRUE))
    where <- paste(stats::na.omit(c(files$comparison[[1L]], files$unit[[1L]], files$grouping[[1L]])),
                   collapse = " / ")
    label <- paste0(if (nzchar(where)) paste0(where, " - ") else "",
                    unname(.KIND_LABELS[kind]), if (nzchar(stem)) paste0(" [", stem, "]") else "")
    json <- c(.table_json(df), list(
      id = paste0("t", length(out) + 1L), kind = kind, label = label,
      comparison = files$comparison[[1L]], unit = files$unit[[1L]], grouping = files$grouping[[1L]],
      n_total = n_total, truncated = n_total > nrow(df)))
    out[[length(out) + 1L]] <- list(json = json, kind = kind, files = files, raw = df)
  }
  # Analyses run in the explorer (cached targets) are part of the results too.
  for (a in .interactive_tables(bundle$dir)) {
    df <- .cap_rows(a$table, max_rows)
    json <- c(.table_json(df), list(
      id = paste0("t", length(out) + 1L), kind = "explorer",
      label = paste0("Explorer - ", a$label), comparison = NA_character_, unit = NA_character_,
      grouping = NA_character_, n_total = nrow(a$table), truncated = nrow(a$table) > nrow(df)))
    out[[length(out) + 1L]] <- list(json = json, kind = "explorer", files = NULL, raw = df)
  }
  say("embedded ", length(out), " tables")
  out
}

# ---- Figures --------------------------------------------------------------------------

.report_figures <- function(bundle, sample_counts, say) {
  figs <- list()
  add <- function(tab, title, group, p) {
    spec <- tryCatch(.figure_spec_json(p), error = function(e) NULL)
    if (is.null(spec)) return(invisible())
    figs[[length(figs) + 1L]] <<- list(id = paste0("f", length(figs) + 1L), tab = tab,
                                       title = title, group = group, spec_json = spec)
  }
  if (!is.null(sample_counts) && nrow(sample_counts)) {
    add("overview", "Cells per sample", NULL, .safe_fig(sparks_fig_cells_per_sample(sample_counts)))
  }
  add("qc", "Cells before / after filtering", NULL, .safe_fig(sparks_fig_qc_stages(bundle$qc$stages)))
  add("qc", "Doublet calls", NULL, .safe_fig(sparks_fig_doublets(bundle$qc$doublets)))
  cat <- bundle$catalog
  groups <- unique(cat[cat$kind %in% c("deg_counts", "gsea", "proportions", "escape") &
                         !is.na(cat$grouping), c("comparison", "unit", "grouping")])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    where <- paste(g$comparison, g$unit, g$grouping, sep = " / ")
    rows <- cat$comparison %in% g$comparison & cat$unit %in% g$unit & cat$grouping %in% g$grouping
    dc <- sparks_table(bundle, "deg_counts", g$comparison, g$unit, g$grouping)
    if (nrow(dc)) add("de", "Significant genes per group", where, .safe_fig(sparks_fig_deg_counts(dc)))
    gs <- sparks_table(bundle, "gsea", g$comparison, g$unit, g$grouping)
    if (nrow(gs)) add("enrichment", "Gene-set enrichment", where, .safe_fig(sparks_fig_gsea(gs)))
    es <- sparks_table(bundle, "escape", g$comparison, g$unit, g$grouping)
    if (nrow(es) && all(c("pathway", "group", "mean_score") %in% names(es))) {
      mat <- tapply(es$mean_score, list(es$pathway, es$group), mean)
      add("enrichment", "Escape signature scores", where,
          .safe_fig(sparks_fig_heatmap(as.matrix(mat), diverging = FALSE)))
    }
    pr <- cat[rows & cat$kind == "proportions", , drop = FALSE]
    for (j in utils::head(seq_len(nrow(pr)), 4L)) {
      p <- .read_result_table(pr$path[[j]])
      by <- sub("\\.(txt|tsv)$", "", sub("^Proportions_", "", pr$name[[j]]))
      add("composition", paste0("Composition: ", by), where,
          .safe_fig(sparks_fig_composition(p)))
    }
  }
  say("built ", length(figs), " figures")
  figs
}

.safe_fig <- function(p) if (inherits(p, "plotly")) p else NULL

# ---- Images, files, run info ------------------------------------------------------------

.IMAGE_PRIORITY <- c("UMAP", "FeaturePlot", "Heatmap", "VlnPlot", "Expression", "Proportions",
                     "QC", "Correlation", "Pseudobulk", "DEG", "GSEA")

.report_images <- function(bundle, budget_mb, max_mb, say) {
  cat <- bundle$catalog
  img <- cat[cat$ext %in% c("png", "jpg", "jpeg", "svg", "gif", "webp"), , drop = FALSE]
  if (!nrow(img)) return(list())
  rank <- match(img$category, .IMAGE_PRIORITY)
  rank[is.na(rank)] <- length(.IMAGE_PRIORITY) + 1L
  img <- img[order(rank, img$size), , drop = FALSE]
  used <- 0
  mime <- c(png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg", svg = "image/svg+xml",
            gif = "image/gif", webp = "image/webp")
  res <- lapply(seq_len(nrow(img)), function(i) {
    mb <- img$size[[i]] / 1048576
    src <- NULL
    if (budget_mb > 0 && mb <= max_mb && used + mb <= budget_mb) {
      src <- paste0("data:", mime[[img$ext[[i]]]], ";base64,",
                    base64enc::base64encode(img$path[[i]]))
      used <<- used + mb
    }
    list(id = paste0("i", i), label = sub("\\.[A-Za-z]+$", "", img$name[[i]]),
         group = paste(stats::na.omit(c(img$unit[[i]], img$grouping[[i]])), collapse = " / "),
         category = img$category[[i]] %||% "", rel = img$rel[[i]], bytes = img$size[[i]],
         src = src)
  })
  say("images: ", sum(!vapply(res, function(r) is.null(r$src), logical(1))), " embedded of ", length(res))
  res
}

.report_files <- function(bundle) {
  cat <- bundle$catalog
  lapply(seq_len(nrow(cat)), function(i) list(rel = cat$rel[[i]], kind = cat$kind[[i]],
                                              size = cat$size[[i]], ext = cat$ext[[i]],
                                              category = cat$category[[i]] %||% ""))
}

.report_run <- function(bundle) {
  cfg_path <- file.path(bundle$dir, "config_used.yaml")
  log_tail <- if (isTRUE(bundle$run$has_log)) utils::tail(readLines(bundle$run$log, warn = FALSE), 200L) else character()
  list(
    config_yaml = if (file.exists(cfg_path)) paste(readLines(cfg_path, warn = FALSE), collapse = "\n") else "",
    warnings = .arr(utils::head(bundle$run$warnings, 200L)),
    session_info = .arr(utils::head(bundle$run$session_info, 120L)),
    log_tail = .arr(log_tail)
  )
}

.report_summary <- function(bundle, datasets, tables) {
  ds <- datasets$datasets
  cells <- if (length(ds)) max(vapply(ds, function(d) as.numeric(d$n_total), numeric(1))) else NA_real_
  samples <- if (!is.null(datasets$sample_counts)) nrow(datasets$sample_counts) else NA_integer_
  list(
    cells = if (is.na(cells)) NULL else cells,
    samples = if (is.na(samples)) NULL else samples,
    genes = if (is.na(datasets$n_genes)) NULL else datasets$n_genes,
    clusters = if (is.na(datasets$clusters)) NULL else datasets$clusters,
    tables = length(tables),
    files = nrow(bundle$catalog),
    comparisons = .arr(unique(stats::na.omit(bundle$units$comparison)))
  )
}

# ---- Serialisation & page ----------------------------------------------------------------

.report_json <- function(data, fig_specs) {
  json <- as.character(jsonlite::toJSON(data, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = NA, force = TRUE))
  for (i in seq_along(fig_specs)) {
    token <- sprintf("\"@@SPARKS_FIG_%d@@\"", i)
    pos <- regexpr(token, json, fixed = TRUE)
    if (pos > 0L) {
      json <- paste0(substr(json, 1L, pos - 1L), fig_specs[[i]],
                     substr(json, pos + nchar(token), nchar(json)))
    }
  }
  # Keep the JSON inert inside a <script> element.
  json <- gsub("</", "<\\/", json, fixed = TRUE)
  gsub("<!--", "<\\!--", json, fixed = TRUE)
}

.plotly_js_source <- function() {
  deps <- suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly())))
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) stop("Could not locate the plotly.js bundle.", call. = FALSE)
  dep <- dep[[1L]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) stop("Could not locate the plotly.js bundle.", call. = FALSE)
  js <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  gsub("</script", "<\\/script", js, fixed = TRUE)
}

.report_asset <- function(name) {
  path <- system.file("report", name, package = "SPARKS")
  if (!nzchar(path)) stop("Report asset not found: ", name, call. = FALSE)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.report_html <- function(meta, json) {
  esc <- htmltools::htmlEscape
  css <- .report_asset("sparks_report.css")
  js  <- .report_asset("sparks_report.js")

  # The inlined plotly bundle (~3.4 MB) exceeded the srcdoc size limit of the
  # outer <iframe> below; the attribute was truncated mid-script and left
  # window.Plotly undefined. Load Plotly from the jsDelivr CDN instead — the
  # embedded HTML stays small and Plotly runs in the outer frame's normal
  # script scope.
  plotly_cdn <- "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.25.2/plotly.min.js"

  inner <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", esc(meta$title), "</title>\n",
    "<script src=\"", plotly_cdn, "\"></script>\n",
    "<style>\n", css, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", esc(meta$title), "</h1>",
    if (isTRUE(meta$demo)) "<span class=\"badge demo\" title=\"Generated by sparks_demo_results()\">Synthetic demo data</span>" else "",
    "<span class=\"spacer\"></span>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Toggle colour theme\">Theme</button></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    # NOTE: .plotly_js_source() removed from <body> — Plotly now comes from <head>.
    "<script id=\"sparks-data\" type=\"application/json\">", json, "</script>\n",
    "<script>\n", js, "\n</script>\n</body>\n</html>\n"
  )

  # Outer shell: a full-viewport iframe with its own clean global scope.
  paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", esc(meta$title), "</title>\n</head>\n<body style=\"margin:0;padding:0;overflow:hidden;\">\n",
    "<iframe\n",
    "  srcdoc=\"", htmltools::htmlEscape(inner, attribute = TRUE), "\"\n",
    "  sandbox=\"allow-scripts\"\n",
    "  style=\"width:100vw;height:100vh;border:0;display:block;\"\n",
    "  title=\"SPARKS report\"\n",
    "></iframe>\n",
    "</body>\n</html>\n"
  )
}

#' Write the pipeline's consolidated report
#'
#' Used by the targets plan and by [sparks()]: writes the unified interactive
#' report, and falls back to the older static consolidated page only when the
#' packages the interactive one needs (plotly, jsonlite, htmltools, base64enc)
#' are missing.
#'
#' @param results_dir Pipeline results directory.
#' @param dependencies Ignored; lets a targets command declare upstream targets.
#' @param output_file Output path.
#' @return The written file path.
#' @keywords internal
.write_unified_report <- function(results_dir, dependencies = NULL,
                                  output_file = file.path(results_dir, "SPARKS_all_results.html")) {
  invisible(dependencies)
  miss <- .report_missing_packages()
  if (length(miss)) {
    message("   [report] Missing package(s) for the interactive report: ",
            paste(miss, collapse = ", "), ". Writing the static consolidated report instead.")
    return(.write_consolidated_report(results_dir, output_file = output_file))
  }
  args <- c(list(results_dir = results_dir, output_file = output_file, quiet = TRUE),
            .report_options(results_dir))
  do.call(sparks_report, args)
  normalizePath(output_file, winslash = "/", mustWork = TRUE)
}

#' Report options from the run's configuration
#'
#' Reads the optional \code{report:} section of \code{config_used.yaml} (written
#' next to the results by [sparks()]). Only known, correctly typed keys are
#' used; anything else is ignored so a typo never blocks the report.
#'
#' @param results_dir Pipeline results directory.
#' @return Named list of arguments for [sparks_report()].
#' @keywords internal
.report_options <- function(results_dir) {
  path <- file.path(results_dir, "config_used.yaml")
  if (!file.exists(path) || !requireNamespace("yaml", quietly = TRUE)) return(list())
  cfg <- tryCatch(yaml::read_yaml(path), error = function(e) NULL)
  rep <- if (is.list(cfg)) cfg$report else NULL
  if (!is.list(rep)) return(list())
  num <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  out <- list()
  for (k in c("max_cells", "max_genes", "max_table_rows")) {
    if (num(rep[[k]])) out[[k]] <- as.integer(rep[[k]])
  }
  for (k in c("image_budget_mb", "max_image_mb")) {
    if (num(rep[[k]])) out[[k]] <- as.numeric(rep[[k]])
  }
  if (is.logical(rep$embed_images) && length(rep$embed_images) == 1L && !is.na(rep$embed_images)) {
    out$embed_images <- rep$embed_images
  }
  if (is.character(rep$title) && length(rep$title) == 1L && nzchar(rep$title)) out$title <- rep$title
  out
}
