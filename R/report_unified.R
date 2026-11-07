#' Build the unified, dynamic SPARKS HTML report
#'
#' Produces **one self-contained HTML file** that gathers every result of a
#' pipeline run into tabs: an overview, QC (thresholds, input formats, cell
#' retention, per-sample QC summaries, doublets, cell cycle / sex, Alevin
#' barcode diagnostics), an interactive embedding explorer (lasso selection,
#' colour by metadata or gene), gene violin and dot plots, markers /
#' differential expression with a live volcano plot and per-group summaries,
#' enrichment (GSEA, escape, skip diagnostics), composition and group sizes,
#' every result table, the pipeline's figures, the file index and run
#' information (configuration, diagnostics, log, session info).
#'
#' The file needs no server, no network and no companion files: the page
#' styles and script are embedded in the package code, plotly.js is inlined
#' once from the installed \pkg{plotly} package, and the data are embedded as
#' JSON. Only a small index is parsed when the page opens; each table, dataset
#' and figure image sits in its own inert data block that is decoded the first
#' time it is shown, and each tab is drawn lazily when opened.
#'
#' Figures are drawn from compact "panels" extracted from the saved Seurat
#' objects ([sparks_panel()]): a down-sampled embedding, a handful of metadata
#' columns and the expression of a limited gene set. Counts and QC summaries
#' are computed on the full objects. Anything larger -- full differential
#' expression on an arbitrary selection, re-clustering, any gene -- is what the
#' Shiny explorer ([sparks_explorer()]) is for.
#'
#' @param results_dir Pipeline results directory.
#' @param output_file Output HTML path.
#' @param title Report title (default derived from the folder name).
#' @param max_cells Cells embedded per dataset (random down-sampling above).
#' @param max_genes Genes embedded per dataset for the violin / dot plots.
#' @param max_table_rows Row cap per embedded table (the files on disk are
#'   never truncated).
#' @param embed_images Embed the pipeline's PNG / SVG figures as base64.
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
  if (dir.exists(output_file)) {
    stop("output_file points to a directory: ", output_file, call. = FALSE)
  }
  say <- function(...) if (!isTRUE(quiet)) message("[sparks_report] ", ...)
  bundle <- bundle %||% sparks_read_results(results_dir)
  say("indexed ", nrow(bundle$catalog), " files")
  make_dir(dirname(output_file))
  out_dir <- normalizePath(dirname(output_file), winslash = "/", mustWork = TRUE)

  blobs <- new.env(parent = emptyenv())
  blobs$items <- list()
  add_blob <- function(id, json) {
    blobs$items[[length(blobs$items) + 1L]] <- list(id = id, json = json)
    invisible(id)
  }

  datasets <- .report_datasets(bundle, max_cells, max_genes, say)
  for (d in datasets$datasets) add_blob(paste0("sd-", d$head$id), .json(d$body))
  tables <- .report_tables(bundle, max_table_rows, say)
  for (t in tables) add_blob(paste0("sd-", t$json$id), .json(t$data))
  figures <- .report_figures(bundle, datasets$sample_counts, say)
  images <- .report_images(bundle, if (isTRUE(embed_images)) image_budget_mb else 0,
                           if (isTRUE(embed_images)) max_image_mb else 0, out_dir, say)

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
    settings = .report_settings(bundle$config),
    qc = .report_qc(bundle),
    notes = .report_notes(bundle, out_dir),
    datasets = lapply(datasets$datasets, `[[`, "head"),
    tables = lapply(tables, `[[`, "json"),
    figures = lapply(figures, function(f) {
      f$spec <- structure(f$spec_json, class = "json")
      f$spec_json <- NULL
      f
    }),
    images = lapply(images, `[[`, "head"),
    files = .report_files(bundle, out_dir),
    run = .report_run(bundle)
  )
  for (im in images) if (!is.null(im$src)) add_blob(paste0("si-", im$head$id), im$src)

  html <- .report_html(data$meta, .json(data), blobs$items)
  tmp <- paste0(output_file, ".part")
  con <- file(tmp, open = "wb")
  ok <- FALSE
  on.exit({
    if (!ok) unlink(tmp)
  }, add = TRUE)
  tryCatch(writeBin(charToRaw(enc2utf8(html)), con), finally = close(con))
  if (!file.rename(tmp, output_file)) {
    if (!file.copy(tmp, output_file, overwrite = TRUE)) {
      stop("Could not write the report to ", output_file, call. = FALSE)
    }
    unlink(tmp)
  }
  ok <- TRUE
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
  deg_counts = "DEG counts", deg_summary = "Marker / DEG summary",
  cluster_markers = "Positive cluster markers",
  top_expressed = "Top expressed genes", pseudobulk_deg = "Pseudobulk DEG (edgeR)",
  pseudobulk_counts = "Pseudobulk counts", avg_expr = "Average expression",
  proportions = "Cell proportions", prop_test = "Proportion tests",
  gsea = "GSEA", escape = "Escape signature scores",
  correlation = "Correlations", other_table = "Other tables",
  explorer = "Explorer analyses"
)

# JSON for an embedded block: jsonlite with the report's conventions, then made
# inert inside a <script> element. "<" only ever occurs inside JSON strings, where
# the \u003c escape is valid, so this cannot break parsing (unlike "<\!--").
.json <- function(x) {
  json <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", na = "null",
                                        digits = NA, force = TRUE, json_verbatim = TRUE))
  gsub("<", "\\u003c", json, fixed = TRUE)
}

.arr <- function(x) I(as.vector(x))

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

.report_gene_pick <- function(bundle, obj, comparison, max_genes) {
  cfg_genes <- unlist(bundle$config$genes$genes_to_plot, use.names = FALSE)
  sig_genes <- unlist(lapply(bundle$config$gene_signatures, `[[`, "genes"), use.names = FALSE)
  markers <- sparks_table(bundle, c("markers_all", "cluster_markers"), comparison = comparison)
  top <- character()
  if (nrow(markers) && all(c("gene", "cluster", "avg_log2FC") %in% names(markers))) {
    first <- markers[markers$.file == markers$.file[[1L]], , drop = FALSE]
    first$avg_log2FC <- suppressWarnings(as.numeric(first$avg_log2FC))
    top <- unlist(lapply(split(first, first$cluster), function(m) {
      utils::head(m$gene[order(-m$avg_log2FC)], 3L)
    }), use.names = FALSE)
  }
  vf <- tryCatch(SeuratObject::VariableFeatures(obj), error = function(e) character())
  genes <- unique(c(cfg_genes, top, sig_genes, utils::head(vf, max_genes)))
  genes <- genes[genes %in% rownames(obj)]
  list(genes = utils::head(genes, max_genes),
       default = utils::head(unique(c(intersect(cfg_genes, genes), top)), 12L))
}

.report_datasets <- function(bundle, max_cells, max_genes, say) {
  objs <- .discover_result_objects(bundle)
  out <- list(); sample_counts <- list(); n_genes <- NA_integer_; clusters <- NA_integer_
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
    meta <- obj[[]]
    if (objs$kind[[i]] == "merged") {
      if ("orig.ident" %in% colnames(meta)) {
        sample_counts[[length(sample_counts) + 1L]] <- as.data.frame(
          table(sample = as.character(meta$orig.ident)), responseName = "cells",
          stringsAsFactors = FALSE)
      }
      if (is.na(n_genes)) n_genes <- nrow(obj)
      if ("seurat_clusters" %in% colnames(meta)) {
        clusters <- sum(clusters, length(unique(meta$seurat_clusters)), na.rm = TRUE)
      }
    }
    id <- paste0("ds", length(out) + 1L)
    head <- list(
      id = id, label = objs$label[[i]], kind = objs$kind[[i]],
      comparison = objs$comparison[[i]], n_total = panel$n_total,
      n_cells = length(panel$cells), reduction = panel$reduction,
      n_genes_embedded = length(panel$genes),
      facts = .report_facts(obj, panel)
    )
    body <- list(
      cells = .arr(panel$cells), x = .arr(panel$x), y = .arr(panel$y),
      meta = lapply(panel$meta, function(m) {
        if (identical(m$type, "cat")) list(type = "cat", levels = .arr(m$levels), codes = .arr(m$codes))
        else list(type = "num", values = .arr(m$values))
      }),
      genes = lapply(panel$genes, function(g) list(i = .arr(g$i), v = .arr(g$v))),
      default_genes = .arr(pick$default)
    )
    out[[length(out) + 1L]] <- list(head = head, body = body)
    rm(obj, panel, meta); invisible(gc(FALSE))
  }
  sc <- if (length(sample_counts)) {
    s <- do.call(rbind, sample_counts)
    stats::aggregate(cells ~ sample, data = s, FUN = sum)
  } else NULL
  list(datasets = out, sample_counts = sc, n_genes = n_genes, clusters = clusters)
}

# Facts computed on the *full* object (the panel is down-sampled): sizes, the
# per-sample table (condition, input format), per-sample QC summaries and exact
# level counts of every categorical column shown in the report.
.report_facts <- function(obj, panel) {
  meta <- obj[[]]
  has <- function(col) col %in% colnames(meta)
  samples <- NULL
  if (has("orig.ident")) {
    s <- as.character(meta$orig.ident)
    first <- !duplicated(s)
    samples <- list(
      sample = .arr(s[first]),
      cells = .arr(as.integer(table(factor(s, levels = s[first])))),
      condition = if (has("condition")) .arr(as.character(meta$condition[first])) else NULL,
      quant_format = if (has("quant_format")) .arr(as.character(meta$quant_format[first])) else NULL
    )
  }
  qc_cols <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"), colnames(meta))
  grp_col <- if (has("orig.ident")) "orig.ident" else if (has("condition")) "condition" else NULL
  qc <- NULL
  if (length(qc_cols) && !is.null(grp_col)) {
    g <- as.character(meta[[grp_col]])
    rows <- lapply(qc_cols, function(m) {
      v <- suppressWarnings(as.numeric(meta[[m]]))
      do.call(rbind, lapply(split(v, factor(g, levels = unique(g))), function(x) {
        x <- x[is.finite(x)]
        q <- if (length(x)) stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE) else rep(NA_real_, 3L)
        data.frame(cells = length(x), median = q[[2L]], q1 = q[[1L]], q3 = q[[3L]])
      }))
    })
    groups <- unique(g)
    qc_df <- data.frame(
      group = rep(groups, times = length(qc_cols)),
      metric = rep(qc_cols, each = length(groups)),
      do.call(rbind, rows), stringsAsFactors = FALSE, row.names = NULL)
    qc_df[c("median", "q1", "q3")] <- lapply(qc_df[c("median", "q1", "q3")], round, 2L)
    names(qc_df)[1L] <- grp_col
    qc <- .table_json(qc_df)
  }
  cat_cols <- names(panel$meta)[vapply(panel$meta, function(m) identical(m$type, "cat"), logical(1))]
  extra <- intersect(c("Phase", "Assigned_Sex", "scDblFinder.class", "condition"), colnames(meta))
  counts <- lapply(stats::setNames(nm = unique(c(cat_cols, extra))), function(col) {
    v <- meta[[col]]
    f <- if (is.factor(v)) droplevels(v) else factor(as.character(v))
    tab <- table(f, useNA = "ifany")
    lv <- names(tab)
    lv[is.na(lv)] <- "NA"
    list(levels = .arr(lv), counts = .arr(as.integer(tab)))
  })
  list(
    cells = ncol(obj), genes = nrow(obj),
    samples_n = if (has("orig.ident")) length(unique(meta$orig.ident)) else NULL,
    conditions_n = if (has("condition")) length(unique(meta$condition)) else NULL,
    clusters_n = if (has("seurat_clusters")) length(unique(meta$seurat_clusters)) else NULL,
    singler = .arr(grep("^singleR_labels", colnames(meta), value = TRUE)),
    escape_assays = .arr(grep("^escape\\.", names(obj@assays), value = TRUE)),
    samples = samples, qc = qc, counts = counts
  )
}

# ---- Tables --------------------------------------------------------------------------

.TABLE_KINDS <- c("markers_all", "deg_sig", "cluster_markers", "deg_counts",
                  "top_expressed", "pseudobulk_deg", "pseudobulk_counts", "avg_expr",
                  "proportions", "prop_test", "gsea", "escape",
                  "correlation", "other_table")

# Kinds whose per-group files are stacked into one table per analysis.
.STACKED_KINDS <- c("deg_sig", "gsea", "escape")

.num_col <- function(df, cands) {
  hit <- intersect(cands, names(df))[1L]
  if (is.na(hit)) return(NULL)
  suppressWarnings(as.numeric(df[[hit]]))
}

.cap_rows <- function(df, max_rows) {
  n <- nrow(df)
  if (n <= max_rows) return(df)
  lfc <- .num_col(df, c("avg_log2FC", "logFC", "log2FoldChange"))
  pv <- .num_col(df, c("p_val_adj", "padj", "FDR", "adj.P.Val"))
  grp <- intersect(c("cluster", "group", "Group"), names(df))[1L]
  if (!is.null(lfc) && !is.null(pv)) {
    df <- df[order(pv, -abs(lfc), na.last = TRUE), , drop = FALSE]
    if (!is.na(grp)) {
      quota <- ceiling(max_rows / max(1L, length(unique(df[[grp]]))))
      keep <- unlist(lapply(split(seq_len(nrow(df)), df[[grp]]), utils::head, quota),
                     use.names = FALSE)
      df <- df[sort(keep), , drop = FALSE]
    }
    return(utils::head(df, max_rows))
  }
  if (!is.null(pv)) return(df[order(pv, na.last = TRUE)[seq_len(max_rows)], , drop = FALSE])
  num <- vapply(df, is.numeric, logical(1))
  if (sum(num) >= 2L) {
    v <- apply(as.matrix(df[, num, drop = FALSE]), 1L, stats::var, na.rm = TRUE)
    return(df[order(-v, na.last = TRUE)[seq_len(max_rows)], , drop = FALSE])
  }
  utils::head(df, max_rows)
}

.table_json <- function(df) {
  df <- df[, !startsWith(names(df), "."), drop = FALSE]
  data <- lapply(df, function(col) {
    if (is.numeric(col)) .arr(signif(col, 5L)) else .arr(as.character(col))
  })
  list(cols = .arr(names(df)), data = data, n = nrow(df))
}

# Per-group summary of a full marker / DEG table (computed before any row cap):
# genes tested, genes with adjusted p < 0.05, median |log2FC|.
.deg_summary <- function(df) {
  grp <- intersect(c("cluster", "group", "Group", "contrast"), names(df))[1L]
  pv <- .num_col(df, c("p_val_adj", "padj", "FDR", "adj.P.Val"))
  lfc <- .num_col(df, c("avg_log2FC", "avg_logFC", "logFC", "log2FoldChange"))
  if (is.null(pv) || is.null(lfc)) return(NULL)
  g <- if (is.na(grp)) rep("All", nrow(df)) else as.character(df[[grp]])
  idx <- split(seq_len(nrow(df)), factor(g, levels = unique(g)))
  data.frame(
    group = names(idx),
    tested = vapply(idx, length, integer(1)),
    padj_below_0.05 = vapply(idx, function(i) sum(!is.na(pv[i]) & pv[i] < 0.05), integer(1)),
    up = vapply(idx, function(i) sum(!is.na(pv[i]) & pv[i] < 0.05 & lfc[i] > 0, na.rm = TRUE), integer(1)),
    down = vapply(idx, function(i) sum(!is.na(pv[i]) & pv[i] < 0.05 & lfc[i] < 0, na.rm = TRUE), integer(1)),
    median_abs_log2FC = vapply(idx, function(i) {
      x <- abs(lfc[i]); x <- x[is.finite(x)]
      if (length(x)) round(stats::median(x), 3L) else NA_real_
    }, numeric(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

.strip_kind_prefix <- function(name) {
  stem <- sub("\\.(txt|tsv|csv)$", "", name)
  for (p in names(.RESULT_KIND_PREFIXES)) {
    if (startsWith(stem, p)) return(substring(stem, nchar(p) + 1L))
  }
  stem
}

.report_tables <- function(bundle, max_rows, say) {
  cat <- bundle$catalog
  cat <- cat[cat$kind %in% .TABLE_KINDS & cat$ext %in% c("txt", "tsv", "csv"), , drop = FALSE]
  out <- list()
  push <- function(df, kind, label, files, n_total, extra = list()) {
    capped <- .cap_rows(df, max_rows)
    tj <- .table_json(capped)
    id <- paste0("t", length(out) + 1L)
    head <- c(list(
      id = id, kind = kind, label = label, cols = tj$cols, n = tj$n,
      comparison = files$comparison[[1L]] %||% NA_character_,
      unit = files$unit[[1L]] %||% NA_character_,
      grouping = files$grouping[[1L]] %||% NA_character_,
      n_total = n_total, truncated = n_total > nrow(capped)), extra)
    out[[length(out) + 1L]] <<- list(json = head, data = tj$data)
    id
  }
  method <- ifelse(cat$kind == "escape", sub("^Escape_Scores_([^_]+)_.*$", "\\1", cat$name), "")
  key <- ifelse(cat$kind %in% .STACKED_KINDS,
                paste(cat$kind, cat$comparison, cat$unit, cat$grouping, method, sep = "|"),
                cat$rel)
  for (k in if (nrow(cat)) unique(key) else character()) {
    files <- cat[key == k, , drop = FALSE]
    frames <- lapply(files$path, .read_result_table)
    frames <- Filter(function(f) nrow(f) > 0L, frames)
    if (!length(frames)) next
    df <- if (length(frames) == 1L) frames[[1L]] else .rbind_fill(frames)
    kind <- files$kind[[1L]]
    stem <- if (kind == "escape") sub("^Escape_Scores_([^_]+)_.*$", "\\1", files$name[[1L]]) else
      if (kind %in% .STACKED_KINDS && nrow(files) > 1L) "" else .strip_kind_prefix(files$name[[1L]])
    where <- paste(stats::na.omit(c(files$comparison[[1L]], files$unit[[1L]], files$grouping[[1L]])),
                   collapse = " / ")
    label <- paste0(if (nzchar(where)) paste0(where, " - ") else "",
                    unname(.KIND_LABELS[kind]), if (nzchar(stem)) paste0(" [", stem, "]") else "")
    id <- push(df, kind, label, files, nrow(df))
    if (kind %in% c("markers_all", "deg_sig", "cluster_markers", "pseudobulk_deg")) {
      sm <- .deg_summary(df)
      if (!is.null(sm) && nrow(sm)) {
        push(sm, "deg_summary", paste0(label, " - summary"), files, nrow(sm),
             list(summary_of = id))
      }
    }
  }
  # Analyses run in the explorer (cached targets) are part of the results too.
  for (a in .interactive_tables(bundle$dir)) {
    files <- data.frame(comparison = NA_character_, unit = NA_character_, grouping = NA_character_)
    push(a$table, "explorer", paste0("Explorer - ", a$label), files, nrow(a$table))
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
  alevin <- tryCatch(.report_alevin(bundle), error = function(e) NULL)
  for (a in alevin) add("qc", a$title, "Alevin barcode diagnostics", .safe_fig(a$fig))

  cat <- bundle$catalog
  groups <- unique(cat[cat$kind %in% c("deg_counts", "gsea", "proportions", "escape") &
                         !is.na(cat$grouping), c("comparison", "unit", "grouping")])
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    where <- paste(stats::na.omit(c(g$comparison, g$unit, g$grouping)), collapse = " / ")
    rows <- cat$comparison %in% g$comparison & cat$unit %in% g$unit & cat$grouping %in% g$grouping
    dc <- sparks_table(bundle, "deg_counts", g$comparison, g$unit, g$grouping)
    if (nrow(dc)) add("de", "Significant genes per group", where, .safe_fig(sparks_fig_deg_counts(dc)))
    gs <- sparks_table(bundle, "gsea", g$comparison, g$unit, g$grouping)
    if (nrow(gs)) add("enrichment", "Gene-set enrichment", where, .safe_fig(sparks_fig_gsea(gs)))
    es <- sparks_table(bundle, "escape", g$comparison, g$unit, g$grouping)
    if (nrow(es) && all(c("pathway", "group", "mean_score") %in% names(es))) {
      # One heatmap per scoring method: UCell and AUCell scores are not comparable.
      method <- if ("method" %in% names(es)) as.character(es$method) else
        sub("^.*Escape_Scores_([^_]+)_.*$", "\\1", es$.file)
      for (m in unique(method)) {
        e <- es[method == m, , drop = FALSE]
        mat <- tapply(suppressWarnings(as.numeric(e$mean_score)), list(e$pathway, e$group), mean)
        add("enrichment", paste0("Escape signature scores (", m, ")"), where,
            .safe_fig(sparks_fig_heatmap(as.matrix(mat), diverging = FALSE)))
      }
    }
    pr <- cat[rows & cat$kind == "proportions", , drop = FALSE]
    for (j in utils::head(seq_len(nrow(pr)), 4L)) {
      p <- .read_result_table(pr$path[[j]])
      add("composition", paste0("Composition: ", .strip_kind_prefix(pr$name[[j]])), where,
          .safe_fig(sparks_fig_composition(p)))
    }
  }
  say("built ", length(figs), " figures")
  figs
}

.safe_fig <- function(p) if (inherits(p, "plotly")) p else NULL

# ---- Alevin barcode diagnostics ---------------------------------------------------------

# Best-effort: needs the run's config (config_used.yaml) and the original Alevin
# folders. Curves are thinned on a log scale and scatters are capped per sample,
# so a run with millions of raw barcodes still yields a small report.
.report_alevin <- function(bundle, max_points = 4000L) {
  cfg <- bundle$config
  if (is.null(cfg)) return(list())
  st <- tryCatch(load_sample_table(cfg), error = function(e) NULL)
  if (is.null(st) || !nrow(st)) return(list())
  data_dir <- cfg$pipeline$data_dir %||% ""
  roots <- unique(c(data_dir,
                    if (!is.null(cfg$pipeline$config_dir)) file.path(cfg$pipeline$config_dir, data_dir),
                    file.path(dirname(bundle$dir), data_dir)))
  roots <- roots[nzchar(roots) & dir.exists(roots)]
  if (!length(roots)) return(list())
  read <- lapply(seq_len(nrow(st)), function(i) {
    folder <- st$folder_id[[i]]
    sample_dir <- file.path(roots, folder)
    sample_dir <- sample_dir[dir.exists(sample_dir)][1L]
    if (is.na(sample_dir)) return(NULL)
    fmt <- if ("quant_format" %in% names(st) && !is.na(st$quant_format[[i]]) &&
               nzchar(st$quant_format[[i]])) st$quant_format[[i]] else cfg$input$format %||% "auto"
    if (identical(fmt, "auto")) {
      fmt <- tryCatch(.detect_quant_format(sample_dir), error = function(e) NA_character_)
    }
    if (!identical(fmt, "alevin")) return(NULL)
    ad <- file.path(sample_dir, "alevin")
    freq_f <- file.path(ad, "raw_cb_frequency.txt")
    feat_f <- file.path(ad, "featureDump.txt")
    if (!file.exists(freq_f) || !file.exists(feat_f)) return(NULL)
    freq <- tryCatch(utils::read.delim(freq_f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                                       col.names = c("CB", "originalFreq"), comment.char = ""),
                     error = function(e) NULL)
    feat <- tryCatch(utils::read.delim(feat_f, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                       check.names = FALSE, comment.char = ""),
                     error = function(e) NULL)
    if (is.null(freq) || is.null(feat) || !"CB" %in% names(feat)) return(NULL)
    aliases <- c(MappingRate = "mappingRate", CorrectedReads = "collapsedFreq",
                 DedupRate = "dedupRate", NumGenesOverMean = "nbrGenesAboveMean",
                 MappedReads = "nbrMappedUMI", DeduplicatedReads = "totalUMICount",
                 NumGenesExpressed = "nbrGenesAboveZero")
    hit <- intersect(names(aliases), names(feat))
    names(feat)[match(hit, names(feat))] <- unname(aliases[hit])
    keep <- intersect(c("CB", "collapsedFreq", "totalUMICount", "nbrGenesAboveZero"), names(feat))
    wl_f <- file.path(ad, "whitelist.txt")
    wl <- if (file.exists(wl_f)) utils::read.delim(wl_f, header = FALSE, stringsAsFactors = FALSE)[[1L]] else character()
    freq$originalFreq <- suppressWarnings(as.numeric(freq$originalFreq))
    freq <- freq[is.finite(freq$originalFreq) & freq$originalFreq > 0, , drop = FALSE]
    freq <- freq[order(-freq$originalFreq), , drop = FALSE]
    freq$rank <- seq_len(nrow(freq))
    m <- merge(freq, feat[, keep, drop = FALSE], by = "CB")
    list(sample = paste(st$protocol[[i]], folder, sep = "_"), freq = freq, merged = m,
         whitelist = wl)
  })
  read <- Filter(Negate(is.null), read)
  if (!length(read)) return(list())
  .require_plotly()
  cols <- .level_colors(vapply(read, `[[`, character(1), "sample"))
  # Knee plot: log-spaced ranks.
  knee <- plotly::plot_ly()
  for (r in read) {
    n <- nrow(r$freq)
    if (!n) next
    idx <- unique(round(exp(seq(0, log(n), length.out = min(n, 1500L)))))
    knee <- plotly::add_trace(knee, type = "scatter", mode = "lines", x = r$freq$rank[idx],
                              y = r$freq$originalFreq[idx], name = r$sample,
                              line = list(color = cols[[r$sample]], width = 2),
                              hovertemplate = paste0(r$sample, "<br>rank %{x}<br>%{y} reads<extra></extra>"))
    if (length(r$whitelist)) {
      cut <- sum(r$freq$CB %in% r$whitelist)
      if (cut > 0) {
        knee <- plotly::add_trace(knee, type = "scatter", mode = "markers", x = cut,
                                  y = r$freq$originalFreq[[min(cut, n)]], showlegend = FALSE,
                                  marker = list(color = cols[[r$sample]], size = 9, symbol = "line-ns-open"),
                                  hovertemplate = paste0(r$sample, ": ", cut, " whitelisted barcodes<extra></extra>"))
      }
    }
  }
  knee <- plotly::layout(knee, xaxis = list(type = "log"), yaxis = list(type = "log"))
  knee <- .fig_config(.fig_base(knee, NULL, "Barcode rank", "Raw read frequency", TRUE), "alevin_knee")
  scatter <- function(xcol, ycol, xlab, ylab, file) {
    p <- plotly::plot_ly()
    any_pts <- FALSE
    for (r in read) {
      m <- r$merged
      if (!all(c(xcol, ycol) %in% names(m))) next
      x <- suppressWarnings(as.numeric(m[[xcol]])); y <- suppressWarnings(as.numeric(m[[ycol]]))
      ok <- which(is.finite(x) & x > 0 & is.finite(y) & y > 0)
      if (!length(ok)) next
      if (length(ok) > max_points) {
        old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
        set.seed(1L); ok <- sort(sample(ok, max_points))
        if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
      }
      inwl <- m$CB[ok] %in% r$whitelist
      for (w in c(TRUE, FALSE)) {
        k <- ok[inwl == w]
        if (!length(k)) next
        p <- plotly::add_trace(p, type = "scattergl", mode = "markers", x = x[k], y = y[k],
                               name = paste0(r$sample, if (w) " (whitelist)" else " (other)"),
                               legendgroup = r$sample,
                               marker = list(color = cols[[r$sample]], size = 4,
                                             opacity = if (w) 0.8 else 0.25),
                               hovertemplate = paste0(r$sample, "<br>", xlab, " %{x}<br>", ylab,
                                                      " %{y}<extra></extra>"))
        any_pts <- TRUE
      }
    }
    if (!any_pts) return(NULL)
    p <- plotly::layout(p, xaxis = list(type = "log"), yaxis = list(type = "log"))
    .fig_config(.fig_base(p, NULL, xlab, ylab, TRUE), file)
  }
  list(
    list(title = "Raw barcode frequency (knee plot)", fig = knee),
    list(title = "Barcode frequency vs UMIs",
         fig = scatter("collapsedFreq", "totalUMICount", "Collapsed frequency", "Total UMIs", "alevin_umi")),
    list(title = "UMIs vs detected genes",
         fig = scatter("totalUMICount", "nbrGenesAboveZero", "Total UMIs", "Detected genes", "alevin_genes"))
  )
}

# ---- Images, files, run info ------------------------------------------------------------

.IMAGE_PRIORITY <- c("UMAP", "FeaturePlot", "Heatmap", "VlnPlot", "Expression", "Proportions",
                     "QC", "Correlation", "Pseudobulk", "DEG", "GSEA", "Escape")

.href <- function(path, out_dir) {
  tryCatch(.report_file_href(path, out_dir), error = function(e) NA_character_)
}

.report_images <- function(bundle, budget_mb, max_mb, out_dir, say) {
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
    list(head = list(
      id = paste0("i", i), label = sub("\\.[A-Za-z]+$", "", img$name[[i]]),
      group = paste(stats::na.omit(c(img$comparison[[i]], img$unit[[i]], img$grouping[[i]])),
                    collapse = " / "),
      category = img$category[[i]] %||% "", rel = img$rel[[i]],
      href = .href(img$path[[i]], out_dir), bytes = img$size[[i]],
      embedded = !is.null(src)), src = src)
  })
  say("images: ", sum(vapply(res, function(r) isTRUE(r$head$embedded), logical(1))),
      " embedded of ", length(res))
  res
}

.report_files <- function(bundle, out_dir) {
  cat <- bundle$catalog
  if (!nrow(cat)) return(list())
  list(rel = .arr(cat$rel), kind = .arr(cat$kind), size = .arr(cat$size),
       category = .arr(ifelse(is.na(cat$category), "", cat$category)),
       href = .arr(vapply(cat$path, .href, character(1), out_dir = out_dir, USE.NAMES = FALSE)))
}

.report_run <- function(bundle) {
  cfg_path <- file.path(bundle$dir, "config_used.yaml")
  log_lines <- if (isTRUE(bundle$run$has_log)) readLines(bundle$run$log, warn = FALSE) else character()
  diag <- grep("\\[(WARNING|ERROR)\\]", log_lines, value = TRUE)
  skips <- grep("\\[SKIP", log_lines, value = TRUE)
  list(
    config_yaml = if (file.exists(cfg_path)) paste(readLines(cfg_path, warn = FALSE), collapse = "\n") else "",
    warnings = .arr(utils::head(bundle$run$warnings, 500L)),
    diagnostics = list(
      type = .arr(c(rep("Warning / error", length(diag)), rep("Skipped / unavailable", length(skips)))),
      message = .arr(c(diag, skips))),
    session_info = .arr(utils::head(bundle$run$session_info, 200L)),
    log_tail = .arr(utils::tail(log_lines, 300L)),
    has_log = isTRUE(bundle$run$has_log)
  )
}

# Settings that explain the results: QC thresholds, input handling and which
# optional analyses were switched on (so empty sections can say why).
.report_settings <- function(cfg) {
  if (!is.list(cfg)) return(NULL)
  val <- function(x) if (is.null(x) || !length(x)) NULL else if (length(x) == 1L) x else .arr(unlist(x))
  list(
    qc = list(min_features = val(cfg$qc$min_features), max_features = val(cfg$qc$max_features),
              max_counts = val(cfg$qc$max_counts), max_mt_percent = val(cfg$qc$max_mt_percent),
              min_cells = val(cfg$qc$min_cells)),
    gene_removal_pattern = val(cfg$species$gene_removal_pattern),
    mt_pattern = val(cfg$species$mt_pattern),
    input = list(format = val(cfg$input$format), save_raw_matrix = val(cfg$input$save_raw_matrix),
                 raw_matrix_format = val(cfg$input$raw_matrix_format)),
    run = list(deg = val(cfg$deg$run), pseudobulk = val(cfg$pseudobulk$run),
               gsea = val(cfg$gsea$run), escape = val(cfg$escape$run),
               proportions = val(cfg$proportions$run), sex_scoring = val(cfg$sex_scoring$run),
               cell_cycle = val(cfg$cell_cycle$run),
               correlation = !is.null(cfg$genes$corr_genes_x) && !is.null(cfg$genes$corr_genes_y)),
    gsea_collection = val(cfg$gsea$collection)
  )
}

# Per-sample cell retention (before / after QC filtering and doublet removal).
.report_qc <- function(bundle) {
  st <- bundle$qc$stages
  out <- list()
  if (is.data.frame(st) && nrow(st)) {
    st$sample <- ifelse(is.na(st$sample), st$label, st$sample)
    b <- st[st$stage == "before_filtering", c("sample", "cells"), drop = FALSE]
    a <- st[st$stage == "after_filtering", c("sample", "cells"), drop = FALSE]
    if (nrow(b) || nrow(a)) {
      s <- unique(c(b$sample, a$sample))
      before <- b$cells[match(s, b$sample)]
      after <- a$cells[match(s, a$sample)]
      out$retention <- .table_json(data.frame(
        sample = s, before = before, after = after,
        retained_pct = ifelse(!is.na(before) & before > 0, round(100 * after / before, 1), NA_real_),
        stringsAsFactors = FALSE))
    }
  }
  d <- bundle$qc$doublets
  if (is.data.frame(d) && nrow(d)) {
    out$doublets <- .table_json(d[, intersect(c("sample", "class", "cells"), names(d)), drop = FALSE])
  }
  out
}

# Free-text diagnostics written next to the results (GSEA skip reasons).
.report_notes <- function(bundle, out_dir) {
  cat <- bundle$catalog
  st <- cat[cat$kind == "gsea_status", , drop = FALSE]
  if (!nrow(st)) return(list())
  list(gsea_status = list(
    file = .arr(st$rel),
    text = .arr(vapply(st$path, function(p) paste(readLines(p, warn = FALSE), collapse = " "),
                       character(1), USE.NAMES = FALSE))))
}

.report_summary <- function(bundle, datasets, tables) {
  merged <- Filter(function(d) identical(d$head$kind, "merged"), datasets$datasets)
  if (!length(merged)) merged <- datasets$datasets
  cells <- if (length(merged)) sum(vapply(merged, function(d) as.numeric(d$head$n_total), numeric(1))) else NA_real_
  samples <- if (!is.null(datasets$sample_counts)) nrow(datasets$sample_counts) else NA_integer_
  conds <- unique(unlist(lapply(merged, function(d) d$head$facts$samples$condition), use.names = FALSE))
  list(
    cells = if (is.na(cells)) NULL else cells,
    samples = if (is.na(samples)) NULL else samples,
    conditions = if (length(conds)) length(conds) else NULL,
    genes = if (is.na(datasets$n_genes)) NULL else datasets$n_genes,
    clusters = if (is.na(datasets$clusters)) NULL else datasets$clusters,
    tables = length(tables),
    files = nrow(bundle$catalog),
    comparisons = .arr(unique(stats::na.omit(bundle$units$comparison[
      bundle$units$comparison != "Sample_QC"])))
  )
}

# ---- Page ---------------------------------------------------------------------------------

.plotly_js_source <- function() {
  deps <- suppressWarnings(suppressMessages(htmltools::findDependencies(plotly::plot_ly())))
  dep <- Filter(function(d) identical(d$name, "plotly-main"), deps)
  if (!length(dep)) stop("Could not locate the plotly.js bundle.", call. = FALSE)
  dep <- dep[[1L]]
  path <- file.path(dep$src$file, dep$script)
  if (!file.exists(path)) path <- system.file(dep$src$file, dep$script, package = dep$package %||% "plotly")
  if (!nzchar(path) || !file.exists(path)) stop("Could not locate the plotly.js bundle.", call. = FALSE)
  js <- rawToChar(readBin(path, "raw", n = file.info(path)$size))
  Encoding(js) <- "UTF-8"
  # Keep the bundle inert for the HTML parser: no "</script" and no "<!--".
  js <- gsub("</script", "<\\/script", js, fixed = TRUE, useBytes = TRUE)
  gsub("<!--", "<\\!--", js, fixed = TRUE, useBytes = TRUE)
}

.report_html <- function(meta, json, blobs) {
  esc <- htmltools::htmlEscape
  blob_html <- vapply(blobs, function(b) {
    type <- if (startsWith(b$id, "si-")) "text/plain" else "application/json"
    paste0("<script type=\"", type, "\" id=\"", b$id, "\">", b$json, "</script>")
  }, character(1))
  paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<meta name=\"generator\" content=\"SPARKS ", esc(meta$version), "\">\n",
    "<title>", esc(meta$title), "</title>\n",
    "<style>\n", .SPARKS_REPORT_CSS, "\n</style>\n</head>\n<body>\n",
    "<header class=\"top\"><div class=\"top-inner\">\n",
    "<div class=\"title-row\"><h1>", esc(meta$title), "</h1>",
    if (isTRUE(meta$demo)) "<span class=\"badge demo\" title=\"Generated by sparks_demo_results()\">Synthetic demo data</span>" else "",
    "<span class=\"spacer\"></span>",
    "<button id=\"theme-btn\" class=\"btn\" type=\"button\" aria-label=\"Toggle colour theme\">Theme</button></div>\n",
    "<nav id=\"tabs\" class=\"tabs\" role=\"tablist\" aria-label=\"Report sections\"></nav>\n",
    "</div></header>\n<main id=\"panels\"><p class=\"empty\" id=\"boot\">Loading report\u2026</p></main>\n",
    "<noscript><p style=\"padding:16px\">This report needs JavaScript enabled.</p></noscript>\n",
    "<script id=\"sparks-data\" type=\"application/json\">", json, "</script>\n",
    paste(blob_html, collapse = "\n"), "\n",
    "<script>\n", .plotly_js_source(), "\n</script>\n",
    "<script>\n", .SPARKS_REPORT_JS, "\n</script>\n</body>\n</html>\n"
  )
}

#' Write the pipeline's consolidated report
#'
#' Used by the targets plan and by [sparks()]: writes the unified interactive
#' report, and falls back to the static consolidated page only when the
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
