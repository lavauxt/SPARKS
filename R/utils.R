#' Null coalescing operator
#'
#' Returns \code{rhs} when \code{lhs} is \code{NULL} \emph{or has length zero}
#' (\code{character(0)}, \code{list()}, ... -- an empty YAML sequence \code{[]}
#' parses to \code{list()}). This is deliberately broader than \code{\%||\%} in
#' base R (>= 4.4.0), rlang and SeuratObject, which only test for \code{NULL}:
#' the pipeline relies on it so that an empty list in a config file falls back
#' to the default. Keep the difference in mind if you attach SPARKS together
#' with one of those packages.
#' @name %||%
#' @param lhs Left hand side
#' @param rhs Right hand side
#' @export
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
}

.validate_seed <- function(seed) {
  if (length(seed) != 1L || is.na(seed) || !is.numeric(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be a single integer between 0 and .Machine$integer.max.")
  }
  as.integer(seed)
}

.set_seed_preserving_state <- function(seed) {
  seed <- .validate_seed(seed)

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  old_kind <- RNGkind()

  set.seed(seed)

  function() {
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}

#' Create directory recursively if it doesn't exist
#' @param path Character. Directory path
#' @return Invisible path
#' @export
make_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

#' Safely execute an expression, returning a fallback on error
#' @param expr Expression to evaluate
#' @param label Character. Label to print if error occurs
#' @param fallback Value to return on error
#' @return Result of expr or fallback
#' @export
safe_run <- function(expr, label = "Task", fallback = NULL) {
  tryCatch({
    expr
  }, error = function(e) {
    message("   [ERROR] ", label, " failed: ", e$message)
    fallback
  })
}

#' Sanitize string for filenames
#' @param x Character vector
#' @return Sanitized character vector
.safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

.find_report_template <- function(configured, config_dir, filename) {
  if (!is.null(configured)) return(configured)

  namespace_path <- tryCatch(
    getNamespaceInfo(environment(.find_report_template), "path"),
    error = function(e) ""
  )
  candidates <- c(
    file.path(config_dir, filename),
    file.path(getwd(), filename),
    file.path(getwd(), "inst", "rmd", filename),
    file.path(namespace_path, "rmd", filename),
    file.path(namespace_path, "inst", "rmd", filename)
  )
  candidates <- candidates[nzchar(candidates)]
  found <- candidates[file.exists(candidates)]
  if (length(found) > 0L) found[[1L]] else NULL
}

.with_seurat_aggregate_advisory_muffled <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("As of Seurat v5, we recommend using AggregateExpression",
              conditionMessage(m), fixed = TRUE)) {
      invokeRestart("muffleMessage")
    }
  })
}

.report_result_file_pattern <- function() {
  "^(AllMarkers|DEG_Counts_|AvgExpr_|Proportions_|PseudobulkDEG_|PseudobulkCounts_|GSEA_|Correlations_|Escape_Heatmap_|Escape_Scores_|.*_FeaturePlot_|VlnPlot_)"
}

.is_report_result_file <- function(paths) {
  grepl(.report_result_file_pattern(), basename(paths))
}

.snapshot_result_files <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[.is_report_result_file(files)]
  files <- normalizePath(files, winslash = "/", mustWork = FALSE)
  info <- file.info(files)
  snapshot <- data.frame(size = info$size, mtime = info$mtime,
                         checksum = unname(tools::md5sum(files)))
  rownames(snapshot) <- files
  snapshot
}

.changed_result_files <- function(snapshot, root) {
  current <- .snapshot_result_files(root)
  if (nrow(current) == 0L) return(character(0))
  old <- match(rownames(current), rownames(snapshot))
  is_new <- is.na(old)
  changed <- rep(FALSE, length(old))
  changed[!is_new] <- current$size[!is_new] != snapshot$size[old[!is_new]] |
    current$mtime[!is_new] != snapshot$mtime[old[!is_new]] |
    current$checksum[!is_new] != snapshot$checksum[old[!is_new]]
  normalizePath(rownames(current)[is_new | changed], winslash = "/",
                mustWork = FALSE)
}

.report_file_href <- function(path, report_dir) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    stop("path must be one existing file path.", call. = FALSE)
  }
  if (!is.character(report_dir) || length(report_dir) != 1L ||
      is.na(report_dir) || !nzchar(report_dir)) {
    stop("report_dir must be one existing directory path.", call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  report_dir <- normalizePath(report_dir, winslash = "/", mustWork = TRUE)
  path_parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
  dir_parts <- strsplit(report_dir, "/", fixed = TRUE)[[1L]]
  comparable_path <- if (.Platform$OS.type == "windows") tolower(path_parts) else
    path_parts
  comparable_dir <- if (.Platform$OS.type == "windows") tolower(dir_parts) else
    dir_parts
  common <- 0L
  limit <- min(length(path_parts), length(dir_parts))
  while (common < limit &&
         identical(comparable_path[common + 1L], comparable_dir[common + 1L])) {
    common <- common + 1L
  }
  if (common == 0L) {
    stop("Report and linked file must be on the same filesystem root.",
         call. = FALSE)
  }
  remaining <- length(path_parts) - common
  relative <- c(rep("..", length(dir_parts) - common),
                if (remaining > 0L) path_parts[(common + 1L):length(path_parts)])
  utils::URLencode(paste(relative, collapse = "/"), reserved = FALSE)
}

.report_file_link <- function(path, report_dir, label = basename(path)) {
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    stop("Report file links require package 'htmltools'.", call. = FALSE)
  }
  htmltools::tags$a(label, href = .report_file_href(path, report_dir),
                    download = basename(path))
}

.report_result_file_pattern <- function() {
  "^(AllMarkers|DEG_Counts_|AvgExpr_|Proportions_|PseudobulkDEG_|PseudobulkCounts_|GSEA_|Correlations_|Escape_Heatmap_|Escape_Scores_|.*_FeaturePlot_|VlnPlot_)"
}

.is_report_result_file <- function(paths) {
  grepl(.report_result_file_pattern(), basename(paths))
}

.append_report_download_section <- function(report_path, files, title) {
  if (!file.exists(report_path) || dir.exists(report_path)) {
    stop("report_path must identify an existing file.", call. = FALSE)
  }
  if (!is.character(title) || length(title) != 1L || is.na(title) ||
      !nzchar(title)) {
    stop("title must be one non-empty string.", call. = FALSE)
  }
  if (!is.character(files)) {
    stop("files must be a character vector of file paths.", call. = FALSE)
  }
  files <- unique(files[!is.na(files) & file.exists(files)])
  file_type <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("png", "jpg", "jpeg", "svg", "pdf")) "Plots" else
      if (ext %in% c("csv", "tsv", "txt")) "Tables" else "Other files"
  }
  groups <- split(files, vapply(files, file_type, character(1)))
  categories <- lapply(names(groups), function(category) {
    links <- lapply(groups[[category]], function(path) {
      rel_path <- .report_file_href(path, dirname(report_path))
      label <- paste0(basename(path), " — ", dirname(rel_path))
      htmltools::tags$li(
        `data-search` = tolower(paste(label, rel_path)),
        .report_file_link(path, dirname(report_path), label)
      )
    })
    htmltools::tags$details(
      class = "sparks-download-category",
      htmltools::tags$summary(paste0(category, " (", length(links), ")")),
      htmltools::tags$ul(links)
    )
  })
  section <- if (length(files)) htmltools::renderTags(htmltools::tags$section(
    id = "sparks-downloads",
    class = "sparks-download-ribbon",
    htmltools::tags$style(htmltools::HTML(
      "#sparks-downloads{position:sticky;top:0;z-index:1000;background:#fff;border:1px solid #cbd5e1;border-radius:.5rem;padding:.45rem .75rem;margin:0 0 1rem;box-shadow:0 2px 8px #0002}#sparks-downloads>details>summary{cursor:pointer;font-weight:600;list-style:none}#sparks-downloads>details>summary::-webkit-details-marker{display:none}#sparks-downloads>details>summary:before{content:'▸';display:inline-block;width:1.2em}#sparks-downloads>details[open]>summary:before{content:'▾'}#sparks-downloads .sparks-download-panel{max-height:65vh;overflow:auto;padding:.75rem 0 .25rem}#sparks-downloads input{box-sizing:border-box;width:min(100%,32rem);padding:.45rem;margin:0 0 .5rem}#sparks-downloads .sparks-download-category{margin:.25rem 0}#sparks-downloads .sparks-download-category>summary{cursor:pointer;font-weight:600}#sparks-downloads ul{margin:.25rem 0;padding-left:1.5rem}#sparks-downloads li{padding:.15rem 0;overflow-wrap:anywhere}"
    )),
    htmltools::tags$details(
      htmltools::tags$summary(paste0(title, " (", length(files), ")")),
      htmltools::tags$div(
        class = "sparks-download-panel",
        htmltools::tags$label(
          `for` = "sparks-download-search", "Filter files"
        ),
        htmltools::tags$input(
          id = "sparks-download-search", type = "search",
          placeholder = "Search file names or folders",
          `aria-label` = "Filter downloadable files"
        ),
        categories
      )
    ),
    htmltools::tags$script(htmltools::HTML(
      "document.addEventListener('input',function(e){if(e.target.id!=='sparks-download-search')return;var q=e.target.value.toLowerCase();document.querySelectorAll('#sparks-downloads li[data-search]').forEach(function(li){li.hidden=!li.dataset.search.includes(q)});document.querySelectorAll('#sparks-downloads .sparks-download-category').forEach(function(g){var shown=g.querySelectorAll('li:not([hidden])').length;g.hidden=shown===0;if(q&&shown)g.open=true;});});"
    ))
  ))$html else ""
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  html <- paste(readLines(report_path, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
  html <- gsub("(?s)<section id=\"sparks-downloads\".*?</section>", "", html,
               perl = TRUE)
  if (!grepl("</body>", html, fixed = TRUE)) {
    stop("report_path does not contain a closing </body> tag.", call. = FALSE)
  }
  if (!grepl("(?i)<body(?:\\s[^>]*)?>", html, perl = TRUE)) {
    stop("report_path does not contain an opening <body> tag.", call. = FALSE)
  }
  html <- sub("(?i)(<body(?:\\s[^>]*)?>)", paste0("\\1\n", section),
              html, perl = TRUE)
  writeLines(html, report_path, useBytes = TRUE)
  invisible(report_path)
}

#' Save a plot to PNG using ggsave or base R for pheatmap/gtable
#' @param p Plot object (ggplot or gtable)
#' @param filename Character. Output file path
#' @param width Numeric
#' @param height Numeric
#' @param dpi Numeric
#' @return NULL
#' @export
save_png <- function(p, filename, width = 8, height = 6, dpi = 300) {
  if (inherits(p, "gtable") || inherits(p, "Heatmap")) {
    grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
    # Close *this* device even if grid.draw() errors. Otherwise every failed
    # draw leaks an open png device (R allows 63) that later plots draw into.
    dev_id <- grDevices::dev.cur()
    on.exit(grDevices::dev.off(dev_id), add = TRUE)
    grid::grid.draw(p)
  } else {
    suppressMessages(
      ggplot2::ggsave(filename, plot = p, width = width, height = height, dpi = dpi, bg = "white")
    )
  }
  invisible(NULL)
}

#' Write table helper
#' @param x Data frame or matrix
#' @param file Character. Output file path
#' @param sep Character. Delimiter
#' @param quote Logical
#' @param row_names Logical
#' @return NULL
#' @export
write_table <- function(x, file, sep = "\t", quote = FALSE, row_names = FALSE) {
  utils::write.table(x, file = file, sep = sep, quote = quote, row.names = row_names)
  invisible(NULL)
}

#' Save Data Frame to TSV
#' @param x Data frame
#' @param file Character. Output file path
#' @return NULL
#' @export
save_tsv <- function(x, file) {
  write_table(x, file, sep = "\t", quote = FALSE, row_names = FALSE)
}

#' Get valid groups containing a minimum number of cells
#' @param meta Data frame (Seurat metadata)
#' @param col Character. Metadata column name
#' @param min_cells Integer
#' @return Character vector of valid groups
#' @export
get_valid_groups <- function(meta, col, min_cells = 3L) {
  if (!col %in% colnames(meta)) return(character(0))
  counts <- table(meta[[col]])
  names(counts)[counts >= min_cells]
}

#' Get Average Expression wrapper
#' @param seurat_obj Seurat object
#' @param layer Character. "data" or "counts"
#' @return Matrix of average expression or NULL
#' @export
get_avg_expr <- function(seurat_obj, layer = "data", assay = "SCT") {
  safe_run(.with_seurat_aggregate_advisory_muffled(
    Seurat::AverageExpression(seurat_obj, assays = assay,
                              layer = layer, return.seurat = FALSE)[[assay]]
  ), label = "AverageExpression")
}

#' Filter requested genes to those present in the object
#' @param genes Character vector
#' @param seurat_obj Seurat object
#' @param context Character for logging
#' @return Character vector of available genes, or NULL
.filter_present_genes <- function(genes, seurat_obj, context = "", assay = NULL) {
  if (is.null(genes) || length(genes) == 0L) return(NULL)
  available <- if (!is.null(assay)) {
    if (!assay %in% names(seurat_obj@assays)) character(0)
    else rownames(seurat_obj[[assay]])
  } else {
    rownames(seurat_obj)
  }
  genes <- unique(as.character(genes))
  valid <- intersect(genes, available)
  missing <- setdiff(genes, valid)
  if (length(missing) > 0L) {
    message("   [SKIP ", context, "] Requested genes absent",
            if (!is.null(assay)) paste0(" from assay '", assay, "'") else "",
            ": ", paste(missing, collapse = ", "))
  }
  if (length(valid) == 0L) {
    message("   [SKIP] ", context, ": none of the requested genes found in data.")
    return(NULL)
  }
  valid
}

#' Get regex pattern for junk genes dependent on species
#' @param species Character. "Mouse" or "Human"
#' @return Regex string
#' @export
get_junk_pattern <- function(species = "Mouse") {
  if (tolower(species) == "human") {
    # BUG FIX: the old pattern's bare "RNA" branch matched anything merely
    # *starting* with "RNA" -- RNASEH2A/B/C, RNASE1-13, RNASET2, RNASEK,
    # RNASEL, etc. are real, biologically meaningful genes, not junk, and were
    # being silently dropped from generate_expression_heatmap()'s "top
    # expressed genes" heatmap for every human run. Dropped the ill-specified
    # "RNR|RNA" branch and aligned this with the vetted pattern already used
    # in the human config template's species$gene_removal_pattern.
    "^(MT-|RPS|RPL|HBA[12]?|HBB$)"
  } else {
    "^(mt-|Rps|Rpl|Rrn|Rn|Hb|Gm).*|.*Rik$"
  }
}

#' Wrapper for FindAllMarkers to catch errors safely
#' @param seurat_obj Seurat object
#' @param only_pos Logical
#' @param min_pct Numeric
#' @param logfc_threshold Numeric
#' @return Data frame of markers or NULL
.find_all_markers_safe <- function(seurat_obj, only_pos = TRUE, min_pct = 0.25, logfc_threshold = 0.25) {
  safe_run(
    Seurat::FindAllMarkers(seurat_obj, only.pos = only_pos, min.pct = min_pct,
                           logfc.threshold = logfc_threshold, verbose = FALSE),
    label = "FindAllMarkers"
  )
}

#' Check if object has a specific reduction
#' @param seurat_obj Seurat object
#' @param reduction Character
#' @return Logical
.has_reduction <- function(seurat_obj, reduction) {
  !is.null(seurat_obj@reductions[[reduction]])
}

#' Check if object has scaled data
#' @param seurat_obj Seurat object
#' @return Logical
.has_scale_data <- function(seurat_obj) {
  assay <- Seurat::DefaultAssay(seurat_obj)
  has_data <- tryCatch({
    mat <- Seurat::LayerData(seurat_obj, assay = assay, layer = "scale.data")
    !is.null(mat) && nrow(mat) > 0
  }, error = function(e) FALSE)
  has_data
}

#' Adaptive label style for DoHeatmap group-bar text, based on group count
#'
#' With many groups (e.g. singleR_labels_fine, ~40-90 categories), horizontal
#' or 45-degree text overlaps and becomes unreadable. This scales the angle
#' toward vertical and shrinks the font as the group count grows.
#' @param n Integer. Number of distinct groups/columns being labeled
#' @return Named list: angle, hjust, size
.heatmap_label_params <- function(n) {
  if (n > 20) {
    list(angle = 90, hjust = 0, size = 2.5)
  } else if (n > 8) {
    list(angle = 90, hjust = 0, size = 3)
  } else if (n > 4) {
    list(angle = 90, hjust = 0, size = 3.5)
  } else {
    list(angle = 0, hjust = 0.5, size = 4)
  }
}

.save_rdata <- function(obj, dir, name) {
  make_dir(dir)
  saveRDS(obj, file = file.path(dir, paste0(name, ".rds")))
  invisible(NULL)
}

.setup_group_dirs <- function(results_dir, comp_group) {
  dirs <- list(
    base       = file.path(results_dir, comp_group),
    qc         = file.path(results_dir, comp_group, "QC"),
    rdata      = file.path(results_dir, comp_group, "RData"),
    raw_matrix = file.path(results_dir, comp_group, "RawMatrix")
  )
  lapply(dirs, make_dir)
  dirs
}

.make_analysis_dirs <- function(group_dir) {
  dirs <- list(
    UMAP        = file.path(group_dir, "UMAP"),
    DEG         = file.path(group_dir, "DEG"),
    Correlation = file.path(group_dir, "Correlation"),
    Pseudobulk  = file.path(group_dir, "Pseudobulk"),
    Expression  = file.path(group_dir, "Expression"),
    Proportions = file.path(group_dir, "Proportions"),
    VlnPlot     = file.path(group_dir, "VlnPlot"),
    Heatmap     = file.path(group_dir, "Heatmap"),
    GSEA        = file.path(group_dir, "GSEA")
  )
  lapply(dirs[names(dirs) != "GSEA"], make_dir)
  dirs
}

.merge_samples <- function(obj_list, folder_ids, add_cell_ids = TRUE) {
  obj_list <- Filter(Negate(is.null), obj_list)
  if (length(obj_list) == 0L) stop("No valid samples to merge.")
  if (length(obj_list) == 1L) return(obj_list[[1L]])
  merge_args <- list(x = obj_list[[1L]], y = obj_list[-1L])
  if (isTRUE(add_cell_ids)) merge_args$add.cell.ids <- as.character(folder_ids)
  obj <- do.call(merge, merge_args)
  SeuratObject::JoinLayers(obj)
}

#' Environment in which the report templates are evaluated
#'
#' A child of the global environment (so the Rmd sees whatever the user has
#' attached) that also carries SPARKS's own \code{\%||\%}: the operator is
#' exported but is not on the search path when the package is only loaded
#' (\code{SPARKS::sparks()}), and base R only provides it from 4.4.0.
#' @return A new environment
.report_env <- function() {
  env <- new.env(parent = globalenv())
  assign("%||%", `%||%`, envir = env)
  env
}

.render_report <- function(template, output_file, params) {
  rmarkdown::render(
    input = template,
    output_file = output_file,
    params = params,
    envir = .report_env(),
    quiet = TRUE
  )
}

#' Generate an HTML QC report for a processed Seurat object
#'
#' The report is always written to \code{out_dir} (the group's QC folder).
#'
#' @param seurat_obj A Seurat object after full processing.
#' @param comp_group Character. Comparison group name.
#' @param out_dir Character. Output directory – the HTML is saved here.
#' @param author Character. Author name shown in the report.
#' @param title Character. Report title.
#' @param rmd_template Character. Path to the qc_report.Rmd template.
#' @param cfg Named list. Full pipeline config (passed as a param so the Rmd
#'   can display QC thresholds). Optional – defaults to NULL.
#' @return Invisibly, the path to the generated HTML file.
#' @export
generate_qc_report <- function(seurat_obj, comp_group, out_dir,
                               author       = "Pipeline",
                               title        = NULL,
                               rmd_template = NULL,
                               cfg          = NULL,
                               run_files = NULL,
                               log_path = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is needed. Please install it: install.packages('rmarkdown')")

  if (is.null(title)) title <- paste("QC Report -", comp_group)

  if (is.null(rmd_template) || !file.exists(rmd_template)) {
    message("   [WARNING] generate_qc_report: template not found at '",
            rmd_template %||% "<NULL>", "'. Skipping HTML report.")
    return(invisible(NULL))
  }

  make_dir(out_dir)
  output_file <- normalizePath(
    file.path(out_dir, paste0("QC_report_", comp_group, ".html")),
    mustWork = FALSE
  )

  .render_report(
    template = rmd_template,
    output_file = output_file,
    params = list(
      seurat_obj = seurat_obj,
      comp_group = comp_group,
      author     = author,
      title      = title,
      out_dir    = normalizePath(out_dir, mustWork = FALSE),
      cfg        = cfg,
      run_files = run_files,
      log_path = log_path
    )
  )
  doublet_files <- c(run_files, list.files(
    file.path(dirname(dirname(out_dir)), "Sample_QC"), recursive = TRUE,
    pattern = "^Doublet_Removed_.*[.]txt$", full.names = TRUE
  ))
  doublet_files <- doublet_files[!is.na(doublet_files) &
    grepl("^Doublet_Removed_.*[.]txt$",
                                       basename(doublet_files))]
  .append_report_download_section(output_file, doublet_files,
                                  "Per-sample doublet statistics")

  message("   QC report saved to: ", output_file)
  invisible(output_file)
}

#' Generate an interactive HTML results report for a comparison group
#'
#' Renders \code{results_report.Rmd} into \code{out_dir} as
#' \code{Results_report_{comp_group}.html}.
#'
#' @param seurat_obj A Seurat object after full processing.
#' @param comp_group Character. Comparison group name.
#' @param out_dir Character. Output directory — the HTML is saved here.
#' @param groupings Character vector. Metadata columns used for grouping
#'   (e.g. "seurat_clusters", "singleR_labels_main"). Passed to the Rmd so it
#'   can render one interactive UMAP / table per grouping.
#' @param author Character. Author name shown in the report.
#' @param title Character. Report title.
#' @param rmd_template Character. Path to results_report.Rmd.
#' @param cfg Named list. Full pipeline config.
#' @return Invisibly, the path to the generated HTML file.
#' @export
generate_results_report <- function(seurat_obj, comp_group, out_dir,
                                    groupings    = NULL,
                                    author       = "Pipeline",
                                    title        = NULL,
                                    rmd_template = NULL,
                                    cfg          = NULL,
                                    run_files = NULL,
                                    log_path = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is needed. Please install it: install.packages('rmarkdown')")

  if (is.null(title)) title <- paste("Results Report -", comp_group)

  if (is.null(rmd_template) || !file.exists(rmd_template)) {
    message("   [WARNING] generate_results_report: template not found at '",
            rmd_template %||% "<NULL>", "'. Skipping results report.")
    return(invisible(NULL))
  }

  make_dir(out_dir)

  output_file <- normalizePath(
    file.path(out_dir, paste0("Results_report_", comp_group, ".html")),
    mustWork = FALSE
  )

  .render_report(
    template = rmd_template,
    output_file = output_file,
    params = list(
      seurat_obj = seurat_obj,
      comp_group = comp_group,
      author     = author,
      title      = title,
      out_dir    = normalizePath(out_dir, mustWork = FALSE),
      cfg        = cfg,
      groupings  = groupings,
      run_files = run_files,
      log_path = log_path
    )
  )
  report_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
  report_files <- report_files[!file.info(report_files)$isdir %in% TRUE]
  report_files <- report_files[.is_report_result_file(report_files)]
  .append_report_download_section(output_file, report_files,
                                  "Download result files")

  message("   Results report saved to: ", output_file)
  invisible(output_file)
}
