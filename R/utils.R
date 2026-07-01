# ── General Utilities ─────────────────────────────────────────────────────────

#' Null coalescing operator
#' @name %||%
#' @param lhs Left hand side
#' @param rhs Right hand side
#' @export
`%||%` <- function(lhs, rhs) {
  if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
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
#' @export
.safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]", "_", x)
}

# ── File I/O Utilities ────────────────────────────────────────────────────────

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
    grid::grid.draw(p)
    grDevices::dev.off()
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

# ── Seurat Data Extractors & Checks ───────────────────────────────────────────

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
get_avg_expr <- function(seurat_obj, layer = "data") {
  safe_run({
    suppressWarnings(
      Seurat::AverageExpression(seurat_obj, assays = "SCT",
                                layer = layer, return.seurat = FALSE)[["SCT"]]
    )
  }, label = "AverageExpression")
}

#' Filter requested genes to those present in the object
#' @param genes Character vector
#' @param seurat_obj Seurat object
#' @param context Character for logging
#' @return Character vector of available genes, or NULL
#' @export
.filter_present_genes <- function(genes, seurat_obj, context = "") {
  if (is.null(genes) || length(genes) == 0L) return(NULL)
  valid <- intersect(as.character(genes), rownames(seurat_obj))
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
    "^(MT-|RPS|RPL|RNR|RNA)"
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
#' @export
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
#' @export
.has_reduction <- function(seurat_obj, reduction) {
  !is.null(seurat_obj@reductions[[reduction]])
}

#' Check if object has scaled data
#' @param seurat_obj Seurat object
#' @return Logical
#' @export
.has_scale_data <- function(seurat_obj) {
  assay <- Seurat::DefaultAssay(seurat_obj)
  has_data <- tryCatch({
    mat <- Seurat::LayerData(seurat_obj, assay = assay, layer = "scale.data")
    !is.null(mat) && nrow(mat) > 0
  }, error = function(e) FALSE)
  has_data
}

.save_rdata <- function(obj, dir, name) {
  make_dir(dir)
  saveRDS(obj, file = file.path(dir, paste0(name, ".rds")))
  invisible(NULL)
}

.setup_group_dirs <- function(results_dir, comp_group) {
  dirs <- list(
    base  = file.path(results_dir, comp_group),
    qc    = file.path(results_dir, comp_group, "QC"),
    rdata = file.path(results_dir, comp_group, "RData")
  )
  lapply(dirs, make_dir)
  dirs
}

.make_analysis_dirs <- function(group_dir) {
  dirs <- list(
    UMAP    = file.path(group_dir, "UMAP"),
    DEG     = file.path(group_dir, "DEG"),
    VlnPlot = file.path(group_dir, "VlnPlot"),
    Heatmap = file.path(group_dir, "Heatmap")
  )
  lapply(dirs, make_dir)
  dirs
}

.merge_samples <- function(obj_list, folder_ids) {
  obj_list <- Filter(Negate(is.null), obj_list)
  if (length(obj_list) == 0L) stop("No valid samples to merge.")
  if (length(obj_list) == 1L) return(obj_list[[1L]])
  obj <- merge(obj_list[[1L]], y = obj_list[-1L],
               add.cell.ids = as.character(folder_ids))
  SeuratObject::JoinLayers(obj)
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
                               cfg          = NULL) {

  if (!requireNamespace("rmarkdown", quietly = TRUE))
    stop("Package 'rmarkdown' is needed. Please install it: install.packages('rmarkdown')")

  if (is.null(title)) title <- paste("QC Report -", comp_group)

  # BUG FIX #9: the original code had an unconditional stop() inside the
  # "template not found" branch, which prevented any report from ever being
  # generated.  We now emit a warning and return invisibly so the pipeline
  # continues.  The caller (main.R) already validates the template path before
  # reaching this function, so this branch is only a safety net.
  if (is.null(rmd_template) || !file.exists(rmd_template)) {
    message("   [WARNING] generate_qc_report: template not found at '",
            rmd_template %||% "<NULL>", "'. Skipping HTML report.")
    return(invisible(NULL))
  }

  make_dir(out_dir)

  # BUG FIX #12: rmarkdown::render() changes the working directory to the
  # location of the input Rmd during knitting.  If output_file is a *relative*
  # path it is therefore resolved against the Rmd's directory, not the caller's
  # working directory — so the HTML ended up in the package/inst folder instead
  # of the QC output folder.  normalizePath(..., mustWork = FALSE) converts the
  # path to an absolute one before render() is called, guaranteeing the file
  # always lands in out_dir regardless of where the Rmd lives.
  output_file <- normalizePath(
    file.path(out_dir, paste0("QC_report_", comp_group, ".html")),
    mustWork = FALSE
  )

  rmarkdown::render(
    input       = rmd_template,
    output_file = output_file,
    params      = list(
      seurat_obj = seurat_obj,
      comp_group = comp_group,
      author     = author,
      title      = title,
      out_dir    = normalizePath(out_dir, mustWork = FALSE),
      cfg        = cfg
    ),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )

  message("   QC report saved to: ", output_file)
  invisible(output_file)
}

# ──────────────────────────────────────────────────────────────────────────────
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
                                    cfg          = NULL) {

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

  rmarkdown::render(
    input       = rmd_template,
    output_file = output_file,
    params      = list(
      seurat_obj = seurat_obj,
      comp_group = comp_group,
      author     = author,
      title      = title,
      out_dir    = normalizePath(out_dir, mustWork = FALSE),
      cfg        = cfg,
      groupings  = groupings
    ),
    envir = new.env(parent = globalenv()),
    quiet = FALSE
  )

  message("   Results report saved to: ", output_file)
  invisible(output_file)
}
