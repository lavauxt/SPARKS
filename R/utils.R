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