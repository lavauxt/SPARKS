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
    envir = .report_env(),
    quiet = FALSE
  )

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
    envir = .report_env(),
    quiet = FALSE
  )

  message("   Results report saved to: ", output_file)
  invisible(output_file)
}
