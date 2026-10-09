#' Recursively merge an override config onto a base config
#'
#' YAML mappings (named lists) are merged key by key; everything else --
#' scalars, vectors and YAML \emph{sequences} (unnamed lists such as
#' \code{subsets}, \code{gene_signatures} or \code{singler$labels}) -- is
#' replaced wholesale by the override. \code{NULL} values in the override are
#' kept as \code{NULL} (a key set to \code{~} stays present).
#'
#' \code{utils::modifyList()} cannot do this: for two unnamed lists it recurses,
#' finds no names in the override and returns the \emph{base} list unchanged. An
#' override \code{subsets:} / \code{gene_signatures:} block was therefore silently
#' ignored whenever the base template defined one, and \code{subsets: []} could
#' never clear the template's subsets.
#' @param base Named list
#' @param override Named list
#' @return Merged named list
#' @keywords internal
.merge_config <- function(base, override) {
  is_mapping <- function(x) {
    is.list(x) && !is.null(names(x)) && all(nzchar(names(x)))
  }
  # An empty list where the base has a mapping overrides nothing (e.g. an empty
  # override file); where the base has a *sequence* it clears it (`subsets: []`).
  if (is_mapping(base) && is.list(override) && length(override) == 0L) return(base)
  if (!is_mapping(base) || !is_mapping(override)) return(override)
  for (key in names(override)) {
    base[key] <- list(
      if (key %in% names(base)) .merge_config(base[[key]], override[[key]])
      else override[[key]]
    )
  }
  base
}

#' Load and merge a base config and an optional override config
#' @param base_config_path Character. Path to base YAML (template)
#' @param override_config_path Character. Path to override YAML (your file)
#' @return Merged named list
#' @export
load_pipeline_config <- function(base_config_path, override_config_path = NULL) {

  if (!file.exists(base_config_path)) stop("Base config not found: ", base_config_path)
  cfg <- yaml::yaml.load_file(base_config_path)

  if (!is.null(override_config_path) && file.exists(override_config_path)) {
    override_cfg <- yaml::yaml.load_file(override_config_path)
    if (is.null(override_cfg)) override_cfg <- list()   # empty override file
    if (!is.list(override_cfg))
      stop("Override config must be a YAML mapping: ", override_config_path)
    # Keys the override sets to null/~ are kept as NULL (not deleted), e.g.
    # `groupings_main: ~` (force auto-detection) must not drop the key and trip
    # the required-keys check below.
    # BUG FIX: this used to be utils::modifyList(cfg, override_cfg, keep.null =
    # TRUE), which silently ignores any override *list* (subsets,
    # gene_signatures, singler$labels, ...) when the base also defines one --
    # see .merge_config().
    cfg <- .merge_config(cfg, override_cfg)
  }


  cfg$pipeline$config_dir <- dirname(normalizePath(base_config_path, mustWork = FALSE))
  if (!is.null(override_config_path) && file.exists(override_config_path)) {
    cfg$pipeline$override_config_dir <- dirname(
      normalizePath(override_config_path, mustWork = FALSE)
    )
  }


  required <- c("pipeline", "qc", "processing", "deg", "plot",
                "species", "singler", "groupings_main")
  missing  <- setdiff(required, names(cfg))
  if (length(missing) > 0L)
    stop("Config missing required keys: ", paste(missing, collapse = ", "))

  # pipeline$independent_samples: TRUE analyses every row of the sample table
  # on its own (own normalisation, clustering, UMAP, annotation, output folder)
  # instead of merging the samples of a comparison_group.
  cfg$pipeline$independent_samples <- cfg$pipeline$independent_samples %||% FALSE

  cfg$processing$cluster_col        <- cfg$processing$cluster_col        %||% "seurat_clusters"
  cfg$processing$condition_col      <- cfg$processing$condition_col      %||% "condition"
  # process_single_sample() stores the condition in the metadata column
  # "condition" and the DEG / proportion / plotting helpers read that column by
  # name, so any other value only fails later with an obscure error.
  if (!identical(cfg$processing$condition_col, "condition")) {
    warning("processing$condition_col = '", cfg$processing$condition_col,
            "' is not supported (the pipeline always stores the condition in ",
            "'condition'); using 'condition'.", call. = FALSE)
    cfg$processing$condition_col <- "condition"
  }
  cfg$processing$reduction          <- cfg$processing$reduction          %||% "umap"
  cfg$processing$cluster_resolution <- cfg$processing$cluster_resolution %||% 0.5
  cfg$processing$npcs               <- cfg$processing$npcs               %||% 50L
  cfg$processing$seed               <- cfg$processing$seed               %||% 2861L
  cfg$processing$seed <- tryCatch(
    .validate_seed(cfg$processing$seed),
    error = function(e) stop("Invalid processing$seed: ", conditionMessage(e), call. = FALSE)
  )
  cfg$processing$n_elbow_dims       <- cfg$processing$n_elbow_dims       %||% 30L
  cfg$processing$sct_assay          <- cfg$processing$sct_assay          %||% "SCT"
  cfg$processing$vars_to_regress    <- cfg$processing$vars_to_regress    %||% "percent.mt"
  # Re-express the SCT corrected counts/data at a common depth (Seurat's
  # PrepSCTFindMarkers) right after SCTransform, so that cross-sample DE and
  # average-expression tables are comparable. See run_seurat_processing().
  cfg$processing$prep_sct_findmarkers <- cfg$processing$prep_sct_findmarkers %||% TRUE


  if (is.null(cfg$processing$pca_dims)) {
    from <- cfg$processing$pca_dims_from %||% 1L
    to   <- cfg$processing$pca_dims_to   %||% 20L
    cfg$processing$pca_dims <- seq.int(from, to)
  }

  cfg$sex_scoring$run      <- cfg$sex_scoring$run      %||% FALSE
  cfg$sex_scoring$regress  <- cfg$sex_scoring$regress  %||% FALSE
  cfg$sex_scoring$markers  <- cfg$sex_scoring$markers  %||% list(female = c(), male = c())

  cfg$cell_cycle$run     <- cfg$cell_cycle$run     %||% FALSE
  cfg$cell_cycle$regress <- cfg$cell_cycle$regress %||% FALSE

  cfg$species$mt_pattern           <- cfg$species$mt_pattern           %||% "^mt-"
  cfg$species$gene_removal_pattern <- cfg$species$gene_removal_pattern %||% "^(mt-|Rps|Rpl|Rrn|Rn|Hb|Gm).*|.*Rik$"
  cfg$species$genes_to_remove      <- cfg$species$genes_to_remove      %||% c()

  # Input format: auto | alevin | 10x | 10x_h5 | starsolo | h5ad | mtx.
  # See .load_counts() and export_raw_matrix() in processing.R.
  cfg$input$format            <- cfg$input$format            %||% "auto"
  cfg$input$starsolo_feature  <- cfg$input$starsolo_feature  %||% "Gene"
  cfg$input$h5ad_raw_slot     <- cfg$input$h5ad_raw_slot     %||% "auto"
  cfg$input$save_raw_matrix   <- cfg$input$save_raw_matrix   %||% TRUE
  cfg$input$raw_matrix_format <- cfg$input$raw_matrix_format %||% "mtx"

  cfg$qc$min_features   <- cfg$qc$min_features   %||% 200L
  cfg$qc$max_features   <- cfg$qc$max_features   %||% 6000L
  cfg$qc$max_counts     <- cfg$qc$max_counts     %||% 30000L
  cfg$qc$max_mt_percent <- cfg$qc$max_mt_percent %||% 10
  cfg$qc$min_cells      <- cfg$qc$min_cells      %||% 3L

  # deg$run: FALSE disables every between-condition comparison (DEG tables,
  # DEG-count UMAP, chi-squared / scProportionTest). Clustering, annotation,
  # cluster markers and descriptive plots are unaffected.
  cfg$deg$run                   <- cfg$deg$run                   %||% TRUE
  cfg$deg$logfc_threshold       <- cfg$deg$logfc_threshold       %||% 0.25
  cfg$deg$min_pct               <- cfg$deg$min_pct               %||% 0.1
  cfg$deg$min_p_val_adj         <- cfg$deg$min_p_val_adj         %||% 0.05
  cfg$deg$min_cells_per_group   <- cfg$deg$min_cells_per_group   %||% 10L
  cfg$deg$min_deg_display       <- cfg$deg$min_deg_display       %||% 5L
  cfg$deg$avg_expression_layers <- cfg$deg$avg_expression_layers %||% list("data")
  cfg$deg$table_sep             <- cfg$deg$table_sep             %||% "\t"
  cfg$deg$table_quote           <- cfg$deg$table_quote           %||% FALSE
  cfg$deg$table_row_names       <- cfg$deg$table_row_names       %||% FALSE

  cfg$plot$top_genes_heatmap_n  <- cfg$plot$top_genes_heatmap_n  %||% 10L
  cfg$plot$umap_width_standard  <- cfg$plot$umap_width_standard  %||% 14
  cfg$plot$umap_height_standard <- cfg$plot$umap_height_standard %||% 7
  cfg$plot$umap_width_fine      <- cfg$plot$umap_width_fine      %||% 20
  cfg$plot$umap_height_fine     <- cfg$plot$umap_height_fine     %||% 10
  cfg$plot$legend_nrow_fine     <- cfg$plot$legend_nrow_fine     %||% 6L
  cfg$plot$elbow_width          <- cfg$plot$elbow_width          %||% 8
  cfg$plot$elbow_height         <- cfg$plot$elbow_height         %||% 5
  cfg$plot$qc_vln_width         <- cfg$plot$qc_vln_width         %||% 12
  cfg$plot$qc_vln_height        <- cfg$plot$qc_vln_height        %||% 5
  cfg$plot$qc_scatter_width     <- cfg$plot$qc_scatter_width     %||% 10
  cfg$plot$qc_scatter_height    <- cfg$plot$qc_scatter_height    %||% 5
  cfg$plot$deg_color_main       <- cfg$plot$deg_color_main       %||% "red"
  cfg$plot$deg_umap_color_low   <- cfg$plot$deg_umap_color_low   %||% "lightgrey"
  cfg$plot$deg_umap_point_size  <- cfg$plot$deg_umap_point_size  %||% 0.5
  cfg$plot$deg_umap_alpha       <- cfg$plot$deg_umap_alpha       %||% 0.6
  cfg$plot$deg_umap_label_size  <- cfg$plot$deg_umap_label_size  %||% 4
  cfg$plot$deg_umap_width       <- cfg$plot$deg_umap_width       %||% 8
  cfg$plot$deg_umap_height      <- cfg$plot$deg_umap_height      %||% 7

  cfg$singler$unassigned_label    <- cfg$singler$unassigned_label    %||% "Unassigned"
  cfg$singler$min_cells_per_group <- cfg$singler$min_cells_per_group %||% 10L
  cfg$singler$labels              <- cfg$singler$labels              %||% list()
  # BUG FIX: sapply() on an empty labels list returns list() (not
  # character(0)). c(character_vector, list()) then silently coerces the
  # whole result to a list -- groupings_main below, and the analogous
  # groupings built in .run_subset() -- which quietly breaks downstream
  # %in%/paste() comparisons whenever singler$labels is empty (e.g. someone
  # disables SingleR for a quick test run). Confirmed empirically. vapply
  # with an explicit empty-case guard keeps this a character vector always.
  cfg$singler$label_names <- if (length(cfg$singler$labels) > 0L) {
    vapply(cfg$singler$labels, `[[`, character(1L), "name")
  } else {
    character(0)
  }

  cfg$escape$run      <- cfg$escape$run      %||% FALSE
  cfg$escape$method   <- cfg$escape$method   %||% "ssGSEA"
  cfg$escape$library  <- cfg$escape$library  %||% "H"
  cfg$escape$min_size <- cfg$escape$min_size %||% 5

  cfg$labeling$min_subset_cells          <- cfg$labeling$min_subset_cells          %||% 50L
  cfg$labeling$unassigned_suffix         <- cfg$labeling$unassigned_suffix         %||% "_Unassigned"
  cfg$labeling$marker_positive_threshold <- cfg$labeling$marker_positive_threshold %||% 0.1
  cfg$subsets            <- cfg$subsets            %||% list()

  if (is.null(cfg$genes)) cfg$genes <- list()
  cfg$genes$genes_to_plot <- cfg$genes$genes_to_plot %||% NULL
  cfg$genes$corr_genes_x  <- cfg$genes$corr_genes_x  %||% NULL
  cfg$genes$corr_genes_y  <- cfg$genes$corr_genes_y  %||% NULL

  cfg$gene_signatures <- cfg$gene_signatures %||% list()
  if (length(cfg$gene_signatures) > 0L) {
    bad <- vapply(cfg$gene_signatures, function(s) {
      is.null(s$name) || is.null(s$genes) || length(s$genes) == 0L
    }, logical(1))
    if (any(bad)) {
      warning("gene_signatures entries missing 'name' or 'genes' will be skipped: ",
              paste(which(bad), collapse = ", "))
      cfg$gene_signatures <- cfg$gene_signatures[!bad]
    }
  }

  cfg$parallel$enable      <- cfg$parallel$enable      %||% FALSE
  cfg$parallel$workers     <- cfg$parallel$workers     %||% 4L
  cfg$parallel$strategy    <- cfg$parallel$strategy    %||% "multisession"
  cfg$parallel$max_size_gb <- cfg$parallel$max_size_gb %||% 8.0

  return(cfg)
}

#' Load sample table from TSV or inline YAML list
#' @param cfg Named list from load_pipeline_config
#' @return Data frame with columns: folder_id, protocol, comparison_group
#' @export
load_sample_table <- function(cfg) {
  st <- cfg$pipeline$sample_table

  if (is.null(st)) {
    stop("pipeline$sample_table is missing from the configuration.")
  }

  if (is.character(st) && length(st) == 1L) {
    sample_table_path <- st
    if (!file.exists(sample_table_path)) {
      config_dirs <- c(
        cfg$pipeline$override_config_dir %||% character(0),
        cfg$pipeline$config_dir %||% character(0)
      )
      config_dirs <- config_dirs[
        vapply(config_dirs, is.character, logical(1)) &
          nzchar(config_dirs)
      ]
      ancestor_dirs <- function(paths) {
        paths <- paths[!is.na(paths) & nzchar(paths)]
        ancestors <- character(0)
        for (path in paths) {
          current <- normalizePath(path, winslash = "/", mustWork = FALSE)
          repeat {
            ancestors <- c(ancestors, current)
            parent <- dirname(current)
            if (identical(parent, current)) break
            current <- parent
          }
        }
        unique(ancestors)
      }

      candidate_dirs <- ancestor_dirs(config_dirs)
      candidate_paths <- file.path(candidate_dirs, sample_table_path)
      sample_table_path <- candidate_paths[file.exists(candidate_paths)][1L]
    }

    if (!is.na(sample_table_path) && file.exists(sample_table_path)) {
      df <- utils::read.delim(
        sample_table_path,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      stop("pipeline$sample_table must be a valid file path or a YAML list.")
    }
  } else if (is.list(st)) {
    df <- do.call(dplyr::bind_rows, lapply(st, as.data.frame, stringsAsFactors = FALSE))
  } else {
    stop("pipeline$sample_table must be a valid file path or a YAML list.")
  }

  if (nrow(df) == 0L) {
    stop("The loaded sample table is empty.")
  }

  required_cols <- c("folder_id", "protocol", "comparison_group")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0L) {
    stop("Sample table missing columns: ", paste(missing, collapse = ", "))
  }

  df$folder_id        <- as.character(df$folder_id)
  df$protocol         <- as.character(df$protocol)
  df$comparison_group <- as.character(df$comparison_group)
  if ("quant_format" %in% colnames(df)) {
    df$quant_format <- as.character(df$quant_format)
  }

  return(df)
}
