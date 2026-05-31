#' Load and merge a base config and an optional override config
#' @param base_config_path Character. Path to base YAML (template)
#' @param override_config_path Character. Path to override YAML (your file)
#' @return Merged named list
#' @export
load_pipeline_config <- function(base_config_path, override_config_path = NULL) {
  
  # 1. Load the base template
  if (!file.exists(base_config_path)) stop("Base config not found: ", base_config_path)
  cfg <- yaml::yaml.load_file(base_config_path)
  
  # 2. Merge overrides if provided
  if (!is.null(override_config_path) && file.exists(override_config_path)) {
    override_cfg <- yaml::yaml.load_file(override_config_path)
    # modifyList recursively replaces base values with override values
    cfg <- utils::modifyList(cfg, override_cfg)
  }

  # ── Required top-level keys ────────────────────────────────────────────────
  required <- c("pipeline", "qc", "processing", "deg", "plot",
                "species", "singler", "genes", "groupings_main")
  missing  <- setdiff(required, names(cfg))
  if (length(missing) > 0L)
    stop("Config missing required keys: ", paste(missing, collapse = ", "))

  # ── Processing Defaults ───────────────────────────────────────────────────
  cfg$processing$cluster_col        <- cfg$processing$cluster_col        %||% "seurat_clusters"
  cfg$processing$condition_col      <- cfg$processing$condition_col      %||% "condition"
  cfg$processing$reduction          <- cfg$processing$reduction          %||% "umap"
  cfg$processing$pca_dims           <- cfg$processing$pca_dims           %||% 1:20
  cfg$processing$cluster_resolution <- cfg$processing$cluster_resolution %||% 0.5
  cfg$processing$npcs               <- cfg$processing$npcs               %||% 50L
  cfg$processing$n_elbow_dims       <- cfg$processing$n_elbow_dims       %||% 30L
  cfg$processing$sct_assay          <- cfg$processing$sct_assay          %||% "SCT"
  cfg$processing$vars_to_regress    <- cfg$processing$vars_to_regress    %||% "percent.mt"

  # ── Sex Scoring Defaults ─────────────────────────────────────────────
  cfg$sex_scoring$run      <- cfg$sex_scoring$run      %||% FALSE
  cfg$sex_scoring$regress  <- cfg$sex_scoring$regress  %||% FALSE
  cfg$sex_scoring$markers  <- cfg$sex_scoring$markers  %||% list(female = c(), male = c())
  
  # ── Cell Cycle Defaults ───────────────────────────────────────────────────
  cfg$cell_cycle$run     <- cfg$cell_cycle$run     %||% FALSE
  cfg$cell_cycle$regress <- cfg$cell_cycle$regress %||% FALSE

  # ── Species Defaults ─────────────────────────────────────────────────────
  cfg$species$mt_pattern           <- cfg$species$mt_pattern           %||% "^mt-"
  cfg$species$gene_removal_pattern <- cfg$species$gene_removal_pattern %||% "^(mt-|Rps|Rpl|Rrn|Rn|Hb|Gm).*|.*Rik$"
  cfg$species$genes_to_remove      <- cfg$species$genes_to_remove      %||% c()

  # ── QC Defaults ───────────────────────────────────────────────────────────
  cfg$qc$min_features   <- cfg$qc$min_features   %||% 200L
  cfg$qc$max_features   <- cfg$qc$max_features   %||% 6000L
  cfg$qc$max_counts     <- cfg$qc$max_counts     %||% 30000L
  cfg$qc$max_mt_percent <- cfg$qc$max_mt_percent %||% 10
  cfg$qc$min_cells      <- cfg$qc$min_cells      %||% 3L

  # ── DEG Defaults ──────────────────────────────────────────────────────────
  cfg$deg$logfc_threshold     <- cfg$deg$logfc_threshold     %||% 0.25
  cfg$deg$min_pct             <- cfg$deg$min_pct             %||% 0.1
  cfg$deg$min_p_val_adj       <- cfg$deg$min_p_val_adj       %||% 0.05
  cfg$deg$min_cells_per_group <- cfg$deg$min_cells_per_group %||% 10L
  cfg$deg$min_deg_display     <- cfg$deg$min_deg_display     %||% 5L
  cfg$deg$avg_expression_layers <- cfg$deg$avg_expression_layers %||% list("data")
  cfg$deg$table_sep           <- cfg$deg$table_sep           %||% "\t"
  cfg$deg$table_quote         <- cfg$deg$table_quote         %||% FALSE
  cfg$deg$table_row_names     <- cfg$deg$table_row_names     %||% FALSE

  # ── Plotting Defaults ─────────────────────────────────────────────────────
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

  # ── SingleR Defaults ──────────────────────────────────────────────────────
  cfg$singler$unassigned_label    <- cfg$singler$unassigned_label    %||% "Unassigned"
  cfg$singler$min_cells_per_group <- cfg$singler$min_cells_per_group %||% 10L
  cfg$singler$labels              <- cfg$singler$labels              %||% list()
  cfg$singler$label_names         <- sapply(cfg$singler$labels, `[[`, "name")

# ── Escape Defaults ───────────────────────────────────────────────────────
  cfg$escape$run      <- cfg$escape$run      %||% FALSE
  cfg$escape$method   <- cfg$escape$method   %||% "ssGSEA"
  cfg$escape$library  <- cfg$escape$library  %||% "H"       # MSigDB Hallmark
  cfg$escape$min_size <- cfg$escape$min_size %||% 5

  # ── Labeling & Gene Defaults ─────────────────────────────────────────────
  cfg$labeling$min_subset_cells          <- cfg$labeling$min_subset_cells          %||% 50L
  cfg$labeling$unassigned_suffix         <- cfg$labeling$unassigned_suffix         %||% "_Unassigned"
  cfg$labeling$marker_positive_threshold <- cfg$labeling$marker_positive_threshold %||% 0.1
  cfg$subsets            <- cfg$subsets            %||% list()
  cfg$genes$corr_genes_x <- cfg$genes$corr_genes_x %||% NULL
  cfg$genes$corr_genes_y <- cfg$genes$corr_genes_y %||% NULL

# ── Parallelization Defaults ──────────────────────────────────────────
  cfg$parallel$enable           <- cfg$parallel$enable           %||% FALSE
  cfg$parallel$workers          <- cfg$parallel$workers          %||% 4L
  cfg$parallel$strategy         <- cfg$parallel$strategy         %||% "multisession"
  cfg$parallel$max_size_gb      <- cfg$parallel$max_size_gb      %||% 8.0

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

  if (is.character(st) && length(st) == 1L && file.exists(st)) {
    df <- utils::read.delim(st, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  } 

  else if (is.list(st)) {
    df <- do.call(dplyr::bind_rows, lapply(st, as.data.frame, stringsAsFactors = FALSE))
  } 
  else {
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

  return(df)
}