# ──────────────────────────────────────────────────────────────────────────────
#' Run the complete scRNA-seq pipeline
#' @param base_config_path Path to the template YAML
#' @param override_config_path Path to your custom override YAML
#' @param sample_metadata Data frame
#' @export
sparks <- function(base_config_path, override_config_path = NULL, sample_metadata = NULL) {
  cfg <- load_pipeline_config(base_config_path, override_config_path)

  if (isTRUE(cfg$parallel$enable)) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("The 'future' package is required for parallelization. Please run: install.packages('future')")
    }

    message("   [Parallelization] Setting up ", cfg$parallel$workers,
            " workers using strategy: ", cfg$parallel$strategy)

    future::plan(strategy = cfg$parallel$strategy, workers = cfg$parallel$workers)
    options(future.globals.maxSize = cfg$parallel$max_size_gb * 1024^3)

    on.exit({
      message("   [Parallelization] Reverting to sequential execution and closing workers...")
      future::plan("sequential")
    }, add = TRUE)
  }

  if (is.null(sample_metadata)) {
    message("   [INFO] No sample_metadata argument provided. Loading from YAML config...")
    sample_metadata <- load_sample_table(cfg)
  }

  cfg$sample_metadata <- sample_metadata
  make_dir(cfg$pipeline$results_dir)

  for (comp_group in unique(sample_metadata$comparison_group)) {
    message("\n##################################")
    message("  Comparison Group: ", comp_group)
    message("##################################\n")

    group_meta <- sample_metadata[sample_metadata$comparison_group == comp_group, ]
    dirs       <- .setup_group_dirs(cfg$pipeline$results_dir, comp_group)

    # ── Per-sample processing ──────────────────────────────────────────────────
    protocol_objects <- lapply(seq_len(nrow(group_meta)), function(i) {
      safe_run(
        process_single_sample(
          folder_id            = group_meta$folder_id[i],
          protocol             = group_meta$protocol[i],
          file_prefix          = paste0(group_meta$protocol[i], "_",
                                        group_meta$folder_id[i]),
          qc_dir               = dirs$qc,
          data_path            = cfg$pipeline$data_dir,
          gene_removal_pattern = cfg$species$gene_removal_pattern,
          mt_pattern           = cfg$species$mt_pattern,
          genes_to_remove      = cfg$species$genes_to_remove,
          min_features         = cfg$qc$min_features,
          max_features         = cfg$qc$max_features,
          max_counts           = cfg$qc$max_counts,
          max_mt_percent       = cfg$qc$max_mt_percent,
          min_cells            = cfg$qc$min_cells,
          cfg                  = cfg
        ),
        label = paste0("process_single_sample: ", group_meta$folder_id[i])
      )
    })

    valid_objects <- Filter(Negate(is.null), protocol_objects)
    if (length(valid_objects) == 0L) {
      message("[SKIP] No valid samples for group: ", comp_group)
      next
    }

    # ── Merge ──────────────────────────────────────────────────────────────────
    valid_idx  <- !sapply(protocol_objects, is.null)
    merged_obj <- .merge_samples(valid_objects,
                                 folder_ids = group_meta$folder_id[valid_idx])
    cond_col    <- cfg$processing$condition_col
    cond_vec    <- as.character(merged_obj@meta.data[[cond_col]])
    cond_levels <- unique(as.character(group_meta$protocol))
    merged_obj[[cond_col]] <- factor(cond_vec, levels = cond_levels)
    # Keep orig.ident as per-sample IDs so Harmony corrects across individual
    # samples, not just conditions.  'sample' is a convenience condition alias.
    merged_obj$sample <- factor(cond_vec, levels = cond_levels)
    save_cell_counts(merged_obj, paste0("before_SCT_", comp_group), dirs$qc)

    # BUG FIX #5: ensure RNA assay has log-normalised data exactly once.
    # The original code had a conditional NormalizeData guard here AND an
    # unconditional NormalizeData call inside the sex/cell-cycle block,
    # causing a redundant double normalisation whenever sex or cell-cycle
    # scoring was enabled.  Run it once unconditionally here; the block below
    # no longer repeats it.
    message("   [INFO] Running NormalizeData on RNA assay...")
    merged_obj <- Seurat::NormalizeData(merged_obj, assay = "RNA", verbose = FALSE)

    # ── Regression variables ────────────────────────────────────────────────────
    regress_vars <- cfg$processing$vars_to_regress

    # ── Sex & Cell Cycle scoring ───────────────────────────────────────────────
    if (cfg$sex_scoring$run || cfg$cell_cycle$run) {

      # --- Sex Scoring ---
      if (cfg$sex_scoring$run) {
        merged_obj <- run_sex_scoring(
          seurat_obj   = merged_obj,
          config_block = cfg$sex_scoring,
          assay        = "RNA"
        )
        if (cfg$sex_scoring$regress && "Sex.Difference" %in% colnames(merged_obj[[]])) {
          message("   [INFO] Adding Sex.Difference to SCTransform regression.")
          regress_vars <- unique(c(regress_vars, "Sex.Difference"))
        }
      }

      # --- Cell Cycle Scoring ---
      if (cfg$cell_cycle$run) {
        message("   [INFO] Running Cell Cycle Scoring...")

        s.genes   <- Seurat::cc.genes.updated.2019$s.genes
        g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

        # BUG FIX #6: original code compared species_target to lowercase
        # "mouse" but the YAML stores "Mouse" (capital M), so this branch was
        # never entered and mouse cell-cycle genes were never title-cased.
        if (tolower(cfg$pipeline$species_target) == "mouse") {
          s.genes   <- stringr::str_to_title(tolower(s.genes))
          g2m.genes <- stringr::str_to_title(tolower(g2m.genes))
        }

        merged_obj <- Seurat::CellCycleScoring(
          object       = merged_obj,
          s.features   = s.genes,
          g2m.features = g2m.genes,
          assay        = "RNA"
        )

        if (cfg$cell_cycle$regress) {
          message("   [INFO] Adding S.Score, G2M.Score to SCTransform regression.")
          regress_vars <- unique(c(regress_vars, "S.Score", "G2M.Score"))
        }
      }
    }

    # ── Seurat processing ──────────────────────────────────────────────────────
    merged_obj <- run_seurat_processing(
      seurat_obj      = merged_obj,
      dims_pca        = cfg$processing$pca_dims,
      resolution      = cfg$processing$cluster_resolution,
      npcs            = cfg$processing$npcs,
      vars_to_regress = regress_vars,
      split_by        = "orig.ident"
    )

    # ── Elbow Plot ─────────────────────────────────────────────────────────────
    save_png(
      Seurat::ElbowPlot(merged_obj,
        ndims = min(cfg$processing$n_elbow_dims, ncol(merged_obj) - 1L)),
      file.path(dirs$qc, paste0("ElbowPlot_", comp_group, ".png")),
      width  = cfg$plot$elbow_width,
      height = cfg$plot$elbow_height
    )
    save_cell_counts(merged_obj, paste0("after_clustering_", comp_group), dirs$qc)

    # ── SingleR ────────────────────────────────────────────────────────────────
    merged_obj <- run_singler_annotation(
      seurat_obj     = merged_obj,
      species_target = cfg$pipeline$species_target,
      ref_celldex1   = cfg$species$ref_primary,
      ref_celldex2   = cfg$species$ref_secondary,
      singler_cfg    = cfg$singler
    )

    # ── Escape Enrichment ──────────────────────────────────────────────────────
    if (isTRUE(cfg$escape$run)) {
      merged_obj <- run_escape_enrichment(
        seurat_obj = merged_obj,
        species    = cfg$pipeline$species_target,
        library    = cfg$escape$library,
        method     = cfg$escape$method,
        min_size   = cfg$escape$min_size
      )
    }

    # ── QC HTML Report ─────────────────────────────────────────────────────────
    # BUG FIX #7: template lookup was using a fragile relative path
    # (results_dir/../templates/qc_report.Rmd) that is almost never correct.
    # Now we search in the config file's own directory first, then in getwd(),
    # then honour an explicit cfg$report$rmd_template override.
    if (requireNamespace("rmarkdown", quietly = TRUE)) {
      template_path <- cfg$report$rmd_template %||% {
        candidates <- c(
          file.path(cfg$pipeline$config_dir, "qc_report.Rmd"),
          file.path(getwd(),                 "qc_report.Rmd")
        )
        found <- candidates[file.exists(candidates)]
        if (length(found) > 0L) found[1L] else NULL
      }

      if (is.null(template_path)) {
        message("   [WARNING] qc_report.Rmd not found – skipping HTML report.",
                "\n            Place qc_report.Rmd alongside your config YAML or set",
                "\n            cfg$report$rmd_template in your override config.")
      } else {
        generate_qc_report(
          seurat_obj   = merged_obj,
          comp_group   = comp_group,
          out_dir      = dirs$qc,          # HTML saved inside QC folder
          author       = cfg$report$author %||% "Pipeline User",
          title        = cfg$report$title  %||% paste("QC Report -", comp_group),
          rmd_template = template_path,
          cfg          = cfg
        )
      }
    } else {
      message("   [SKIP] rmarkdown not installed – QC report not generated.")
    }

    # ── Main groupings ─────────────────────────────────────────────────────────
    groupings_main <- cfg$groupings_main
    if (is.null(groupings_main)) {
      singler_names  <- vapply(cfg$singler$labels, function(x) x$name, character(1))
      groupings_main <- unique(c(cfg$processing$cluster_col, singler_names))
    }

    # ── Main analysis unit ─────────────────────────────────────────────────────
    run_analysis_unit(
      seurat_obj   = merged_obj,
      display_name = "Main",
      groupings    = groupings_main,
      genes_list   = cfg$genes$genes_to_plot,
      base_dir     = cfg$pipeline$results_dir,
      suffix       = comp_group,
      deg_color    = cfg$plot$deg_color_main,
      cfg          = cfg
    )

    # ── Subsets ────────────────────────────────────────────────────────────────
    for (s in cfg$subsets) {
      safe_run(
        .run_subset(main_obj = merged_obj, subset_cfg = s,
                    base_dir = cfg$pipeline$results_dir,
                    suffix   = comp_group, cfg = cfg),
        label = paste0("subset: ", s$display_name)
      )
    }

    .save_rdata(merged_obj, dirs$rdata, paste0("Merged_", comp_group))
    rm(merged_obj); gc()
  }

  writeLines(capture.output(utils::sessionInfo()),
             file.path(cfg$pipeline$results_dir, "session_info.txt"))

  w <- capture.output(warnings())
  if (length(w) > 0L && !all(grepl("no warnings", w, ignore.case = TRUE))) {
    writeLines(w, file.path(cfg$pipeline$results_dir, "warnings.txt"))
    message("Pipeline finished with warnings. See warnings.txt")
  } else {
    message("Pipeline finished successfully.")
  }
}

# ──────────────────────────────────────────────────────────────────────────────
#' Run all analyses for a Seurat object (main or subset)
#' @param seurat_obj Seurat object
#' @param display_name Character
#' @param groupings Character vector. Metadata columns to loop over
#' @param genes_list Character vector
#' @param base_dir Character
#' @param suffix Character
#' @param deg_color Character
#' @param cfg Named list
#' @return NULL
#' @export
run_analysis_unit <- function(seurat_obj, display_name, groupings, genes_list,
                               base_dir, suffix, deg_color, cfg) {
  message("\n=== Analysis Unit: ", display_name, " ===")

  valid_groupings <- groupings[groupings %in% colnames(seurat_obj@meta.data)]
  skipped         <- setdiff(groupings, valid_groupings)
  if (length(skipped) > 0L)
    message("   [SKIP] columns not found: ", paste(skipped, collapse = ", "))

  Seurat::DefaultAssay(seurat_obj) <- "SCT"

  if (inherits(seurat_obj[["SCT"]], "Assay5")) {
    message("   [INFO] Assay5 detected. Skipping PrepSCTFindMarkers (not required).")
  } else {
    seurat_obj <- safe_run(
      Seurat::PrepSCTFindMarkers(seurat_obj, verbose = TRUE),
      label    = "PrepSCTFindMarkers",
      fallback = seurat_obj
    )
  }

  for (grp in valid_groupings) {
    run_grouping_analysis(
      seurat_obj   = seurat_obj,
      group_col    = grp,
      file_prefix  = display_name,
      genes_list   = genes_list,
      base_dir     = base_dir,
      suffix       = suffix,
      deg_color    = deg_color,
      cfg          = cfg
    )
  }

  fp_dir <- file.path(base_dir, suffix, display_name, "FeaturePlot")
  make_dir(fp_dir)
  generate_feature_plots(seurat_obj, genes_list, fp_dir, display_name,
                         reduction = cfg$processing$reduction)
}

# ──────────────────────────────────────────────────────────────────────────────
#' Run all analyses for one grouping column
#' @param seurat_obj Seurat object. Must already have PrepSCTFindMarkers applied.
#' @param group_col Character
#' @param file_prefix Character
#' @param genes_list Character vector
#' @param base_dir Character
#' @param suffix Character
#' @param deg_color Character
#' @param cfg Named list
#' @return NULL
#' @export
run_grouping_analysis <- function(seurat_obj, group_col, file_prefix,
                                   genes_list, base_dir, suffix,
                                   deg_color, cfg) {
  message("   -> Grouping: ", group_col)

  group_dir <- file.path(base_dir, suffix, file_prefix, group_col)
  dirs      <- .make_analysis_dirs(group_dir)

  # ── UMAP ────────────────────────────────────────────────────────────────────
  singler_entry <- Filter(function(l) l$name == group_col, cfg$singler$labels)
  is_fine       <- length(singler_entry) > 0L && isTRUE(singler_entry[[1L]]$is_fine)
  umap_w        <- if (is_fine) cfg$plot$umap_width_fine  else cfg$plot$umap_width_standard
  umap_h        <- if (is_fine) cfg$plot$umap_height_fine else cfg$plot$umap_height_standard

  cond_vals <- unique(seurat_obj@meta.data[[cfg$processing$condition_col]])
  p_umap <- Seurat::DimPlot(seurat_obj,
    reduction = cfg$processing$reduction,
    group.by  = group_col,
    label     = TRUE, repel = TRUE,
    split.by  = cfg$processing$condition_col) +
    ggplot2::ggtitle(paste0(file_prefix, " | ", group_col, " | ",
                             paste(cond_vals, collapse = " vs "), " | ", suffix))

  if (is_fine)
    p_umap <- p_umap +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = ggplot2::guide_legend(
        nrow         = cfg$plot$legend_nrow_fine,
        override.aes = list(size = 3)
      ))

  save_png(p_umap,
    file.path(dirs$UMAP, paste0("UMAP_", file_prefix, "_", group_col, ".png")),
    width = umap_w, height = umap_h)

  # ── Violin plots ────────────────────────────────────────────────────────────
  generate_violin_plots(seurat_obj, genes_list, dirs$VlnPlot, file_prefix,
                         group_by_col = group_col)

  # ── Proportions ─────────────────────────────────────────────────────────────
  run_proportion_analysis(seurat_obj, group_col, dirs$DEG, file_prefix)
  run_scproportion_test(seurat_obj,   group_col, dirs$DEG, file_prefix)

  # ── DEG ─────────────────────────────────────────────────────────────────────
  all_markers <- run_deg_analysis(seurat_obj,
    logfc_threshold     = cfg$deg$logfc_threshold,
    min_pct             = cfg$deg$min_pct,
    group_by_col        = group_col,
    min_cells_per_group = cfg$deg$min_cells_per_group)

  valid_groups <- get_valid_groups(seurat_obj@meta.data, group_col,
                                   min_cells = cfg$deg$min_cells_per_group)
  deg_counts   <- process_and_save_deg(all_markers, valid_groups,
    dirs$DEG, file_prefix,
    group_by_col    = group_col,
    padj_threshold  = cfg$deg$min_p_val_adj,
    table_sep       = cfg$deg$table_sep,
    table_quote     = cfg$deg$table_quote,
    table_row_names = cfg$deg$table_row_names)

  generate_deg_umap(seurat_obj, deg_counts, dirs$UMAP, file_prefix,
    group_by_col    = group_col,
    reduction       = cfg$processing$reduction,
    color_high      = deg_color,
    color_low       = cfg$plot$deg_umap_color_low,
    min_deg_display = cfg$deg$min_deg_display,
    point_size      = cfg$plot$deg_umap_point_size,
    alpha           = cfg$plot$deg_umap_alpha,
    label_size      = cfg$plot$deg_umap_label_size,
    width           = cfg$plot$deg_umap_width,
    height          = cfg$plot$deg_umap_height)

  # ── Average expression ──────────────────────────────────────────────────────
  for (layer in cfg$deg$avg_expression_layers) {
    save_average_expression(seurat_obj, dirs$DEG, file_prefix,
      group_by_col    = group_col,
      layer           = layer,
      table_sep       = cfg$deg$table_sep,
      table_quote     = cfg$deg$table_quote,
      table_row_names = cfg$deg$table_row_names)
  }

  # ── Pseudobulk ──────────────────────────────────────────────────────────────
  save_pseudobulk_counts(seurat_obj, dirs$DEG, file_prefix,
                         group_by_col = group_col,
                         table_sep    = cfg$deg$table_sep)

  # ── Heatmaps ────────────────────────────────────────────────────────────────
  generate_cluster_markers_and_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix)

  generate_cluster_zscore_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix,
                                  species_target = cfg$pipeline$species_target,
                                  top_n          = cfg$plot$top_genes_heatmap_n)

  generate_expression_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix,
                               species_target = cfg$pipeline$species_target)
  generate_top_expressed_genes(seurat_obj, group_col, dirs$Heatmap, file_prefix)

  # ── Escape Enrichment Plots ─────────────────────────────────────────────────
  if (isTRUE(cfg$escape$run)) {
    generate_escape_plots(seurat_obj, method = cfg$escape$method,
                          group_col = group_col, out_dir = group_dir,
                          prefix = file_prefix)
  }

  # ── Correlations ────────────────────────────────────────────────────────────
  cx <- cfg$genes$corr_genes_x
  cy <- cfg$genes$corr_genes_y
  if (!is.null(cx) && !is.null(cy) && length(cx) > 0L && length(cy) > 0L) {
    run_gene_correlations(
      seurat_obj,
      grouping_col = group_col,
      genes_x      = cx,
      genes_y      = cy,
      out_dir      = group_dir,
      prefix       = paste0(file_prefix, "_", group_col),
      cond_col     = "condition",
      method       = "pearson",
      assay        = "RNA"
    )
  } else {
    message("[SKIP Correlation] Genes to correlate not set")
  }
}

# ──────────────────────────────────────────────────────────────────────────────
.run_subset <- function(main_obj, subset_cfg, base_dir, suffix, cfg) {
  message("\n=== Subset: ", subset_cfg$display_name, " ===")

  match_data <- as.character(main_obj@meta.data[[subset_cfg$match_col]])
  idx <- if (isTRUE(subset_cfg$exact_match)) {
    match_data %in% subset_cfg$pattern
  } else {
    grepl(paste(subset_cfg$pattern, collapse = "|"), match_data, ignore.case = TRUE)
  }

  if (sum(idx) < cfg$labeling$min_subset_cells) {
    message("   [SKIP] Only ", sum(idx), " cells (min: ",
            cfg$labeling$min_subset_cells, ").")
    return(invisible(NULL))
  }
  message("   Cells selected: ", sum(idx))

  sub_obj <- subset(main_obj, cells = colnames(main_obj)[idx])
  Seurat::DefaultAssay(sub_obj) <- "RNA"
  sub_obj <- SeuratObject::JoinLayers(sub_obj)

  # BUG FIX #8: the YAML defines per-subset dimensions as pca_dims_from /
  # pca_dims_to (same convention as the top-level processing block), but the
  # original code passed subset_cfg$pca_dims which is always NULL after YAML
  # parsing.  The function therefore always fell back to the default 1:20,
  # ignoring the user-specified subset dimensions entirely.
  pca_dims_sub <- if (!is.null(subset_cfg$pca_dims_from) &&
                       !is.null(subset_cfg$pca_dims_to)) {
    seq.int(subset_cfg$pca_dims_from, subset_cfg$pca_dims_to)
  } else {
    subset_cfg$pca_dims %||% cfg$processing$pca_dims
  }

  sub_obj <- run_seurat_processing(sub_obj,
    dims_pca   = pca_dims_sub,
    resolution = cfg$processing$cluster_resolution,
    npcs       = cfg$processing$npcs,
    split_by   = "orig.ident")

  singler_names <- cfg$singler$label_names[
    cfg$singler$label_names %in% colnames(sub_obj@meta.data)]

  if (length(subset_cfg$label_rules) > 0L) {
    sub_obj <- label_cells_by_markers(sub_obj,
      new_col_name     = subset_cfg$type_col,
      rules            = subset_cfg$label_rules,
      unassigned_label = paste0(subset_cfg$display_name,
                                cfg$labeling$unassigned_suffix),
      threshold        = cfg$labeling$marker_positive_threshold)
    groupings <- unique(c(cfg$processing$cluster_col, singler_names,
                          subset_cfg$type_col))
  } else {
    sub_obj[[subset_cfg$type_col]] <- sub_obj[[cfg$processing$cluster_col]]
    groupings <- unique(c(cfg$processing$cluster_col, singler_names))
  }

  run_analysis_unit(seurat_obj = sub_obj, display_name = subset_cfg$display_name,
                    groupings  = groupings, genes_list = subset_cfg$genes,
                    base_dir   = base_dir, suffix = suffix,
                    deg_color  = subset_cfg$deg_color, cfg = cfg)

  if (!is.null(subset_cfg$genes) && length(subset_cfg$genes) > 0) {
    message("   --> Generating Subcluster Heatmaps for: ", subset_cfg$display_name)
    hm_dir <- file.path(base_dir, suffix, subset_cfg$display_name, "Heatmap")
    generate_subcluster_heatmaps(
      seurat_obj = sub_obj,
      genes      = subset_cfg$genes,
      out_dir    = hm_dir,
      prefix     = subset_cfg$display_name
    )
  }

  rdata_dir <- file.path(base_dir, suffix, "RData")
  .save_rdata(sub_obj, rdata_dir,
    paste0("Subset_",
           gsub("[^A-Za-z0-9]", "_", paste(subset_cfg$pattern, collapse = "_")),
           "_", suffix))
  message("=== Done: ", subset_cfg$display_name, " ===")
}
