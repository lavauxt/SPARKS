#' Run the complete scRNA-seq pipeline
#' @param base_config_path Path to the template YAML
#' @param override_config_path Path to your custom override YAML
#' @param sample_metadata Data frame
#' @export
sparks <- function(base_config_path, override_config_path = NULL, sample_metadata = NULL) {
  cfg <- load_pipeline_config(base_config_path, override_config_path)
  restore_rng <- .set_seed_preserving_state(cfg$processing$seed)
  on.exit(restore_rng(), add = TRUE)

  make_dir(cfg$pipeline$results_dir)
  initial_result_files <- .snapshot_result_files(cfg$pipeline$results_dir)
  run_log <- file(file.path(cfg$pipeline$results_dir, "pipeline.log"), open = "wt")
  run_warnings <- character(0)
  log_condition <- function(type, condition) {
    writeLines(
      paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " [", type, "] ",
             conditionMessage(condition)),
      run_log
    )
    flush(run_log)
  }
  unlink(file.path(cfg$pipeline$results_dir, "warnings.txt"))
  on.exit({
    if (length(run_warnings) > 0L) {
      writeLines(unique(run_warnings), file.path(cfg$pipeline$results_dir, "warnings.txt"))
    }
    close(run_log)
  }, add = TRUE)

  withCallingHandlers({

  if (!requireNamespace("presto", quietly = TRUE)) {
    message("   [INFO] Optional package 'presto' not installed: Seurat's FindMarkers() / ",
            "FindAllMarkers() fall back to a much slower Wilcoxon test.",
            "\n          Install it with: devtools::install_github('immunogenomics/presto')")
  }

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

  if (isTRUE(cfg$pipeline$independent_samples)) {
    message("   [INFO] pipeline$independent_samples = TRUE: every sample is analysed ",
            "on its own (no merging, no between-sample comparison).")
    sample_metadata$comparison_group <-
      if (anyDuplicated(sample_metadata$folder_id) == 0L) {
        sample_metadata$folder_id
      } else {
        paste(sample_metadata$comparison_group, sample_metadata$folder_id, sep = "_")
      }
  }

  cfg$sample_metadata <- sample_metadata

  for (comp_group in unique(sample_metadata$comparison_group)) {
    message("\n##################################")
    message("  Comparison Group: ", comp_group)
    message("##################################\n")

    group_meta <- sample_metadata[sample_metadata$comparison_group == comp_group, ]
    dirs       <- .setup_group_dirs(cfg$pipeline$results_dir, comp_group)

    protocol_objects <- lapply(seq_len(nrow(group_meta)), function(i) {
      # BUG FIX: an entirely blank quant_format column is read as logical NA,
      # and nzchar(NA_character_) is TRUE -- so NA used to be passed on as the
      # format and every sample failed with "Unknown quant_format 'NA'". A blank
      # (or NA) cell must mean "use input$format".
      qf_cell <- if ("quant_format" %in% colnames(group_meta)) {
        trimws(as.character(group_meta$quant_format[i]))
      } else {
        NA_character_
      }
      qf <- if (!is.na(qf_cell) && nzchar(qf_cell)) qf_cell else NULL
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
          cfg                  = cfg,
          raw_matrix_dir       = dirs$raw_matrix,
          quant_format         = qf
        ),
        label = paste0("process_single_sample: ", group_meta$folder_id[i])
      )
    })

    valid_objects <- Filter(Negate(is.null), protocol_objects)
    if (length(valid_objects) == 0L) {
      message("[SKIP] No valid samples for group: ", comp_group)
      next
    }

    valid_idx  <- !sapply(protocol_objects, is.null)
    merged_obj <- .merge_samples(valid_objects,
                                 folder_ids = group_meta$folder_id[valid_idx])
    cond_col    <- cfg$processing$condition_col
    cond_vec    <- as.character(merged_obj@meta.data[[cond_col]])
    cond_levels <- unique(as.character(group_meta$protocol))
    merged_obj[[cond_col]] <- factor(cond_vec, levels = cond_levels)

    # BUG FIX: this used to alias `sample` to the condition/protocol vector
    # (factor(cond_vec, ...)), i.e. only as many distinct values as there are
    # conditions (typically 2, e.g. WT/KO). save_pseudobulk_counts() counts
    # unique `sample` values against `min_replicates` (default 3) to decide
    # whether to pseudobulk -- with `sample` == condition that count was
    # always the number of *conditions*, never the number of biological
    # replicates, so pseudobulk export was silently skipped on effectively
    # every realistic dataset. `orig.ident` (set in process_single_sample()
    # to protocol_folderid, unique per replicate) is what this should hold.
    merged_obj$sample <- factor(as.character(merged_obj$orig.ident))
    save_cell_counts(merged_obj, paste0("before_SCT_", comp_group), dirs$qc)


    message("   [INFO] Running NormalizeData on RNA assay...")
    merged_obj <- Seurat::NormalizeData(merged_obj, assay = "RNA", verbose = FALSE)

    regress_vars <- cfg$processing$vars_to_regress

    if (cfg$sex_scoring$run || cfg$cell_cycle$run) {

      if (cfg$sex_scoring$run) {
        merged_obj <- run_sex_scoring(
          seurat_obj   = merged_obj,
          config_block = cfg$sex_scoring,
          assay        = "RNA",
          seed         = cfg$processing$seed
        )
        if (cfg$sex_scoring$regress && "Sex.Difference" %in% colnames(merged_obj[[]])) {
          message("   [INFO] Adding Sex.Difference to SCTransform regression.")
          regress_vars <- unique(c(regress_vars, "Sex.Difference"))
        }
      }

      if (cfg$cell_cycle$run) {
        message("   [INFO] Running Cell Cycle Scoring...")

        s.genes   <- Seurat::cc.genes.updated.2019$s.genes
        g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

        if (tolower(cfg$pipeline$species_target) == "mouse") {
          s.genes   <- stringr::str_to_title(tolower(s.genes))
          g2m.genes <- stringr::str_to_title(tolower(g2m.genes))
        }

        available_genes <- rownames(merged_obj[["RNA"]])
        s.genes <- intersect(s.genes, available_genes)
        g2m.genes <- intersect(g2m.genes, available_genes)
        if (length(s.genes) > 0L && length(g2m.genes) > 0L) {
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
        } else {
          message("   [SKIP Cell Cycle] S or G2M marker genes are absent from the RNA assay.")
        }
      }
    }

    merged_obj <- run_seurat_processing(
      seurat_obj      = merged_obj,
      dims_pca        = cfg$processing$pca_dims,
      resolution      = cfg$processing$cluster_resolution,
      npcs            = cfg$processing$npcs,
      seed            = cfg$processing$seed,
      vars_to_regress = regress_vars,
      split_by        = "orig.ident",
      prep_sct_findmarkers = isTRUE(cfg$processing$prep_sct_findmarkers)
    )

    save_png(
      Seurat::ElbowPlot(merged_obj,
        ndims = min(cfg$processing$n_elbow_dims, ncol(merged_obj) - 1L)),
      file.path(dirs$qc, paste0("ElbowPlot_", comp_group, ".png")),
      width  = cfg$plot$elbow_width,
      height = cfg$plot$elbow_height
    )
    save_cell_counts(merged_obj, paste0("after_clustering_", comp_group), dirs$qc)

    merged_obj <- run_singler_annotation(
      seurat_obj     = merged_obj,
      species_target = cfg$pipeline$species_target,
      ref_celldex1   = cfg$species$ref_primary,
      ref_celldex2   = cfg$species$ref_secondary,
      singler_cfg    = cfg$singler
    )

    if (isTRUE(cfg$escape$run)) {
      merged_obj <- run_escape_enrichment(
        seurat_obj = merged_obj,
        species    = cfg$pipeline$species_target,
        library    = cfg$escape$library,
        method     = cfg$escape$method,
        min_size   = cfg$escape$min_size
      )
    }


    groupings_main <- cfg$groupings_main
    if (is.null(groupings_main)) {
      singler_names  <- vapply(cfg$singler$labels, function(x) x$name, character(1))
      groupings_main <- unique(c(cfg$processing$cluster_col, singler_names))
    }

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

  if (requireNamespace("rmarkdown", quietly = TRUE)) {
    qc_template <- .find_report_template(
      cfg$report$rmd_template, cfg$pipeline$config_dir, "qc_report.Rmd"
    )
    results_template <- .find_report_template(
      cfg$report$results_rmd_template, cfg$pipeline$config_dir,
      "results_report.Rmd"
    )
    if (is.null(qc_template)) message("   [SKIP QC report] Template is not available.")
    if (is.null(results_template)) message("   [SKIP Results report] Template is not available.")

    for (comp_group in unique(sample_metadata$comparison_group)) {
      dirs <- .setup_group_dirs(cfg$pipeline$results_dir, comp_group)
      object_path <- file.path(dirs$rdata, paste0("Merged_", comp_group, ".rds"))
      if (!file.exists(object_path)) next
      report_obj <- readRDS(object_path)
      groupings <- cfg$groupings_main
      if (is.null(groupings)) {
        groupings <- unique(c(cfg$processing$cluster_col, cfg$singler$label_names))
      }
      current_files <- .changed_result_files(initial_result_files, dirs$base)

      if (!is.null(qc_template)) {
        safe_run(generate_qc_report(
          seurat_obj = report_obj,
          comp_group = comp_group,
          out_dir = dirs$qc,
          author = cfg$report$author %||% "Pipeline User",
          title = cfg$report$title %||% paste("QC Report -", comp_group),
          rmd_template = qc_template,
          cfg = cfg,
          run_files = intersect(current_files,
                                list.files(dirs$qc, recursive = TRUE,
                                           full.names = TRUE)),
          log_path = file.path(cfg$pipeline$results_dir, "pipeline.log")
        ), label = "QC report")
      }
      if (!is.null(results_template)) {
        report_groupings <- groupings
        if (is.null(report_groupings) || length(report_groupings) == 0L) {
          report_groupings <- unique(c(cfg$processing$cluster_col,
                                       cfg$singler$label_names))
        }
        report_groupings <- report_groupings[
          report_groupings %in% colnames(report_obj@meta.data)
        ]
        safe_run(generate_results_report(
          seurat_obj = report_obj,
          comp_group = comp_group,
          out_dir = dirs$base,
          groupings = report_groupings,
          author = cfg$report$author %||% "Pipeline User",
          title = cfg$report$title %||% paste("Results Report -", comp_group),
          rmd_template = results_template,
          cfg = cfg,
          run_files = current_files,
          log_path = file.path(cfg$pipeline$results_dir, "pipeline.log")
        ), label = "Results report")
      }
      rm(report_obj)
    }
  } else {
    message("   [SKIP HTML reports] Package 'rmarkdown' is not installed.")
  }

  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(cfg$pipeline$results_dir, "session_info.txt"))

  if (length(run_warnings) > 0L) {
    message("Pipeline finished with warnings. See warnings.txt")
  } else {
    unlink(file.path(cfg$pipeline$results_dir, "warnings.txt"))
    message("Pipeline finished successfully.")
  }
  },
  message = function(m) log_condition("MESSAGE", m),
  warning = function(w) {
    run_warnings <<- c(run_warnings, conditionMessage(w))
    log_condition("WARNING", w)
  },
  error = function(e) log_condition("ERROR", e)
  )
}

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

  sct_available <- "SCT" %in% names(seurat_obj@assays)
  if (!sct_available) {
    message("   [SKIP DEG] SCT assay is not available; descriptive processing continues.")
  } else {
    Seurat::DefaultAssay(seurat_obj) <- "SCT"
  }

  gsea_pathways <- NULL
  if (isTRUE(cfg$gsea$run) && requireNamespace("edgeR", quietly = TRUE) &&
      requireNamespace("fgsea", quietly = TRUE)) {
    gsea_pathways <- tryCatch(
      .hallmark_gene_sets(cfg$pipeline$species_target, cfg$gsea$collection),
      error = function(e) {
        message("   [SKIP GSEA] Could not load MSigDB collection '",
                cfg$gsea$collection, "' gene sets: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(gsea_pathways)) {
      message("   [SKIP GSEA] No gene sets loaded for MSigDB collection '",
              cfg$gsea$collection, "'.")
    }
  } else if (isTRUE(cfg$gsea$run) && !requireNamespace("edgeR", quietly = TRUE)) {
    message("   [SKIP GSEA] Ranked GSEA requires Bioconductor package 'edgeR'.")
  } else if (isTRUE(cfg$gsea$run) && !requireNamespace("fgsea", quietly = TRUE)) {
    message("   [SKIP GSEA] Ranked GSEA requires package 'fgsea'.")
  }

  # Older or externally generated objects may still retain an SCTAssay.
  prep_status <- seurat_obj@misc$sparks_prep_sct_findmarkers %||% "unknown"
  if (isTRUE(cfg$processing$prep_sct_findmarkers) &&
      sct_available &&
      identical(prep_status, "unknown") &&
      inherits(seurat_obj[["SCT"]], "SCTAssay")) {
    prep_result <- tryCatch(
      list(object = Seurat::PrepSCTFindMarkers(seurat_obj, verbose = FALSE),
           status = "succeeded_on_analysis_object"),
      error = function(e) {
        message("   [WARNING] PrepSCTFindMarkers failed: ", conditionMessage(e))
        list(object = seurat_obj, status = "failed")
      }
    )
    seurat_obj <- prep_result$object
    prep_status <- prep_result$status
  }
  if (inherits(seurat_obj[["SCT"]], "Assay5")) {
    message("   [INFO] SCT Assay5 cannot retain SCT models; PrepSCTFindMarkers status: ",
            prep_status, ".")
  } else if (!isTRUE(cfg$processing$prep_sct_findmarkers)) {
    message("   [INFO] PrepSCTFindMarkers disabled by configuration.")
  }
  skip_cell_deg <- !sct_available || identical(prep_status, "failed") ||
    (isTRUE(cfg$processing$prep_sct_findmarkers) &&
       identical(prep_status, "not_applicable_no_sct_models"))
  if (skip_cell_deg) {
    message("   [SKIP DEG] SCT preparation failed; cell-level SCT differential expression is unavailable.")
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
      cfg          = cfg,
      gsea_pathways = gsea_pathways,
      skip_cell_deg = skip_cell_deg
    )
  }

  fp_dir <- file.path(base_dir, suffix, display_name, "FeaturePlot")
  make_dir(fp_dir)
  generate_feature_plots(seurat_obj, genes_list, fp_dir, display_name,
                         reduction = cfg$processing$reduction)
}

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
                                   deg_color, cfg, gsea_pathways = NULL,
                                   skip_cell_deg = FALSE) {
  message("   -> Grouping: ", group_col)

  group_dir <- file.path(base_dir, suffix, file_prefix, group_col)
  dirs      <- .make_analysis_dirs(group_dir)

  singler_entry <- Filter(function(l) l$name == group_col, cfg$singler$labels)
  is_fine       <- length(singler_entry) > 0L && isTRUE(singler_entry[[1L]]$is_fine)
  umap_w        <- if (is_fine) cfg$plot$umap_width_fine  else cfg$plot$umap_width_standard
  umap_h        <- if (is_fine) cfg$plot$umap_height_fine else cfg$plot$umap_height_standard

  do_deg <- isTRUE(cfg$deg$run %||% TRUE)
  do_cell_deg <- do_deg && !isTRUE(skip_cell_deg)
  if (!do_deg) {
    message("   [INFO] deg$run = FALSE: skipping between-condition DEG and statistics.")
  } else if (!do_cell_deg) {
    message("   [SKIP DEG] SCT preparation/assay unavailable; cell-level DEG skipped.")
  }

  cond_vals <- unique(seurat_obj@meta.data[[cfg$processing$condition_col]])
  p_umap <- Seurat::DimPlot(seurat_obj,
    reduction = cfg$processing$reduction,
    group.by  = group_col,
    label     = TRUE, repel = TRUE,
    split.by  = cfg$processing$condition_col) +
    ggplot2::ggtitle(paste0(file_prefix, " | ", group_col, " | ",
                             paste(cond_vals, collapse = if (do_cell_deg) " vs " else " / "),
                             " | ", suffix))

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

  generate_violin_plots(seurat_obj, genes_list, dirs$VlnPlot, file_prefix,
                         group_by_col = group_col)

  if (isTRUE(cfg$proportions$run)) {
    run_proportion_analysis(seurat_obj, group_col, dirs$Proportions, file_prefix,
                            run_test = do_cell_deg)
    if (do_cell_deg) {
      run_scproportion_test(seurat_obj, group_col, dirs$Proportions, file_prefix)
    }
  }

  if (do_cell_deg) {
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
  }

  if (isTRUE(cfg$pseudobulk$run) || isTRUE(cfg$gsea$run)) {
    .run_sample_level_analysis(
      seurat_obj = seurat_obj,
      group_col = group_col,
      pseudobulk_dir = file.path(dirs$Pseudobulk, "DEG"),
      gsea_dir = dirs$GSEA,
      file_prefix = file_prefix,
      species = cfg$pipeline$species_target,
      condition_col = cfg$processing$condition_col,
      sample_col = "sample",
      min_replicates = cfg$pseudobulk$min_replicates,
      run_pseudobulk_deg = isTRUE(cfg$pseudobulk$run),
      run_gsea = isTRUE(cfg$gsea$run),
      pathways = gsea_pathways,
      collection = cfg$gsea$collection,
      gsea_min_size = cfg$gsea$min_size,
      gsea_max_size = cfg$gsea$max_size
    )
  }

  for (layer in cfg$expression$layers) {
    save_average_expression(seurat_obj, dirs$Expression, file_prefix,
      group_by_col    = group_col,
      layer           = layer,
      table_sep       = cfg$deg$table_sep,
      table_quote     = cfg$deg$table_quote,
      table_row_names = cfg$deg$table_row_names)
  }

  if (isTRUE(cfg$pseudobulk$save_counts)) {
    save_pseudobulk_counts(seurat_obj, dirs$Pseudobulk, file_prefix,
                           group_by_col  = group_col,
                           table_sep     = cfg$deg$table_sep,
                           condition_col = cfg$processing$condition_col,
                           min_replicates = cfg$pseudobulk$min_replicates)
  }

  generate_cluster_markers_and_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix)

  generate_cluster_zscore_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix,
                                  species_target = cfg$pipeline$species_target,
                                  top_n          = cfg$plot$top_genes_heatmap_n)

  generate_cluster_zscore_heatmap_split_condition(
    seurat_obj,
    group_by_col    = group_col,
    condition_col   = cfg$processing$condition_col,
    out_dir         = dirs$Heatmap,
    prefix          = file_prefix,
    species_target  = cfg$pipeline$species_target,
    top_n           = cfg$plot$top_genes_heatmap_n
  )

  generate_expression_heatmap(seurat_obj, group_col, dirs$Heatmap, file_prefix,
                               species_target = cfg$pipeline$species_target)
  generate_top_expressed_genes(seurat_obj, group_col, dirs$Heatmap, file_prefix)

  for (sig in cfg$gene_signatures) {
    scale_dot <- sig$dot_scale %||% FALSE
    safe_run(
      generate_gene_signature_plots(
        seurat_obj     = seurat_obj,
        genes          = sig$genes,
        out_dir        = dirs$Heatmap,
        prefix         = file_prefix,
        group_by_col   = group_col,
        signature_name = sig$name,
        condition_col  = cfg$processing$condition_col,
        scale_dotplot  = scale_dot          
      ),
      label = paste0("Gene Signature '", sig$name, "': ", file_prefix, " | ", group_col)
    )
  }

  for (sig in cfg$gene_signatures) {
    scale_dot <- sig$dot_scale %||% FALSE  
    safe_run(
      generate_gene_signature_per_group_dotplots(
        seurat_obj          = seurat_obj,
        genes               = sig$genes,
        out_dir             = dirs$Heatmap,
        prefix              = file_prefix,
        group_by_col        = group_col,
        signature_name      = sig$name,
        condition_col       = cfg$processing$condition_col,
        min_cells_per_group = cfg$labeling$min_subset_cells %||% 10L,
        scale_dotplot       = scale_dot         
      ),
      label = paste0("Per-group DotPlot '", sig$name, "': ", file_prefix, " | ", group_col)
    )

      safe_run(
          generate_gene_signature_boxplots(
            seurat_obj          = seurat_obj,
            genes               = sig$genes,
            out_dir             = dirs$Heatmap,
            prefix              = file_prefix,
            group_by_col        = group_col,
            signature_name      = sig$name,
            condition_col       = cfg$processing$condition_col,
            min_cells_per_group = cfg$labeling$min_subset_cells %||% 10L,
            add_jitter          = TRUE   # set to FALSE if you want only boxes
          ),
          label = paste0("Boxplot '", sig$name, "': ", file_prefix, " | ", group_col)
        )
      }

  if (isTRUE(cfg$escape$run)) {
    generate_escape_plots(seurat_obj, method = cfg$escape$method,
                          group_col = group_col, out_dir = group_dir,
                          prefix = file_prefix)
  } else {
    message("   [INFO] Per-cell escape scoring disabled by configuration.")
  }

  cx <- cfg$genes$corr_genes_x
  cy <- cfg$genes$corr_genes_y
  if (!is.null(cx) && !is.null(cy) && length(cx) > 0L && length(cy) > 0L) {
    run_gene_correlations(
      seurat_obj,
      grouping_col = group_col,
      genes_x      = cx,
      genes_y      = cy,
      out_dir      = dirs$Correlation,
      prefix       = paste0(file_prefix, "_", group_col),
      cond_col     = "condition",
      method       = "pearson",
      assay        = "RNA"
    )
  } else {
    message("[SKIP Correlation] Genes to correlate not set")
  }
}

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
  available_subset_layers <- tryCatch(
    SeuratObject::Layers(sub_obj[["RNA"]]),
    error = function(e) character(0)
  )
  if (!"data" %in% available_subset_layers && "counts" %in% available_subset_layers) {
    message("   [INFO] RNA 'data' layer absent in subset; running NormalizeData.")
    sub_obj <- Seurat::NormalizeData(sub_obj, assay = "RNA", verbose = FALSE)
  }
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
    seed       = cfg$processing$seed,
    split_by   = "orig.ident",
    prep_sct_findmarkers = isTRUE(cfg$processing$prep_sct_findmarkers))

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

  # `genes: ~` (or no `genes` key) falls back to genes$genes_to_plot, as the
  # config templates document -- the code used to pass NULL through instead.
  sub_genes <- subset_cfg$genes %||% cfg$genes$genes_to_plot

  run_analysis_unit(seurat_obj = sub_obj, display_name = subset_cfg$display_name,
                    groupings  = groupings, genes_list = sub_genes,
                    base_dir   = base_dir, suffix = suffix,
                    deg_color  = subset_cfg$deg_color, cfg = cfg)

  if (length(sub_genes) > 0) {
    message("   --> Generating Subcluster Heatmaps for: ", subset_cfg$display_name)
    hm_dir <- file.path(base_dir, suffix, subset_cfg$display_name, "Heatmap")
    generate_subcluster_heatmaps(
      seurat_obj = sub_obj,
      genes      = sub_genes,
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
