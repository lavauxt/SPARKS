#' Run DEG analysis (WT vs KO) per group level
#'
#' Uses SCT assay with PrepSCTFindMarkers + FindMarkers.
#' Returns a combined data.frame (like FindAllMarkers) with a `cluster` column,
#' or NULL if no DEGs found.
#'
#' @param seurat_obj Seurat object. SCT assay must be active + PrepSCTFindMarkers already called.
#' @param logfc_threshold Numeric
#' @param min_pct Numeric
#' @param group_by_col Character. Metadata column for cell groups
#' @param min_cells_per_group Integer
#' @return data.frame or NULL
#' @export
run_deg_analysis <- function(seurat_obj, logfc_threshold, min_pct,
                              group_by_col, min_cells_per_group = 10L) {
  protocols <- unique(as.character(seurat_obj$condition))
  if (length(protocols) < 2L) {
    message("   [SKIP] DEG: only 1 condition present.")
    return(NULL)
  }

  seurat_obj$group_condition <- paste(
    seurat_obj@meta.data[[group_by_col]],
    seurat_obj$condition,
    sep = "_"
  )
  Seurat::Idents(seurat_obj) <- "group_condition"
  all_idents <- as.character(unique(Seurat::Idents(seurat_obj)))

  groups <- get_valid_groups(seurat_obj@meta.data, group_by_col,
                              min_cells = min_cells_per_group)
  if (length(groups) == 0L) {
    message("   [SKIP] DEG: no groups with >= ", min_cells_per_group, " cells.")
    return(NULL)
  }

  cond1 <- protocols[1L]
  cond2 <- protocols[2L]

  marker_list <- lapply(stats::setNames(groups, groups), function(grp) {
    id1 <- paste(grp, cond1, sep = "_")
    id2 <- paste(grp, cond2, sep = "_")

    if (!id1 %in% all_idents || !id2 %in% all_idents) {
      message("   [SKIP] DEG '", grp, "': ident not found.")
      return(NULL)
    }

    n1 <- sum(Seurat::Idents(seurat_obj) == id1)
    n2 <- sum(Seurat::Idents(seurat_obj) == id2)
    if (n1 < min_cells_per_group || n2 < min_cells_per_group) {
      message("   [SKIP] DEG '", grp, "': too few cells (",
              n1, " vs ", n2, ").")
      return(NULL)
    }

    markers <- safe_run(
      Seurat::FindMarkers(seurat_obj,
        ident.1         = id1,
        ident.2         = id2,
        assay           = "SCT",
        logfc.threshold = logfc_threshold,
        min.pct         = min_pct,
        test.use        = "wilcox",
        verbose         = FALSE),
      label = paste0("FindMarkers '", grp, "'")
    )
    if (is.null(markers) || nrow(markers) == 0L) return(NULL)

    markers$gene    <- rownames(markers)
    markers$cluster <- grp
    markers
  })

  result <- do.call(rbind, Filter(Negate(is.null), marker_list))
  if (is.null(result) || nrow(result) == 0L) return(NULL)
  result
}

#' Process and save DEG results to TSV/TXT files
#'
#' Writes:
#'   - AllMarkers_{prefix}_by_{col}.txt  (all results)
#'   - DEG_{prefix}_{group}_by_{col}.txt (per group, significant only)
#'   - DEG_Counts_{prefix}_by_{col}.txt  (summary counts)
#'
#' @param all_markers data.frame from run_deg_analysis or NULL
#' @param group_ids Character vector. All expected group IDs (for zero-fill)
#' @param deg_dir Character
#' @param file_prefix Character
#' @param group_by_col Character
#' @param padj_threshold Numeric
#' @param table_sep Character
#' @param table_quote Logical
#' @param table_row_names Logical
#' @return data.frame of DEG counts per group
#' @export
process_and_save_deg <- function(all_markers, group_ids, deg_dir, file_prefix,
                                  group_by_col, padj_threshold = 0.05,
                                  table_sep       = "\t",
                                  table_quote     = FALSE,
                                  table_row_names = FALSE) {
  make_dir(deg_dir)

  .save <- function(df, stem) {
    safe_run(
      write_table(df,
        file.path(deg_dir,
                  paste0(stem, "_", file_prefix, "_by_", group_by_col, ".txt")),
        sep       = table_sep,
        quote     = table_quote,
        row_names = table_row_names),
      label = paste0("write_table: ", stem)
    )
  }

  empty_counts <- data.frame(
    Group     = as.character(group_ids),
    DEG_Count = 0L,
    stringsAsFactors = FALSE
  )

  if (is.null(all_markers) || nrow(all_markers) == 0L) {
    message("   No DEGs found for '", group_by_col, "'")
    .save(empty_counts, "DEG_Counts")
    return(invisible(empty_counts))
  }

  .save(all_markers, "AllMarkers")

  sig <- all_markers[
    !is.na(all_markers$p_val_adj) & all_markers$p_val_adj < padj_threshold, ,
    drop = FALSE
  ]

  for (grp in unique(as.character(sig$cluster))) {
    grp_df <- sig[as.character(sig$cluster) == grp, , drop = FALSE]
    grp_df <- grp_df[order(grp_df$avg_log2FC, decreasing = TRUE), ]
    safe_run(
      write_table(grp_df,
        file.path(deg_dir,
                  paste0("DEG_", file_prefix, "_",
                         .safe_filename(grp), "_by_", group_by_col, ".txt")),
        sep       = table_sep,
        quote     = table_quote,
        row_names = table_row_names),
      label = paste0("write per-group DEG: ", grp)
    )
  }

  freq_tbl             <- table(as.character(sig$cluster))
  deg_counts           <- empty_counts
  idx                  <- match(deg_counts$Group, names(freq_tbl))
  deg_counts$DEG_Count <- ifelse(is.na(idx), 0L, as.integer(freq_tbl[idx]))

  .save(deg_counts, "DEG_Counts")
  deg_counts
}

#' Compute and save average expression per group x condition
#'
#' Writes: AvgExpr_{layer}_{prefix}_by_{col}.txt
#'
#' @param seurat_obj Seurat object
#' @param output_dir Character
#' @param file_prefix Character
#' @param group_by_col Character
#' @param layer Character. "data" or "counts"
#' @param table_sep Character
#' @param table_quote Logical
#' @param table_row_names Logical
#' @return NULL
#' @export
save_average_expression <- function(seurat_obj, output_dir, file_prefix,
                                     group_by_col    = "seurat_clusters",
                                     layer           = "data",
                                     table_sep       = "\t",
                                     table_quote     = FALSE,
                                     table_row_names = FALSE) {
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))
  if (!"condition"  %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))

   Seurat::DefaultAssay(seurat_obj) <- "SCT"

  seurat_obj$group_condition <- paste(
    seurat_obj@meta.data[[group_by_col]],
    seurat_obj$condition,
    sep = "_"
  )
  Seurat::Idents(seurat_obj) <- "group_condition"

  avg <- get_avg_expr(seurat_obj, layer = layer)
  if (is.null(avg)) return(invisible(NULL))

  df      <- as.data.frame(avg)
  df$gene <- rownames(df)
  df      <- df[, c("gene", setdiff(colnames(df), "gene")), drop = FALSE]

  write_table(df,
    file.path(output_dir,
              paste0("AvgExpr_", layer, "_", file_prefix,
                     "_by_", group_by_col, ".txt")),
    sep       = table_sep,
    quote     = table_quote,
    row_names = table_row_names)
  invisible(NULL)
}

#' Compute and save pseudobulk counts per sample x group
#'
#' Writes: PseudobulkCounts_{prefix}_by_{col}.txt
#'
#' @param seurat_obj Seurat object
#' @param output_dir Character
#' @param file_prefix Character
#' @param group_by_col Character
#' @param table_sep Character
#' @param min_replicates Integer. Minimum unique samples required for pseudobulking
#' @return NULL
#' @export
save_pseudobulk_counts <- function(seurat_obj, output_dir, file_prefix,
                                    group_by_col = "seurat_clusters",
                                    table_sep    = "\t",
                                    min_replicates = 3L) {
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))
  
  sample_col <- if ("sample" %in% colnames(seurat_obj@meta.data)) "sample" else if ("orig.ident" %in% colnames(seurat_obj@meta.data)) "orig.ident" else NULL
  if (is.null(sample_col)) return(invisible(NULL))

  n_reps <- length(unique(seurat_obj@meta.data[[sample_col]]))
  if (n_reps < min_replicates) {
    message("   [SKIP] Pseudobulk: found ", n_reps, " replicates, requires at least ", min_replicates)
    return(invisible(NULL))
  }

  seurat_obj$pseudobulk_group <- paste(seurat_obj@meta.data[[group_by_col]],
                                       seurat_obj@meta.data[[sample_col]], sep = "_")
  Seurat::Idents(seurat_obj) <- "pseudobulk_group"

  agg <- safe_run({
    Seurat::AggregateExpression(seurat_obj, assays = "RNA", return.seurat = FALSE)[["RNA"]]
  }, label = "Pseudobulk AggregateExpression")

  if (is.null(agg)) return(invisible(NULL))

  df <- as.data.frame(agg)
  df$gene <- rownames(df)
  df <- df[, c("gene", setdiff(colnames(df), "gene")), drop = FALSE]

  write_table(df,
    file.path(output_dir,
              paste0("PseudobulkCounts_", file_prefix, "_by_", group_by_col, ".txt")),
    sep       = table_sep,
    quote     = FALSE,
    row_names = FALSE)
  invisible(NULL)
}

#' Compute cell type proportion table and chi squared test
#' @param seurat_obj Seurat object
#' @param grouping_col Character
#' @param output_dir Character
#' @param file_prefix Character
#' @return NULL
#' @export
run_proportion_analysis <- function(seurat_obj, grouping_col,
                                     output_dir, file_prefix) {
  if (!grouping_col %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))
  if (!"condition"  %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))

  make_dir(output_dir)

  tbl <- table(
    Group     = seurat_obj@meta.data[[grouping_col]],
    Condition = seurat_obj$condition
  )
  if (nrow(tbl) < 2L || ncol(tbl) < 2L) return(invisible(NULL))

  prop_df       <- as.data.frame(prop.table(tbl, margin = 2L))
  prop_df$Count <- as.integer(tbl)

  save_tsv(prop_df,
    file.path(output_dir,
              paste0("Proportions_", file_prefix, "_", grouping_col, ".txt")))

  chi <- safe_run(stats::chisq.test(tbl),
                  label = paste0("chisq.test_", grouping_col))
  if (!is.null(chi))
    writeLines(capture.output(print(chi)),
               file.path(output_dir,
                         paste0("ChiSq_", file_prefix, "_", grouping_col, ".txt")))

  p <- ggplot2::ggplot(prop_df,
         ggplot2::aes(x = Condition, y = Freq,
                      fill = Group)) +
       ggplot2::geom_col(position = "fill") +
       ggplot2::scale_y_continuous(labels = scales::percent_format()) +
       ggplot2::labs(y = "Proportion",
                     title = paste0(file_prefix, " | proportions by ", grouping_col)) +
       ggplot2::theme_minimal()

  save_png(p,
    file.path(output_dir,
              paste0("Proportion_", file_prefix, "_", grouping_col, ".png")))
  invisible(NULL)
}

#' scProportionTest
#' @param seurat_obj Seurat object
#' @param grouping_col Character
#' @param output_dir Character
#' @param file_prefix Character
#' @return NULL
#' @export
run_scproportion_test <- function(seurat_obj,
                                  grouping_col,
                                  output_dir,
                                  file_prefix) {

  md <- seurat_obj@meta.data

  if (!grouping_col %in% colnames(md)) return(invisible(NULL))
  if (!"sample" %in% colnames(md)) return(invisible(NULL))
  if (!"condition" %in% colnames(md)) return(invisible(NULL))

  conditions <- as.character(unique(md$condition))
  if (length(conditions) != 2) {
    message("   [SKIP] requires exactly 2 conditions")
    return(invisible(NULL))
  }

prop_test_result <- tryCatch({
    obj <- scProportionTest::sc_utils(seurat_obj)
    scProportionTest::permutation_test(
      obj,
      cluster_identity = grouping_col,
      sample_1         = conditions[1],  
      sample_2         = conditions[2],  
      sample_identity  = "condition"    
    )
  }, error = function(e) {
    message("   [ERROR] scProportionTest failed: ", e$message)
    message("   DEBUG: Object class is ", class(seurat_obj))
    message("   DEBUG: Conditions are: ", paste(conditions, collapse=" vs "))
    return(NULL)
  })

  if (!is.null(prop_test_result)) {
    out_file <- file.path(
      output_dir,
      paste0(
        "scProportionTest_",
        file_prefix,
        "_",
        grouping_col,
        ".tsv"
      )
    )

    if (is.data.frame(prop_test_result)) {
      write.table(
        prop_test_result,
        file = out_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    } else {
      saveRDS(prop_test_result, file = sub("\\.tsv$", ".rds", out_file))
    }
  }
  invisible(NULL)
}

#' Compute gene-gene correlation matrices per group and condition, plus global
#'
#' This function computes pairwise correlations between two gene sets
#' (genes_x vs genes_y) within each group and condition, and optionally
#' across all cells (global). Generates heatmaps and exports long-format tables.
#'
#' @param seurat_obj A Seurat object containing expression data.
#' @param grouping_col Character. Metadata column used to define groups.
#' @param genes_x Character vector. First gene set (rows of correlation matrix).
#' @param genes_y Character vector. Second gene set (columns of correlation matrix).
#' @param out_dir Character. Output directory where results will be saved.
#' @param prefix Character. Prefix added to output file names.
#' @param cond_col Character. Metadata column used to define conditions
#'   (default: "condition").
#' @param method Character. Correlation method to use ("pearson" or "spearman").
#' @param global_plot Logical. If TRUE, generate an overall correlation heatmap
#'   using all cells (ignoring grouping_col and cond_col). Default TRUE.
#'
#' @return NULL (invisibly). Writes heatmaps and correlation tables to disk.
#'
#' @export
#'
run_gene_correlations <- function(seurat_obj, 
                                  grouping_col, 
                                  genes_x, 
                                  genes_y,
                                  out_dir, 
                                  prefix,
                                  cond_col = "condition",
                                  method = "pearson",
                                  global_plot = TRUE) {

  message("--- Starting gene correlation analysis (Grouping: ", grouping_col, ") ---")
  
  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  corr_dir <- file.path(out_dir, "Correlation")
  if (!dir.exists(corr_dir)) dir.create(corr_dir, recursive = TRUE)

  groups <- unique(as.character(seurat_obj@meta.data[[grouping_col]]))
  groups <- groups[!is.na(groups) & !grepl("Unassigned", groups)]
  
  conditions <- unique(as.character(seurat_obj@meta.data[[cond_col]]))
  conditions <- conditions[!is.na(conditions)]
  expr <- Seurat::GetAssayData(seurat_obj, assay = "SCT", layer = "data")

  # --- Report original and missing genes for set X ---
  message("   Gene set X (", length(genes_x), " genes provided):")
  present_x <- intersect(genes_x, rownames(expr))
  missing_x <- setdiff(genes_x, rownames(expr))
  if (length(present_x) > 0) {
    message("      -> Found: ", paste(present_x, collapse = ", "))
  } else {
    message("      -> Found: none")
  }
  if (length(missing_x) > 0) {
    message("      -> Missing (not in SCT assay): ", paste(missing_x, collapse = ", "))
  }

  # --- Report original and missing genes for set Y ---
  message("   Gene set Y (", length(genes_y), " genes provided):")
  present_y <- intersect(genes_y, rownames(expr))
  missing_y <- setdiff(genes_y, rownames(expr))
  if (length(present_y) > 0) {
    message("      -> Found: ", paste(present_y, collapse = ", "))
  } else {
    message("      -> Found: none")
  }
  if (length(missing_y) > 0) {
    message("      -> Missing (not in SCT assay): ", paste(missing_y, collapse = ", "))
  }

  # --- Use only present genes ---
  valid_genes_x <- present_x
  valid_genes_y <- present_y

  if (length(valid_genes_x) < 1 || length(valid_genes_y) < 1) {
    message(" -> Skipping Correlation: At least one gene set is empty after filtering.")
    return(invisible(NULL))
  }

  # ========== GLOBAL CORRELATION (all cells, ignoring groups) ==========
  if (global_plot) {
    message("   Computing global correlation across all cells...")
    all_cells <- colnames(expr)
    mat_x_all <- t(as.matrix(expr[valid_genes_x, all_cells, drop = FALSE]))
    mat_y_all <- t(as.matrix(expr[valid_genes_y, all_cells, drop = FALSE]))
    cor_global <- suppressWarnings(cor(mat_x_all, mat_y_all, method = method))
    cor_global[is.na(cor_global)] <- 0
    cor_global_df <- as.data.frame(as.table(cor_global))
    colnames(cor_global_df) <- c("Gene_X", "Gene_Y", "Correlation")
    cor_global_df$Method <- method
    cor_global_df$Type <- "Global"
    
    # Save global correlation table
    global_tsv <- file.path(corr_dir, paste0("Correlations_", method, "_", prefix, "_global.tsv"))
    write.table(cor_global_df, global_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    message("   -> Global correlation table saved: ", basename(global_tsv))
    
    # Generate global heatmap
    tryCatch({
      p_global <- ggplot2::ggplot(cor_global_df, ggplot2::aes(x = Gene_Y, y = Gene_X, fill = Correlation)) +
        ggplot2::geom_tile(color = "darkgray", linewidth = 0.6) + 
        ggplot2::geom_text(
          ggplot2::aes(label = sprintf("%.2f", Correlation)), 
          size = 3.8, 
          color = "black"
        ) +
        ggplot2::scale_fill_gradient2(
          low = "blue", mid = "white", high = "red", midpoint = 0,
          limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
          labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
          guide = ggplot2::guide_colorbar(
            title = "Corr", title.position = "top", title.hjust = 0.5,
            barwidth = ggplot2::unit(0.45, "cm"), barheight = ggplot2::unit(4.2, "cm"),
            ticks.colour = "white", ticks.linewidth = 0.8
          )
        ) +
        ggplot2::coord_fixed() +
        ggplot2::scale_y_discrete(limits = rev(unique(cor_global_df$Gene_X))) +
        ggplot2::scale_x_discrete(limits = unique(cor_global_df$Gene_Y)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = paste0(toupper(method), " Correlation (Global)"),
          subtitle = paste0("all cells (n=", length(all_cells), ")"),
          x = NULL, y = NULL
        ) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5, margin = ggplot2::margin(b = 4)),
          plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, margin = ggplot2::margin(b = 12)),
          axis.text.y = ggplot2::element_text(size = 11, color = "black"),
          axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11, color = "black"),
          panel.grid = ggplot2::element_blank(),
          legend.title = ggplot2::element_text(size = 11, face = "bold"),
          legend.text = ggplot2::element_text(size = 10),
          legend.margin = ggplot2::margin(l = 12)
        )
      
      global_png <- file.path(corr_dir, paste0("Corr_", method, "_", prefix, "_global.png"))
      ggplot2::ggsave(filename = global_png, plot = p_global, width = 5.6, height = 4.6, dpi = 300, bg = "white")
      message("   -> Global heatmap saved: ", basename(global_png))
    }, error = function(e) {
      message("   [WARNING] Global heatmap could not be generated: ", e$message)
    })
  }

  # ========== PER-GROUP AND PER-CONDITION CORRELATIONS ==========
  all_cors <- list()

  for (grp in groups) {
    for (cond in conditions) {
      
      cells <- rownames(seurat_obj@meta.data[
        seurat_obj@meta.data[[grouping_col]] == grp &
        seurat_obj@meta.data[[cond_col]] == cond, 
      ])

      n_cells <- length(cells)
      grp_cond_label <- paste0(grp, "_", cond)

      if (n_cells >= 30) {
        message(paste("Processing Cluster:", grp, "| Cond:", cond, "| n:", n_cells))
        mat_x <- t(as.matrix(expr[valid_genes_x, cells, drop = FALSE]))
        mat_y <- t(as.matrix(expr[valid_genes_y, cells, drop = FALSE]))
        cor_res <- suppressWarnings(cor(mat_x, mat_y, method = method))
        cor_res[is.na(cor_res)] <- 0
        cor_df <- as.data.frame(as.table(cor_res))
        colnames(cor_df) <- c("Gene_X", "Gene_Y", "Correlation")
        cor_df$Cluster <- grp
        cor_df$Condition <- cond
        cor_df$N_Cells <- n_cells
        cor_df$Method <- method
        all_cors[[grp_cond_label]] <- cor_df
        
        # Generate heatmap for this group/condition
        tryCatch({
          clean_label <- gsub("[^A-Za-z0-9]", "_", grp_cond_label)
          out_file <- file.path(corr_dir, paste0("Corr_", method, "_", prefix, "_", clean_label, ".png"))

          p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Gene_Y, y = Gene_X, fill = Correlation)) +
            ggplot2::geom_tile(color = "darkgray", linewidth = 0.6) + 
            ggplot2::geom_text(
              ggplot2::aes(label = sprintf("%.2f", Correlation)), 
              size = 3.8, 
              color = "black"
            ) +
            ggplot2::scale_fill_gradient2(
              low = "blue", mid = "white", high = "red", midpoint = 0,
              limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
              labels = c("-1.0", "-0.5", "0.0", "0.5", "1.0"),
              guide = ggplot2::guide_colorbar(
                title = "Corr", title.position = "top", title.hjust = 0.5,
                barwidth = ggplot2::unit(0.45, "cm"), barheight = ggplot2::unit(4.2, "cm"),
                ticks.colour = "white", ticks.linewidth = 0.8
              )
            ) +
            ggplot2::coord_fixed() +
            ggplot2::scale_y_discrete(limits = rev(unique(cor_df$Gene_X))) +
            ggplot2::scale_x_discrete(limits = unique(cor_df$Gene_Y)) +
            ggplot2::theme_minimal() +
            ggplot2::labs(
              title = paste0(toupper(method), " Corr: ", grp),
              subtitle = paste0("(", cond, ", n=", n_cells, ")"),
              x = NULL, y = NULL
            ) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold", size = 15, hjust = 0.5, margin = ggplot2::margin(b = 4)),
              plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, margin = ggplot2::margin(b = 12)),
              axis.text.y = ggplot2::element_text(size = 11, color = "black"),
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 11, color = "black"),
              panel.grid = ggplot2::element_blank(),
              legend.title = ggplot2::element_text(size = 11, face = "bold"),
              legend.text = ggplot2::element_text(size = 10),
              legend.margin = ggplot2::margin(l = 12)
            )
          ggplot2::ggsave(filename = out_file, plot = p, width = 5.6, height = 4.6, dpi = 300, bg = "white")
        }, error = function(e) {
          message(" -> Heatmap skipped for ", grp_cond_label, ": ", e$message)
        })     
      } else {
        message(paste("Skipping:", grp_cond_label, "| Insufficient cells (N =", n_cells, ")"))
      }
    }
  }

  # Combine all per-group correlation tables
  if (length(all_cors) > 0) {
    final_df <- do.call(rbind, all_cors)
    tsv_file <- file.path(corr_dir, paste0("Correlations_", method, "_", prefix, "_by_group.tsv"))
    write.table(final_df, tsv_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message("--- Per-group correlation table saved to: ", basename(tsv_file), " ---")
  } else {
    message("--- No groups met the cell count threshold (N >= 30). ---")
  }

  message("--- Gene correlation analysis finished ---")
  invisible(NULL)
}