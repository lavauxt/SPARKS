#' Save a FeaturePlot for each gene in the list
#' @param seurat_obj Seurat object
#' @param genes Character vector
#' @param out_dir Character
#' @param prefix Character
#' @param reduction Character
#' @return NULL
#' @export
generate_feature_plots <- function(seurat_obj, genes, out_dir, prefix,
                                   reduction = "umap") {

  genes <- .filter_present_genes(genes, seurat_obj, "FeaturePlot")
  if (is.null(genes)) return(invisible(NULL))

  if (!.has_reduction(seurat_obj, reduction)) {
    message("   [SKIP FeaturePlot] reduction '", reduction, "' not found")
    return(invisible(NULL))
  }

  make_dir(out_dir)

for (gene in genes) {
  p <- Seurat::FeaturePlot(
    seurat_obj,
    features = gene,
    reduction = reduction,
    split.by = "condition",
    pt.size = 1,
    keep.scale = "all"
  ) +
    patchwork::plot_annotation(
      title = paste("Expression of", gene, "in", prefix)
    )

  dynamic_width <- max(
    6,
    length(unique(seurat_obj$condition[!is.na(seurat_obj$condition)])) * 5 + 1
  )

  save_png(p,
    filename = file.path(out_dir, paste0(prefix, "_FeaturePlot_", gene, ".png")),
    width    = dynamic_width,
    height   = 5
  )
}

  invisible(NULL)
}


#' Save VlnPlots split by condition for each gene
#' @param seurat_obj Seurat object
#' @param genes Character vector
#' @param out_dir Character
#' @param prefix Character
#' @param group_by_col Character. Metadata column for x-axis groups
#' @param split_by Character or NULL. Metadata column to split violins
#' @return NULL
#' @export
generate_violin_plots <- function(seurat_obj, genes, out_dir, prefix,
                                   group_by_col = "seurat_clusters",
                                   split_by     = "condition") {

  genes <- .filter_present_genes(genes, seurat_obj, "VlnPlot")
  if (is.null(genes)) return(invisible(NULL))

  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP VlnPlot] '", group_by_col, "' not in metadata")
    return(invisible(NULL))
  }
  use_split <- !is.null(split_by) &&
               split_by %in% colnames(seurat_obj@meta.data) &&
               length(unique(seurat_obj@meta.data[[split_by]])) > 1L

  make_dir(out_dir)
  for (g in genes) {
    p <- safe_run({
      args <- list(
        seurat_obj,
        features = g,
        group.by = group_by_col,
        pt.size = 0
      )
      if (use_split) {
        args$split.by  <- split_by
        args$split.plot <- TRUE
      }
      do.call(Seurat::VlnPlot, args)
    }, label = paste0("VlnPlot   ", g))
    if (!is.null(p)) {
      save_png(
        p,
        file.path(out_dir,
                  paste0("VlnPlot_", prefix, "_", g, ".png"))
      )
    }
  }
  invisible(NULL)
}

#' UMAP coloured by number of DEGs per cluster
#' @param seurat_obj Seurat object
#' @param deg_counts data.frame with columns Group and DEG_Count, or NULL
#' @param umap_dir Character
#' @param file_prefix Character
#' @param group_by_col Character
#' @param color_high Character
#' @param color_low Character
#' @param min_deg_display Integer. Min DEG_Count to label a cluster
#' @param reduction Character
#' @param point_size Numeric
#' @param alpha Numeric
#' @param label_size Numeric
#' @param width Numeric
#' @param height Numeric
#' @return NULL
#' @export
generate_deg_umap <- function(seurat_obj, deg_counts, umap_dir, file_prefix,
                               group_by_col,
                               color_high      = "red",
                               color_low       = "lightgrey",
                               min_deg_display = 5L,
                               reduction       = "umap",
                               point_size      = 0.5,
                               alpha           = 0.6,
                               label_size      = 4,
                               width           = 8,
                               height          = 7) {
  if (is.null(deg_counts) || nrow(deg_counts) == 0L) return(invisible(NULL))
  if (!any(deg_counts$DEG_Count >= min_deg_display)) {
    message("   [SKIP DEG UMAP] no groups with >= ", min_deg_display, " DEGs")
    return(invisible(NULL))
  }
  if (!.has_reduction(seurat_obj, reduction)) {
    message("   [SKIP DEG UMAP] reduction '", reduction, "' not found")
    return(invisible(NULL))
  }
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))

  emb <- as.data.frame(Seurat::Embeddings(seurat_obj, reduction = reduction))
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  emb[[group_by_col]] <- as.character(seurat_obj@meta.data[[group_by_col]])

  deg_counts$Group    <- as.character(deg_counts$Group)
  emb                 <- merge(emb, deg_counts, by.x = group_by_col,
                                by.y = "Group", all.x = TRUE)
  emb$DEG_Count[is.na(emb$DEG_Count)] <- 0L

  centroids <- stats::aggregate(cbind(UMAP_1, UMAP_2) ~ emb[[group_by_col]],
                                 data = emb, FUN = mean)
  colnames(centroids)[1L] <- group_by_col
  centroids <- merge(centroids, deg_counts, by.x = group_by_col,
                     by.y = "Group", all.x = TRUE)
  centroids$DEG_Count[is.na(centroids$DEG_Count)] <- 0L
  centroids_lab <- centroids[centroids$DEG_Count >= min_deg_display, ]

  p <- ggplot2::ggplot(emb,
         ggplot2::aes(x = UMAP_1, y = UMAP_2, colour = DEG_Count)) +
       ggplot2::geom_point(size = point_size, alpha = alpha) +
       ggplot2::scale_colour_gradient(low = color_low, high = color_high) +
       ggplot2::theme_minimal() +
       ggplot2::ggtitle(paste0(file_prefix, " | DEG count per ", group_by_col))

  if (nrow(centroids_lab) > 0L && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_label_repel(
      data    = centroids_lab,
      mapping = ggplot2::aes(label = paste0(.data[[group_by_col]],
                                             " (", DEG_Count, ")")),
      size    = label_size,
      colour  = "black",
      fill    = ggplot2::alpha("white", 0.7)
    )
  }

  make_dir(umap_dir)
  save_png(p, file.path(umap_dir,
    paste0("UMAP_DEG_", file_prefix, "_", group_by_col, ".png")),
    width = width, height = height)
  invisible(NULL)
}

#' FindAllMarkers heatmap for one grouping
#' @param seurat_obj Seurat object
#' @param group_by_col Character
#' @param out_dir Character
#' @param prefix Character
#' @param top_n Integer. Top markers per cluster
#' @param min_pct Numeric
#' @param logfc_threshold Numeric
#' @return NULL
#' @export
generate_cluster_markers_and_heatmap <- function(seurat_obj, group_by_col,
                                                  out_dir, prefix,
                                                  top_n           = 5L,
                                                  min_pct         = 0.25,
                                                  logfc_threshold = 0.25) {
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) return(invisible(NULL))

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  Seurat::Idents(seurat_obj) <- group_by_col

  markers <- .find_all_markers_safe(seurat_obj, only_pos = TRUE,
                                     min_pct = min_pct,
                                     logfc_threshold = logfc_threshold)
  if (is.null(markers)) return(invisible(NULL))

  make_dir(out_dir)
  save_tsv(markers,
    file.path(out_dir, paste0("ClusterMarkers_", prefix, "_", group_by_col, ".txt")))

  top_genes <- markers |>
    (\(df) split(df, df$cluster))() |>
    lapply(function(x) head(x[order(x$avg_log2FC, decreasing = TRUE), "gene"], top_n)) |>
    unlist() |>
    unique()

  if (length(top_genes) == 0L) return(invisible(NULL))

  if (!.has_scale_data(seurat_obj)) {
    seurat_obj <- safe_run(Seurat::ScaleData(seurat_obj, verbose = FALSE),
                           label = "ScaleData for heatmap")
    if (is.null(seurat_obj)) return(invisible(NULL))
  }

  n_groups    <- length(unique(as.character(seurat_obj@meta.data[[group_by_col]])))
  label_style <- .heatmap_label_params(n_groups)

  ht <- safe_run(
    Seurat::DoHeatmap(seurat_obj, features = top_genes,
                      group.by = group_by_col,
                      angle    = label_style$angle,
                      hjust    = label_style$hjust,
                      size     = label_style$size) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
      ggplot2::ggtitle(paste0("Top", top_n, " markers  ",
                               prefix, " | ", group_by_col)),
    label = "DoHeatmap"
  )
  if (!is.null(ht))
    save_png(ht,
      file.path(out_dir,
                paste0("Heatmap_", prefix, "_", group_by_col, ".png")),
      width = min(50, max(14, n_groups * 0.4 + 6)), height = 10)
  invisible(NULL)
}

#' Average-expression heatmap (pheatmap) per group x condition
#'
#' Writes: ExprHeatmap_{prefix}_{col}.png
#'
#' @param seurat_obj Seurat object
#' @param group_by_col Character
#' @param out_dir Character
#' @param prefix Character
#' @param species_target Character
#' @param top_n Integer. Top expressed genes per column
#' @return NULL
#' @export
generate_expression_heatmap <- function(seurat_obj, group_by_col, out_dir,
                                         prefix, species_target = "Mouse",
                                         top_n = 10L) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("   [SKIP ExprHeatmap] pheatmap not installed")
    return(invisible(NULL))
  }
  groups <- get_valid_groups(seurat_obj@meta.data, group_by_col)
  if (length(groups) < 2L) {
    message("   [SKIP ExprHeatmap] fewer than 2 valid groups in '",
            group_by_col, "'")
    return(invisible(NULL))
  }

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  seurat_obj$group_condition <- paste(
    seurat_obj@meta.data[[group_by_col]],
    seurat_obj$condition,
    sep = "_"
  )
  Seurat::Idents(seurat_obj) <- "group_condition"

  avg <- get_avg_expr(seurat_obj, layer = "data")
  if (is.null(avg) || nrow(avg) == 0L) return(invisible(NULL))

  junk_pat  <- get_junk_pattern(species_target)
  avg_filt  <- avg[!grepl(junk_pat, rownames(avg)), , drop = FALSE]
  if (nrow(avg_filt) < 2L) return(invisible(NULL))

  top_genes <- unique(unlist(lapply(seq_len(ncol(avg_filt)), function(i) {
    n <- min(top_n, nrow(avg_filt))
    names(sort(avg_filt[, i], decreasing = TRUE)[seq_len(n)])
  })))
  if (length(top_genes) < 2L) return(invisible(NULL))

  ph <- safe_run(
    pheatmap::pheatmap(
      as.matrix(avg_filt[top_genes, , drop = FALSE]),
      scale         = "row",
      cluster_cols  = FALSE,
      cluster_rows  = TRUE,
      show_rownames = TRUE,
      fontsize_row  = 6,
      main          = paste0("Top", top_n, " Expressed  ",
                              prefix, " | ", group_by_col),
      silent        = TRUE
    ),
    label = "pheatmap expression heatmap"
  )
  if (is.null(ph)) return(invisible(NULL))

  make_dir(out_dir)
  save_png(ph$gtable,
    file.path(out_dir,
              paste0("ExprHeatmap_", prefix, "_", group_by_col, ".png")),
    width = 10, height = 12)
  invisible(NULL)
}

#' Top expressed genes table per group
#'
#' Writes: TopExpressed_{prefix}_{col}.txt
#'
#' @param seurat_obj Seurat object
#' @param group_by_col Character
#' @param out_dir Character
#' @param prefix Character
#' @param top_n Integer
#' @return NULL
#' @export
generate_top_expressed_genes <- function(seurat_obj, group_by_col,
                                          out_dir, prefix, top_n = 10L) {
  valid_cells <- !is.na(seurat_obj@meta.data[[group_by_col]]) &
                 seurat_obj@meta.data[[group_by_col]] != "Unassigned"
  if (sum(valid_cells) < 2L) return(invisible(NULL))

  so_sub <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
  Seurat::DefaultAssay(so_sub) <- "SCT"
  Seurat::Idents(so_sub) <- group_by_col

  avg <- get_avg_expr(so_sub, layer = "data")
  if (is.null(avg)) return(invisible(NULL))

  top_df <- do.call(rbind, lapply(colnames(avg), function(cluster) {
    n   <- min(top_n, nrow(avg))
    top <- sort(avg[, cluster], decreasing = TRUE)[seq_len(n)]
    data.frame(
      Group             = cluster,
      Gene              = names(top),
      AverageExpression = as.numeric(top),
      stringsAsFactors  = FALSE,
      row.names         = NULL
    )
  }))

  make_dir(out_dir)
  save_tsv(top_df,
    file.path(out_dir,
              paste0("TopExpressed_", prefix, "_", group_by_col, ".txt")))
  invisible(NULL)
}

#' Generate Subcluster Heatmaps (Per Cell and Aggregated)
#'
#' @param seurat_obj Seurat object
#' @param genes Vector of genes to plot
#' @param out_dir Output directory path
#' @param prefix Prefix for the file names
#' @export
generate_subcluster_heatmaps <- function(seurat_obj, genes, out_dir, prefix) {
  make_dir(out_dir)
  
  genes <- .filter_present_genes(genes, seurat_obj, "Subcluster Heatmaps")
  if (is.null(genes) || length(genes) < 2) {
    message("   [SKIP] Not enough valid genes for Subcluster Heatmap: ", prefix)
    return(invisible(NULL))
  }

  # ---------------------------------------------------------
  # 1. Per-Cell Heatmap (DoHeatmap)
  # ---------------------------------------------------------
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = genes, verbose = FALSE)
  n_cond_groups <- length(unique(as.character(seurat_obj$condition)))
  cond_label_style <- .heatmap_label_params(n_cond_groups)
  p_cell <- Seurat::DoHeatmap(
    seurat_obj, 
    features = genes, 
    group.by = "condition", 
    angle = cond_label_style$angle,
    hjust = cond_label_style$hjust,
    size = cond_label_style$size
  ) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0, name = "z-score"
    ) +
    ggplot2::guides(color = "none") + 
    ggplot2::ggtitle(paste0(prefix, " - z-score per cell")) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))

  ggplot2::ggsave(
    file.path(out_dir, paste0(prefix, "_Heatmap_PerCell.png")), 
    plot = p_cell, width = 10, height = 8
  )

  # ---------------------------------------------------------
  # 2. Aggregated Heatmap (ggplot2)
  # ---------------------------------------------------------
  mat <- Seurat::GetAssayData(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj), layer = "scale.data")
  mat <- mat[genes, , drop = FALSE]

  conditions <- seurat_obj$condition
  
  cond_levels <- levels(factor(conditions))
  agg_mat <- vapply(cond_levels, function(cond) {
    cells <- names(conditions)[conditions == cond]
    if (length(cells) == 0L) return(rep(NA_real_, nrow(mat)))
    if (length(cells) == 1L) return(as.numeric(mat[, cells]))
    rowMeans(mat[, cells, drop = FALSE], na.rm = TRUE)
  }, FUN.VALUE = numeric(nrow(mat)))
  rownames(agg_mat) <- rownames(mat)
  
  df_long <- as.data.frame(agg_mat)
  df_long$gene <- rownames(df_long)
  df_long <- tidyr::pivot_longer(df_long, cols = -gene, names_to = "Condition", values_to = "z_score")

  df_long$Condition <- factor(df_long$Condition, levels = levels(factor(conditions)))
  df_long$gene      <- factor(df_long$gene, levels = rev(genes))

  p_agg <- ggplot2::ggplot(df_long, ggplot2::aes(x = Condition, y = gene, fill = z_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      name = "z-score",
      breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
      labels = c("-2", "-1", "-0.5", "0", "0.5", "1", "2"),
      limits = c(-2, 2),
      oob = scales::squish,
      guide = ggplot2::guide_colorbar(
        title.position = "top", 
        title.hjust = 0.5
      )
    ) +
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y      = ggplot2::element_text(face = "italic", size = 10),
      panel.grid       = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position  = "right",
      legend.title     = ggplot2::element_text(size = 10, face = "bold")
    ) +
    ggplot2::labs(
      title = paste0(prefix, " - average z-score"),
      x = NULL,
      y = NULL
    )

  calc_height <- max(4, length(genes) * 0.4)

  ggplot2::ggsave(
    file.path(out_dir, paste0(prefix, "_Heatmap_Aggregated.png")), 
    plot = p_agg, width = 5, height = calc_height
  )
}

#' Generate Gene Signature Panel (per-cell heatmap + aggregated heatmap + DotPlot)
#'
#' Reproduces, for an arbitrary named gene list, the same three views used to
#' inspect a curated gene panel: a per-cell z-score heatmap (DoHeatmap), a
#' group-averaged z-score heatmap, and a Seurat DotPlot (average expression +
#' percent expressed). Designed to be called once per grouping column, so it
#' naturally covers every cluster/subcluster grouping when invoked from
#' \code{run_grouping_analysis()} (Main analysis unit AND every subset).
#'
#' The per-cell heatmap and the aggregated heatmap nest \code{condition_col}
#' inside each \code{group_by_col} level (e.g. "B cells | WT", "B cells | KO"),
#' so WT/KO stays visually separated within every cluster/label instead of
#' being averaged together. The DotPlot uses Seurat's own \code{split.by} for
#' the same purpose. This nesting is skipped automatically when
#' \code{group_by_col} and \code{condition_col} are the same column, or when
#' \code{condition_col} has fewer than 2 levels — nesting a value against
#' itself (e.g. an already-computed fold-change/ratio) would be meaningless.
#'
#' @param seurat_obj Seurat object (DefaultAssay should already be set, e.g. "SCT")
#' @param genes Character vector of gene symbols making up the signature
#' @param out_dir Character. Output directory (typically the group's Heatmap folder)
#' @param prefix Character. File-prefix, e.g. display_name (e.g. "Main", "Fibroblasts")
#' @param group_by_col Character. Metadata column to group by (e.g. "seurat_clusters",
#'   "singleR_labels_main", or a subset's type_col).
#' @param signature_name Character. Name of the gene signature (used in filenames/titles)
#' @param condition_col Character. Metadata column holding the WT/KO-style condition,
#'   nested inside each group. Default "condition".
#' @return NULL
#' @export
generate_gene_signature_plots <- function(seurat_obj, genes, out_dir, prefix,
                                          group_by_col, signature_name,
                                          condition_col = "condition") {

  genes <- .filter_present_genes(genes, seurat_obj,
                                  paste0("Gene Signature '", signature_name, "'"))
  if (is.null(genes) || length(genes) < 2L) {
    message("   [SKIP] Not enough valid genes for signature '", signature_name,
            "' (", prefix, " | ", group_by_col, ")")
    return(invisible(NULL))
  }

  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP Gene Signature] '", group_by_col, "' not in metadata")
    return(invisible(NULL))
  }

  valid_cells <- !is.na(seurat_obj@meta.data[[group_by_col]])
  if (sum(valid_cells) < 2L) return(invisible(NULL))
  if (!all(valid_cells)) {
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
  }

  n_groups <- length(unique(as.character(seurat_obj@meta.data[[group_by_col]])))
  if (n_groups < 2L) {
    message("   [SKIP Gene Signature] '", group_by_col,
            "' has fewer than 2 groups for '", signature_name, "'")
    return(invisible(NULL))
  }

  # ── Decide whether condition should be nested inside each group ────────────
  # Skipped when condition_col IS the grouping column (nesting a variable
  # against itself), when it's absent, or when it has < 2 levels — and always
  # skipped for values that already encode a WT-vs-KO comparison internally
  # (e.g. a fold-change/ratio metric), since re-splitting those by condition
  # would not make sense.
  has_condition <- !identical(condition_col, group_by_col) &&
    condition_col %in% colnames(seurat_obj@meta.data) &&
    length(unique(as.character(stats::na.omit(seurat_obj@meta.data[[condition_col]])))) > 1L

  group_vec <- as.character(seurat_obj@meta.data[[group_by_col]])
  if (has_condition) {
    cond_vec   <- as.character(seurat_obj@meta.data[[condition_col]])
    combo      <- paste(group_vec, cond_vec, sep = " | ")
    combo_lvls <- as.vector(outer(sort(unique(group_vec)), sort(unique(cond_vec)),
                                  paste, sep = " | "))
    combo_lvls <- combo_lvls[combo_lvls %in% unique(combo)]
    seurat_obj$.sig_group <- factor(combo, levels = combo_lvls)
  } else {
    seurat_obj$.sig_group <- factor(group_vec, levels = sort(unique(group_vec)))
  }
  n_combo <- nlevels(seurat_obj$.sig_group)

  make_dir(out_dir)
  file_tag <- paste0(prefix, "_", group_by_col, "_", .safe_filename(signature_name))
  label_style <- .heatmap_label_params(n_combo)

  # ── 1. Per-cell z-score heatmap (DoHeatmap) ────────────────────────────────
  seurat_obj <- safe_run(
    Seurat::ScaleData(seurat_obj, features = genes, verbose = FALSE),
    label = paste0("ScaleData (", signature_name, ")"), fallback = seurat_obj
  )

  p_cell <- safe_run({
    Seurat::DoHeatmap(
      seurat_obj,
      features = genes,
      group.by = ".sig_group",
      angle    = label_style$angle,
      hjust    = label_style$hjust,
      size     = label_style$size
    ) +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, name = "z-score"
      ) +
      ggplot2::guides(color = "none") +
      ggplot2::ggtitle(paste0(prefix, " | ", signature_name, " - z-score per cell (",
                              group_by_col,
                              if (has_condition) paste0(" x ", condition_col) else "",
                              ")")) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 10))
  }, label = paste0("DoHeatmap (", signature_name, ")"))

  if (!is.null(p_cell)) {
    calc_width <- min(50, max(10, n_combo * 0.35 + 4))
    save_png(p_cell,
      file.path(out_dir, paste0(file_tag, "_Heatmap_PerCell.png")),
      width = calc_width, height = max(6, length(genes) * 0.5 + 2))
  }

  # ── 2. Aggregated (group-averaged) z-score heatmap ─────────────────────────
  agg_plot <- safe_run({
    mat <- Seurat::GetAssayData(seurat_obj, assay = Seurat::DefaultAssay(seurat_obj),
                                layer = "scale.data")
    mat <- mat[genes, , drop = FALSE]

    combo_vec  <- as.character(seurat_obj$.sig_group)
    combo_lvls <- levels(seurat_obj$.sig_group)
    agg_mat <- vapply(combo_lvls, function(grp) {
      cells <- colnames(seurat_obj)[combo_vec == grp]
      if (length(cells) == 0L) return(rep(NA_real_, nrow(mat)))
      if (length(cells) == 1L) return(as.numeric(mat[, cells]))
      rowMeans(mat[, cells, drop = FALSE], na.rm = TRUE)
    }, FUN.VALUE = numeric(nrow(mat)))
    rownames(agg_mat) <- rownames(mat)

    df_long <- as.data.frame(agg_mat)
    df_long$gene <- rownames(df_long)
    df_long <- tidyr::pivot_longer(df_long, cols = -gene,
                                   names_to = "Group", values_to = "z_score")
    df_long$Group <- factor(df_long$Group, levels = combo_lvls)
    df_long$gene  <- factor(df_long$gene,  levels = rev(genes))

    ggplot2::ggplot(df_long, ggplot2::aes(x = Group, y = gene, fill = z_score)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        name = "z-score",
        breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2),
        labels = c("-2", "-1", "-0.5", "0", "0.5", "1", "2"),
        limits = c(-2, 2),
        oob = scales::squish,
        guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5)
      ) +
      ggplot2::scale_x_discrete(position = "top") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(face = "bold", size = 9,
                                                 angle = label_style$angle,
                                                 hjust = label_style$hjust,
                                                 vjust = 0.5),
        axis.text.y     = ggplot2::element_text(face = "italic", size = 10),
        panel.grid      = ggplot2::element_blank(),
        plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right",
        legend.title    = ggplot2::element_text(size = 10, face = "bold")
      ) +
      ggplot2::labs(
        title = paste0(prefix, " | ", signature_name, " - average z-score (",
                       group_by_col, if (has_condition) paste0(" x ", condition_col) else "",
                       ")"),
        x = NULL, y = NULL
      )
  }, label = paste0("Aggregated heatmap (", signature_name, ")"))

  if (!is.null(agg_plot)) {
    calc_width <- min(50, max(5, n_combo * 0.55 + 2))
    save_png(agg_plot,
      file.path(out_dir, paste0(file_tag, "_Heatmap_Aggregated.png")),
      width = calc_width, height = max(4, length(genes) * 0.4))
  }

  # ── 3. DotPlot (average expression + percent expressed) ────────────────────
  # Updated: red-blue gradient, explicit legend title "Avg Expression"
  p_dot <- safe_run({
    dp <- if (has_condition) {
      # condition split -> use red-blue gradient for expression
      Seurat::DotPlot(
        seurat_obj,
        features = genes,
        group.by = group_by_col,
        split.by = condition_col,
        cols     = c("blue", "red")   # low = blue, high = red
      )
    } else {
      Seurat::DotPlot(seurat_obj,
                      features = genes,
                      group.by = group_by_col,
                      cols     = c("blue", "red"))
    }

    dp +
      Seurat::RotatedAxis() +
      ggplot2::labs(
        title = paste0(prefix, " | ", signature_name, " (",
                       group_by_col, if (has_condition) paste0(" x ", condition_col) else "",
                       ")"),
        x = "Features", y = group_by_col
      ) +
      ggplot2::guides(colour = ggplot2::guide_colorbar(title = "Avg Expression")) +  # explicit legend
      ggplot2::theme(
        plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.y = ggplot2::element_text(size = label_style$size * 2.2)
      )
  }, label = paste0("DotPlot (", signature_name, ")"))

  if (!is.null(p_dot)) {
    # wider to accommodate legend
    save_png(p_dot,
      file.path(out_dir, paste0(file_tag, "_DotPlot.png")),
      width  = max(10, length(genes) * 0.7 + 4),
      height = max(3, n_combo * 0.35 + 2))
  }

  invisible(NULL)
}

#' Generate Z-scored Heatmap of Top Expressed Genes per Cluster
#'
#' @param seurat_obj Seurat object
#' @param group_by_col Character. Metadata column for clusters (e.g. "seurat_clusters")
#' @param out_dir Character. Output directory
#' @param prefix Character. Prefix for file names
#' @param species_target Character. "Mouse" or "Human" to filter junk genes
#' @param top_n Integer. Number of top genes to select per cluster
#' @return NULL
#' @export
generate_cluster_zscore_heatmap <- function(seurat_obj, group_by_col, out_dir, prefix,
                                            species_target = "Mouse", top_n = 10L) {
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("   [SKIP Z-Score Heatmap] pheatmap not installed")
    return(invisible(NULL))
  }
  
  groups <- get_valid_groups(seurat_obj@meta.data, group_by_col)
  if (length(groups) < 2L) {
    message("   [SKIP Z-Score Heatmap] fewer than 2 valid groups in '", group_by_col, "'")
    return(invisible(NULL))
  }

  Seurat::Idents(seurat_obj) <- group_by_col
  avg <- get_avg_expr(seurat_obj, layer = "data")
  if (is.null(avg) || nrow(avg) == 0L) return(invisible(NULL))
  junk_pat <- get_junk_pattern(species_target)
  avg_filt <- avg[!grepl(junk_pat, rownames(avg)), , drop = FALSE]
  if (nrow(avg_filt) < 2L) return(invisible(NULL))
  top_genes <- unique(unlist(lapply(seq_len(ncol(avg_filt)), function(i) {
    n <- min(top_n, nrow(avg_filt))
    names(sort(avg_filt[, i], decreasing = TRUE)[seq_len(n)])
  })))
  
  if (length(top_genes) < 2L) return(invisible(NULL))

  breaks_list <- seq(-2, 2, by = 0.04)
  color_pal   <- grDevices::colorRampPalette(c("blue", "white", "red"))(length(breaks_list))

  ph <- safe_run(
    pheatmap::pheatmap(
      as.matrix(avg_filt[top_genes, , drop = FALSE]),
      scale         = "row",            
      cluster_cols  = TRUE,            
      cluster_rows  = TRUE,
      show_rownames = TRUE,
      fontsize_row  = max(5, 12 - length(top_genes)/15),
      main          = paste0("Top ", top_n, " Avg Genes (Z-Scored) | ", prefix, " | ", group_by_col),
      color         = color_pal,
      breaks        = breaks_list,
      silent        = TRUE
    ),
    label = "pheatmap z-score heatmap"
  )
  
  if (is.null(ph)) return(invisible(NULL))

  make_dir(out_dir)
  
  calc_height <- max(6, length(top_genes) * 0.15)
  calc_width  <- max(6, length(groups) * 0.5 + 2)

  save_png(ph$gtable,
    file.path(out_dir, paste0("ZScore_TopGenes_", prefix, "_", group_by_col, ".png")),
    width = calc_width, height = calc_height)
    
  invisible(NULL)
}

# =============================================================================
# NEW FUNCTION: Z-scored heatmap split by condition
# =============================================================================

#' Generate Z-scored Heatmap Split by Condition
#'
#' For each cluster (or label), compute average expression separately for each
#' condition and produce a heatmap where columns are cluster_condition
#' combinations and rows are genes, with row z-scores.
#'
#' @param seurat_obj Seurat object
#' @param group_by_col Character. Metadata column for clusters (e.g. "seurat_clusters")
#' @param condition_col Character. Metadata column for condition (default "condition")
#' @param out_dir Character. Output directory
#' @param prefix Character. Prefix for file names
#' @param species_target Character. "Mouse" or "Human" to filter junk genes
#' @param top_n Integer. Number of top genes to select per cluster_condition
#' @return NULL
#' @export
generate_cluster_zscore_heatmap_split_condition <- function(seurat_obj,
                                                            group_by_col,
                                                            condition_col = "condition",
                                                            out_dir,
                                                            prefix,
                                                            species_target = "Mouse",
                                                            top_n = 10L) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("   [SKIP Z-Score Heatmap (split by condition)] pheatmap not installed")
    return(invisible(NULL))
  }

  # Check columns
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP] '", group_by_col, "' not in metadata")
    return(invisible(NULL))
  }
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP] '", condition_col, "' not in metadata")
    return(invisible(NULL))
  }

  # Create combined group_condition identity
  seurat_obj$group_condition <- paste(
    seurat_obj@meta.data[[group_by_col]],
    seurat_obj@meta.data[[condition_col]],
    sep = " | "
  )

  # Only keep combinations that have at least 3 cells (optional)
  valid_groups <- get_valid_groups(seurat_obj@meta.data, "group_condition", min_cells = 3L)
  if (length(valid_groups) < 2L) {
    message("   [SKIP Z-Score Heatmap (split)] fewer than 2 valid cluster×condition groups")
    return(invisible(NULL))
  }

  Seurat::Idents(seurat_obj) <- "group_condition"
  avg <- get_avg_expr(seurat_obj, layer = "data")
  if (is.null(avg) || nrow(avg) == 0L) return(invisible(NULL))

  # Filter junk genes
  junk_pat <- get_junk_pattern(species_target)
  avg_filt <- avg[!grepl(junk_pat, rownames(avg)), , drop = FALSE]
  if (nrow(avg_filt) < 2L) return(invisible(NULL))

  # Select top genes per combined group
  top_genes <- unique(unlist(lapply(seq_len(ncol(avg_filt)), function(i) {
    n <- min(top_n, nrow(avg_filt))
    names(sort(avg_filt[, i], decreasing = TRUE)[seq_len(n)])
  })))
  if (length(top_genes) < 2L) return(invisible(NULL))

  # Build heatmap matrix (rows = genes, columns = group_condition)
  mat <- as.matrix(avg_filt[top_genes, , drop = FALSE])

  # Row z-score
  mat_scaled <- t(scale(t(mat)))

  # Define colour breaks
  breaks_list <- seq(-2, 2, by = 0.04)
  color_pal   <- grDevices::colorRampPalette(c("blue", "white", "red"))(length(breaks_list))

  ph <- safe_run(
    pheatmap::pheatmap(
      mat_scaled,
      scale         = "none",          # already scaled
      cluster_cols  = TRUE,
      cluster_rows  = TRUE,
      show_rownames = TRUE,
      fontsize_row  = max(5, 12 - length(top_genes) / 15),
      main          = paste0("Top ", top_n, " Genes (Z-Scored) | ", prefix, " | ",
                             group_by_col, " × ", condition_col),
      color         = color_pal,
      breaks        = breaks_list,
      silent        = TRUE
    ),
    label = "pheatmap z-score (split by condition)"
  )

  if (is.null(ph)) return(invisible(NULL))

  make_dir(out_dir)

  calc_height <- max(6, length(top_genes) * 0.15)
  calc_width  <- max(6, length(valid_groups) * 0.5 + 2)

  save_png(ph$gtable,
    file.path(out_dir,
              paste0("ZScore_TopGenes_", prefix, "_", group_by_col,
                     "_by_", condition_col, ".png")),
    width = calc_width, height = calc_height)

  invisible(NULL)
}

# =============================================================================
# NEW FUNCTION: Per-group DotPlots for a gene signature
# =============================================================================

#' Generate per-group DotPlots for a gene signature
#'
#' For each group in `group_by_col`, subset the cells and produce a Seurat
#' DotPlot showing expression of the signature genes, split by condition.
#' One PNG file is saved per group.
#'
#' @param seurat_obj Seurat object
#' @param genes Character vector. Signature genes.
#' @param out_dir Character. Output directory.
#' @param prefix Character. Prefix for file names.
#' @param group_by_col Character. Metadata column defining groups (clusters, labels).
#' @param signature_name Character. Name of the signature (used in titles/filenames).
#' @param condition_col Character. Metadata column for condition (e.g., "condition").
#' @param min_cells_per_group Integer. Minimum number of cells in a group to plot.
#' @param dot_colors Character vector of two colours for condition split (optional).
#' @return NULL, invisibly.
#' @export
generate_gene_signature_per_group_dotplots <- function(seurat_obj,
                                                       genes,
                                                       out_dir,
                                                       prefix,
                                                       group_by_col,
                                                       signature_name,
                                                       condition_col = "condition",
                                                       min_cells_per_group = 10L,
                                                       dot_colors = NULL) {

  genes <- .filter_present_genes(genes, seurat_obj,
                                 paste0("Per-group DotPlot for '", signature_name, "'"))
  if (is.null(genes) || length(genes) < 1L) {
    message("   [SKIP] No valid genes for signature '", signature_name,
            "' per-group dot plots.")
    return(invisible(NULL))
  }

  # Check columns
  if (!group_by_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP Per-group DotPlot] '", group_by_col, "' not in metadata")
    return(invisible(NULL))
  }
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP Per-group DotPlot] '", condition_col, "' not in metadata")
    return(invisible(NULL))
  }

  # Get groups with enough cells (global, but we'll subset later)
  all_groups <- as.character(unique(seurat_obj@meta.data[[group_by_col]]))
  all_groups <- all_groups[!is.na(all_groups) & all_groups != "Unassigned"]  # optional filter

  make_dir(out_dir)

  # Determine dot colours if not provided
  if (is.null(dot_colors)) {
    n_cond <- length(unique(seurat_obj@meta.data[[condition_col]]))
    dot_colors <- scales::hue_pal()(n_cond)
  }

  for (grp in all_groups) {
    # Subset cells for this group
    cells_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_by_col]] == grp]
    if (length(cells_keep) < min_cells_per_group) {
      message("   [SKIP] Group '", grp, "' has only ", length(cells_keep),
              " cells (min = ", min_cells_per_group, ")")
      next
    }

    # Create temporary object with only these cells
    sub_obj <- tryCatch(
      subset(seurat_obj, cells = cells_keep),
      error = function(e) NULL
    )
    if (is.null(sub_obj) || ncol(sub_obj) == 0L) next

    # Check that both conditions are present in this subset (optional, but we can still plot)
    conds_present <- unique(sub_obj@meta.data[[condition_col]])
    if (length(conds_present) < 2L) {
      message("   [SKIP] Group '", grp, "' has only one condition (",
              paste(conds_present, collapse = ", "), ") – skipping dot plot.")
      next
    }

    # Create DotPlot – red-blue gradient, explicit legend title
    p <- tryCatch({
      Seurat::DotPlot(sub_obj,
                      features = genes,
                      group.by = condition_col,    # conditions on y-axis
                      cols     = c("blue", "red")) +   # low=blue, high=red
        Seurat::RotatedAxis() +
        ggplot2::labs(
          title = paste0(prefix, " | ", signature_name, " in ", grp,
                         " (", paste(conds_present, collapse = " vs "), ")"),
          x = "Genes",
          y = condition_col
        ) +
        ggplot2::guides(colour = ggplot2::guide_colorbar(title = "Avg Expression")) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    }, error = function(e) NULL)

    if (!is.null(p)) {
      safe_grp <- .safe_filename(grp)
      filename <- file.path(out_dir,
                            paste0(prefix, "_", group_by_col, "_",
                                   signature_name, "_", safe_grp, "_DotPlot.png"))
      save_png(p, filename,
               width = max(8, length(genes) * 0.6 + 3),
               height = max(4, 1 + length(conds_present) * 0.5))
    }
  }

  invisible(NULL)
}