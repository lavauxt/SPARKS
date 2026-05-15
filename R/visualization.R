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

  ggplot2::ggsave(
    filename = file.path(
      out_dir,
      paste0(prefix, "_FeaturePlot_", gene, ".png")
    ),
    plot = p,
    width = dynamic_width,
    height = 5
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

  ht <- safe_run(
    Seurat::DoHeatmap(seurat_obj, features = top_genes) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
      ggplot2::ggtitle(paste0("Top", top_n, " markers  ",
                               prefix, " | ", group_by_col)),
    label = "DoHeatmap"
  )
  if (!is.null(ht))
    save_png(ht,
      file.path(out_dir,
                paste0("Heatmap_", prefix, "_", group_by_col, ".png")),
      width = 14, height = 10)
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