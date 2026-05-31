#' Save text summary of cell counts per sample and cluster
#'
#' Writes `{prefix}_EMPTY.txt` when ncol == 0 so downstream file-existence
#' checks always find *something*.
#'
#' @param object A Seurat object or NULL
#' @param prefix Character. Prefix for output file name
#' @param outdir Character. Output directory path
#' @return NULL invisibly
#' @export
save_cell_counts <- function(object, prefix, outdir) {
  make_dir(outdir)
  if (is.null(object) || ncol(object) == 0L) {
    writeLines("No cells remaining.",
               file.path(outdir, paste0(prefix, "_EMPTY.txt")))
    return(invisible(NULL))
  }

  path  <- file.path(outdir, paste0(prefix, "_cell_counts.txt"))
  lines <- paste("Total cells:", ncol(object))

  if ("orig.ident" %in% colnames(object@meta.data)) {
    id_tab <- table(object$orig.ident)
    lines  <- c(lines, "\nCells by orig.ident:",
                paste(names(id_tab), id_tab, sep = ": ", collapse = "\n"))

    if ("seurat_clusters" %in% colnames(object@meta.data)) {
      cluster_tab  <- table(object$seurat_clusters, object$orig.ident)
      active_prots <- colnames(cluster_tab)[colSums(cluster_tab) > 0L]
      lines        <- c(lines, "\nCells per cluster by orig.ident:")
      for (prot in active_prots) {
        counts <- cluster_tab[cluster_tab[, prot] > 0L, prot]
        if (length(counts) > 0L)
          lines <- c(lines,
                     paste0("\norig.ident: ", prot),
                     paste(names(counts), counts, sep = ": ", collapse = "\n"))
      }
    }
  }
  writeLines(lines, path)
  invisible(NULL)
}

#' Generate standard QC violin and scatter plots
#'
#' @param so A Seurat object
#' @param qc_dir Character. Output directory
#' @param prefix Character. Prefix for output files
#' @param cfg Named list. Pipeline config (for plot dimensions)
#' @return NULL invisibly
#' @export
generate_qc_plots <- function(so, qc_dir, prefix, cfg) {
  if (is.null(so) || ncol(so) == 0L) return(invisible(NULL))
  make_dir(qc_dir)

  vln <- safe_run(
    Seurat::VlnPlot(so,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3L, pt.size = 1),
    label = paste0("QC VlnPlot   ", prefix)
  )
  if (!is.null(vln))
    save_png(vln,
      file.path(qc_dir, paste0("VlnPlot_", prefix, ".png")),
      width  = cfg$plot$qc_vln_width  %||% 12,
      height = cfg$plot$qc_vln_height %||% 5)

  scatter <- safe_run({
    p1 <- Seurat::FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
    p2 <- Seurat::FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p1 + p2
  }, label = paste0("QC FeatureScatter  ", prefix))
  if (!is.null(scatter))
    save_png(scatter,
      file.path(qc_dir, paste0("FeatureScatter_", prefix, ".png")),
      width  = cfg$plot$qc_scatter_width  %||% 10,
      height = cfg$plot$qc_scatter_height %||% 5)

  invisible(NULL)
}