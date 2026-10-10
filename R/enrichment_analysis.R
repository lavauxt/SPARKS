#' Run escape pathway enrichment analysis
#' @param seurat_obj Seurat object
#' @param species Character. "Mouse" or "Human"
#' @param library Character. e.g., "H" for Hallmark, "C5" for GO
#' @param method Character. e.g., "ssGSEA", "UCell", "AUCell"
#' @param min_size Integer. Min genes in a pathway.
#' @return Seurat object with escape assay added
#' @export
run_escape_enrichment <- function(seurat_obj, species = "Mouse", library = "H",
                                  method = "ssGSEA", min_size = 5) {
  if (!requireNamespace("escape", quietly = TRUE)) {
    message("   [SKIP Escape] Package 'escape' is not installed.")
    return(seurat_obj)
  }
  if (!requireNamespace("msigdb", quietly = TRUE)) {
    message("   [SKIP Escape] Package 'msigdb' is not installed; no gene sets loaded.")
    return(seurat_obj)
  }

  message("--- Running escape pathway enrichment (", method, ") ---")
  esp_species <- if (tolower(species) == "mouse") "Mus musculus" else "Homo sapiens"
  gene_sets <- safe_run(
    escape::getGeneSets(species = esp_species, library = library),
    label = "escape::getGeneSets"
  )
  if (is.null(gene_sets) || length(gene_sets) == 0L) {
    message("   [SKIP Escape] No gene sets were returned.")
    return(seurat_obj)
  }

  assay_name <- paste0("escape.", method)
  seurat_obj <- safe_run(
    escape::runEscape(seurat_obj, method = method,
                      new.assay.name = assay_name,
                      gene.sets = gene_sets, min.size = min_size),
    label = "escape::runEscape",
    fallback = seurat_obj
  )
  if (assay_name %in% names(seurat_obj@assays)) {
    message("   -> Enrichment successful. Assay '", assay_name, "' added.")
  } else {
    message("   [SKIP Escape] Enrichment did not create assay '", assay_name, "'.")
  }
  seurat_obj
}

#' Generate and save escape heatmaps
#' @param seurat_obj Seurat object
#' @param method Character
#' @param group_col Character
#' @param out_dir Character
#' @param prefix Character
#' @export
generate_escape_plots <- function(seurat_obj, method = "ssGSEA", group_col,
                                  out_dir, prefix) {
  if (!requireNamespace("escape", quietly = TRUE)) {
    message("   [SKIP Escape plot] Package 'escape' is not installed.")
    return(invisible(NULL))
  }
  assay_name <- paste0("escape.", method)
  if (!assay_name %in% names(seurat_obj@assays)) {
    message("   [SKIP Escape plot] Assay '", assay_name, "' is not present.")
    return(invisible(NULL))
  }
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    message("   [SKIP Escape plot] Grouping metadata '", group_col, "' is missing.")
    return(invisible(NULL))
  }

  escape_dir <- file.path(out_dir, "Escape")
  make_dir(escape_dir)
  p <- safe_run(
    escape::heatmapEnrichment(seurat_obj, group.by = group_col,
                              assay = assay_name, scale = TRUE,
                              cluster.rows = TRUE, cluster.columns = TRUE) +
      ggplot2::ggtitle(paste0("Enrichment (", method, ") | ", group_col)),
    label = paste0("escape::heatmapEnrichment ", group_col)
  )
  if (!is.null(p)) {
    save_png(p, file.path(escape_dir,
                          paste0("Escape_Heatmap_", prefix, "_", group_col, ".png")),
             width = 10, height = 8)
  }
  invisible(NULL)
}

.hallmark_gene_sets <- function(species = "Mouse", collection = "H") {
  if (!requireNamespace("msigdb", quietly = TRUE) ||
      !requireNamespace("GSEABase", quietly = TRUE)) {
    stop("Hallmark GSEA requires 'msigdb' and 'GSEABase'.", call. = FALSE)
  }
  org <- if (tolower(species) == "mouse") "mm" else "hs"
  sets <- msigdb::getMsigdb(org = org, id = "SYM")
  collection <- tolower(collection)
  sets <- msigdb::subsetCollection(sets, collection = collection)
  pathways <- lapply(as.list(sets), GSEABase::geneIds)
  names(pathways) <- vapply(as.list(sets), GSEABase::setName, character(1))
  pathways[lengths(pathways) > 0L]
}
