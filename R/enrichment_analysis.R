#' Run escape pathway enrichment analysis
#' @param seurat_obj Seurat object
#' @param species Character. "Mouse" or "Human"
#' @param library Character. e.g., "H" for Hallmark, "C5" for GO
#' @param method Character. e.g., "ssGSEA", "UCell", "AUCell"
#' @param min_size Integer. Min genes in a pathway.
#' @param gene_sets Optional named list of pathways; if NULL, load from escape.
#' @param gene_set_cache Optional environment used to reuse loaded gene sets.
#' @return Seurat object with escape assay added
#' @export
run_escape_enrichment <- function(seurat_obj, species = "Mouse", library = "H",
                                  method = "ssGSEA", min_size = 5,
                                  gene_sets = NULL, gene_set_cache = NULL) {
  if (!requireNamespace("escape", quietly = TRUE)) {
    message("   [SKIP Escape] Package 'escape' is not installed.")
    return(seurat_obj)
  }

  message("--- Running escape pathway enrichment (", method, ") ---")
  if (is.null(gene_sets)) {
    gene_sets <- .get_escape_gene_sets(species, library, gene_set_cache)
  }
  if (is.null(gene_sets) || length(gene_sets) == 0L) {
    message("   [SKIP Escape] No gene sets were returned.")
    return(seurat_obj)
  }

  gene_sets <- lapply(gene_sets, unique)
  expressed_genes <- rownames(SeuratObject::LayerData(
    seurat_obj, assay = "RNA", layer = "counts"
  ))
  gene_sets <- lapply(gene_sets, intersect, y = expressed_genes)
  gene_sets <- gene_sets[lengths(gene_sets) >= min_size]
  if (length(gene_sets) == 0L) {
    message("   [SKIP Escape] No ", library,
            " gene sets retain at least ", min_size,
            " genes present in the RNA assay.")
    return(seurat_obj)
  }

  assay_name <- paste0("escape.", method)
  original_obj <- seurat_obj
  if (assay_name %in% names(seurat_obj@assays)) {
    message("   [INFO Escape] Replacing existing assay '", assay_name, "'.")
    seurat_obj[[assay_name]] <- NULL
  }
  escape_args <- list(
    input.data = seurat_obj,
    method = method,
    new.assay.name = assay_name,
    gene.sets = gene_sets,
    min.size = min_size
  )
  if (identical(tolower(method), "ucell")) {
    escape_args$maxRank <- max(1500L, length(expressed_genes))
  }
  seurat_obj <- safe_run(
    do.call(escape::runEscape, escape_args),
    label = "escape::runEscape",
    fallback = original_obj
  )
  if (assay_name %in% names(seurat_obj@assays)) {
    message("   -> Enrichment successful. Assay '", assay_name, "' added.")
  } else {
    message("   [SKIP Escape] Enrichment did not create assay '", assay_name, "'.")
  }
  seurat_obj
}

.get_escape_gene_sets <- function(species, library, gene_set_cache = NULL) {
  cache_key <- paste(tolower(species), paste(library, collapse = ","), sep = "::")
  if (is.environment(gene_set_cache) &&
      exists(cache_key, envir = gene_set_cache, inherits = FALSE)) {
    return(get(cache_key, envir = gene_set_cache, inherits = FALSE))
  }

  if (!requireNamespace("escape", quietly = TRUE)) return(NULL)
  if (!requireNamespace("msigdb", quietly = TRUE)) {
    message("   [SKIP Escape] Package 'msigdb' is not installed; no gene sets loaded.")
    if (is.environment(gene_set_cache)) assign(cache_key, NULL, envir = gene_set_cache)
    return(NULL)
  }

  esp_species <- if (tolower(species) == "mouse") "Mus musculus" else "Homo sapiens"
  gene_sets <- safe_run(
    escape::getGeneSets(species = esp_species, library = library),
    label = "escape::getGeneSets"
  )
  if (is.environment(gene_set_cache)) {
    assign(cache_key, gene_sets, envir = gene_set_cache)
  }
  gene_sets
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

  score_matrix <- tryCatch(
    SeuratObject::LayerData(seurat_obj, assay = assay_name, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(score_matrix)) {
    message("   [SKIP Escape plot] Score data are not available in assay '",
            assay_name, "'.")
    return(invisible(NULL))
  }
  n_groups <- length(unique(stats::na.omit(
    as.character(seurat_obj@meta.data[[group_col]])
  )))

  escape_dir <- file.path(out_dir, "Escape")
  make_dir(escape_dir)
  p <- safe_run(
    escape::heatmapEnrichment(seurat_obj, group.by = group_col,
                              assay = assay_name, scale = TRUE,
                              cluster.rows = nrow(score_matrix) > 1L,
                              cluster.columns = n_groups > 1L) +
      ggplot2::ggtitle(paste0("Enrichment (", method, ") | ", group_col)),
    label = paste0("escape::heatmapEnrichment ", group_col)
  )
  if (!is.null(p)) {
    plot_path <- file.path(
      escape_dir, paste0("Escape_Heatmap_", prefix, "_", group_col, ".png")
    )
    save_png(p, plot_path, width = 10, height = 8)
    return(invisible(plot_path))
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

.get_gsea_pathways <- function(species, collection, pathway_cache = NULL,
                               pathway_loader = .hallmark_gene_sets) {
  cache_key <- paste(tolower(species), tolower(collection), sep = "::")
  if (is.environment(pathway_cache) &&
      exists(cache_key, envir = pathway_cache, inherits = FALSE)) {
    return(get(cache_key, envir = pathway_cache, inherits = FALSE))
  }

  pathways <- tryCatch(pathway_loader(species, collection), error = function(e) {
    message("   [SKIP GSEA] Could not load MSigDB collection '", collection,
            "' gene sets: ", conditionMessage(e))
    NULL
  })
  if (is.environment(pathway_cache)) {
    assign(cache_key, pathways, envir = pathway_cache)
  }
  pathways
}
