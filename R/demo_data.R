#' Small synthetic Seurat object (for demos, examples and tests)
#'
#' Builds a seeded single-cell dataset with three planted cell populations
#' (T-like, B-like, monocyte-like), four samples in two conditions, a PCA, a
#' UMAP and clusters. Nothing here is real biology: it exists so the report and
#' the explorer can be tried -- and tested -- without running the pipeline.
#'
#' @param n_cells Number of cells.
#' @param n_genes Number of genes (at least 40).
#' @param seed Random seed.
#' @param embed Run PCA / neighbours / clustering / UMAP (set \code{FALSE} for a
#'   bare counts object).
#' @return A Seurat object.
#' @export
#' @examples
#' \dontrun{
#' obj <- sparks_toy_object(300)
#' table(obj$seurat_clusters)
#' }
sparks_toy_object <- function(n_cells = 600L, n_genes = 300L, seed = 1L,
                              embed = TRUE) {
  for (pkg in c("Seurat", "SeuratObject", "Matrix")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("sparks_toy_object() requires package '", pkg, "'.", call. = FALSE)
    }
  }
  n_genes <- max(as.integer(n_genes), 40L)
  n_cells <- max(as.integer(n_cells), 30L)
  withr_seed <- function(seed, code) {
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    force(code)
  }
  withr_seed(seed, {
    markers <- list(
      "T cells"   = c("CD3E", "CD3D", "IL7R", "CD8A", "CCR7", "LTB"),
      "B cells"   = c("MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD74"),
      "Monocytes" = c("LYZ", "CST3", "S100A8", "S100A9", "CD14", "FCN1")
    )
    mt <- c("MT-ND1", "MT-CO1", "MT-CYB")
    named <- c(unlist(markers, use.names = FALSE), mt)
    filler <- sprintf("GENE%04d", seq_len(max(0L, n_genes - length(named))))
    genes <- c(named, filler)
    pop <- sample(names(markers), n_cells, replace = TRUE, prob = c(0.45, 0.30, 0.25))
    sample_id <- sample(paste0("S", 1:4), n_cells, replace = TRUE)
    condition <- ifelse(sample_id %in% c("S1", "S2"), "WT", "KO")
    base <- stats::rlnorm(length(genes), meanlog = -1, sdlog = 1)
    lambda <- matrix(base, nrow = length(genes), ncol = n_cells)
    for (p in names(markers)) {
      rows <- match(markers[[p]], genes)
      lambda[rows, pop == p] <- lambda[rows, pop == p] * 12
    }
    # A mild condition effect so DE between WT / KO has something to find.
    ko_rows <- match(c("CD74", "LYZ"), genes)
    lambda[ko_rows, condition == "KO"] <- lambda[ko_rows, condition == "KO"] * 2.5
    counts <- matrix(stats::rpois(length(lambda), lambda), nrow = length(genes),
                     dimnames = list(genes, sprintf("CELL%04d", seq_len(n_cells))))
    counts <- methods::as(counts, "CsparseMatrix")
    obj <- Seurat::CreateSeuratObject(counts, project = "SPARKS_toy",
                                      min.cells = 0, min.features = 0)
    obj$orig.ident <- factor(sample_id)
    SeuratObject::Idents(obj) <- "orig.ident"
    obj$condition <- factor(condition, levels = c("WT", "KO"))
    obj$true_population <- factor(pop)
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = min(100L, n_genes - 1L),
                                        verbose = FALSE)
    if (isTRUE(embed)) obj <- .toy_embed(obj, seed)
    obj
  })
}

.toy_embed <- function(obj, seed = 1L) {
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  npcs <- min(10L, ncol(obj) - 1L, nrow(obj) - 1L)
  obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  dims <- seq_len(min(8L, npcs))
  obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 0.4, verbose = FALSE)
  if (requireNamespace("uwot", quietly = TRUE)) {
    obj <- Seurat::RunUMAP(obj, dims = dims, seed.use = seed, verbose = FALSE,
                           n.neighbors = min(15L, ncol(obj) - 1L))
  } else {
    # No uwot: fall back to the first two PCs so the object still has a 2-D embedding.
    pcs <- Seurat::Embeddings(obj, "pca")[, 1:2, drop = FALSE]
    colnames(pcs) <- c("umap_1", "umap_2")
    obj[["umap"]] <- SeuratObject::CreateDimReducObject(pcs, key = "umap_",
                                                        assay = "RNA")
  }
  obj
}

#' Write a synthetic SPARKS results directory
#'
#' Produces a results tree with the same layout, file names and columns as a
#' real pipeline run (QC counts, cluster markers, DEG tables, proportions,
#' pseudobulk, GSEA, merged object) built from [sparks_toy_object()], so
#' [sparks_report()] and [sparks_explorer()] can be tried immediately. A
#' \code{DEMO_DATA.txt} marker file flags the directory as synthetic and the
#' report shows a "demo data" badge.
#'
#' @param dir Directory to create / fill.
#' @param n_cells,n_genes,seed Passed to [sparks_toy_object()].
#' @param comparison Name of the comparison group (top-level folder).
#' @return \code{dir}, invisibly.
#' @export
#' @examples
#' \dontrun{
#' demo <- tempfile("sparks_demo_")
#' sparks_demo_results(demo)
#' sparks_report(demo)
#' }
sparks_demo_results <- function(dir, n_cells = 600L, n_genes = 300L, seed = 1L,
                                comparison = "Demo") {
  make_dir(dir)
  obj <- sparks_toy_object(n_cells, n_genes, seed)
  obj$cell_type <- factor(
    c("T cells", "B cells", "Monocytes")[match(as.character(obj$true_population),
                                               c("T cells", "B cells", "Monocytes"))],
    levels = c("T cells", "B cells", "Monocytes"))
  unit <- "All_cells"
  unit_dir <- file.path(dir, comparison, comparison, unit)
  qc_dir <- file.path(dir, comparison, "QC")
  rdata_dir <- file.path(dir, comparison, "RData")
  for (d in c(qc_dir, rdata_dir)) make_dir(d)

  # QC: per-sample cell counts before / after filtering, doublet summaries.
  samples <- levels(obj$orig.ident)
  for (s in samples) {
    n_after <- sum(obj$orig.ident == s)
    n_before <- n_after + as.integer(stats::runif(1, 8, 30))
    n_dbl <- as.integer(stats::runif(1, 3, 12))
    writeLines(c(paste("Total cells:", n_before), "\nCells by orig.ident:",
                 paste0(s, ": ", n_before)),
               file.path(qc_dir, paste0("before_filtering_", s, "_cell_counts.txt")))
    writeLines(c(paste("Total cells:", n_after), "\nCells by orig.ident:",
                 paste0(s, ": ", n_after)),
               file.path(qc_dir, paste0("after_filtering_", s, "_cell_counts.txt")))
    writeLines(c("", "doublet singlet", paste(n_dbl, n_after)),
               file.path(qc_dir, paste0("DoubletStats_", s, ".txt")))
    writeLines(paste("Removed by scDblFinder:", n_dbl),
               file.path(qc_dir, paste0("Doublet_Removed_", s, ".txt")))
  }
  writeLines(c(paste("Total cells:", ncol(obj)), "\nCells by orig.ident:",
               paste0(names(table(obj$orig.ident)), ": ", as.integer(table(obj$orig.ident)))),
             file.path(qc_dir, paste0("after_clustering_", comparison, "_cell_counts.txt")))

  for (grouping in c("seurat_clusters", "cell_type")) {
    gdir <- file.path(unit_dir, grouping)
    deg_dir <- file.path(gdir, "DEG")
    for (d in file.path(gdir, c("DEG", "Proportions", "Pseudobulk", "GSEA", "Expression"))) make_dir(d)
    Seurat::Idents(obj) <- grouping
    markers <- suppressWarnings(Seurat::FindAllMarkers(
      obj, only.pos = FALSE, min.pct = 0.05, logfc.threshold = 0.1, verbose = FALSE))
    if (is.null(markers) || !nrow(markers)) next
    write_table(markers, file.path(deg_dir, paste0("AllMarkers_", unit, "_by_", grouping, ".txt")))
    sig <- markers[!is.na(markers$p_val_adj) & markers$p_val_adj < 0.05, , drop = FALSE]
    counts <- data.frame(Group = levels(droplevels(factor(markers$cluster))),
                         stringsAsFactors = FALSE)
    counts$DEG_Count <- vapply(counts$Group, function(g) sum(sig$cluster == g), integer(1))
    write_table(counts, file.path(deg_dir, paste0("DEG_Counts_", unit, "_by_", grouping, ".txt")))
    for (g in unique(as.character(sig$cluster))) {
      write_table(sig[sig$cluster == g, , drop = FALSE],
                  file.path(deg_dir, paste0("DEG_", unit, "_", .safe_filename(g), "_by_", grouping, ".txt")))
    }
    # Cluster markers (positive only) like generate_cluster_markers().
    pos <- markers[markers$avg_log2FC > 0, , drop = FALSE]
    write_table(pos, file.path(deg_dir, paste0("ClusterMarkers_", unit, "_", grouping, ".txt")))

    # Proportions by sample and by condition.
    for (by in c("orig.ident", "condition")) {
      tbl <- table(obj@meta.data[[grouping]], obj@meta.data[[by]])
      prop <- as.data.frame(prop.table(tbl, margin = 2L))
      prop$Count <- as.integer(tbl)
      names(prop)[1:2] <- c("Var1", "Var2")
      write_table(prop, file.path(gdir, "Proportions",
                                  paste0("Proportions_", unit, "_", by, "_", grouping, ".txt")))
    }

    # Pseudobulk counts + a log-fold-change table (illustrative, not edgeR).
    pb <- Seurat::AggregateExpression(obj, group.by = c(grouping, "orig.ident"),
                                      assays = "RNA", return.seurat = FALSE, verbose = FALSE)$RNA
    pb <- as.data.frame(as.matrix(pb))
    pb$gene <- rownames(pb)
    pb <- pb[, c("gene", setdiff(colnames(pb), "gene"))]
    write_table(pb, file.path(gdir, "Pseudobulk",
                              paste0("PseudobulkCounts_", unit, "_by_", grouping, ".txt")))

    # Toy GSEA table: a "pathway" per planted population, NES from marker logFC.
    gsea <- do.call(rbind, lapply(levels(droplevels(factor(markers$cluster))), function(g) {
      m <- markers[markers$cluster == g, , drop = FALSE]
      m <- m[order(-m$avg_log2FC), , drop = FALSE]
      paths <- c("T_CELL_SIGNALING", "B_CELL_RECEPTOR", "MONOCYTE_ACTIVATION", "OXPHOS", "CELL_CYCLE")
      data.frame(
        pathway = paths,
        pval = stats::runif(length(paths), 1e-6, 0.2),
        padj = stats::runif(length(paths), 1e-4, 0.3),
        ES = stats::rnorm(length(paths), 0.4, 0.3),
        NES = round(stats::rnorm(length(paths), 0, 1.6), 3L),
        size = sample(15:80, length(paths)),
        leadingEdge = vapply(seq_along(paths), function(i)
          paste(utils::head(m$gene, 4L), collapse = ","), character(1)),
        group = g,
        stringsAsFactors = FALSE
      )
    }))
    write_table(gsea, file.path(gdir, "GSEA",
                                paste0("GSEA_toy_", unit, "_", grouping, ".tsv")))
  }

  saveRDS(obj, file.path(rdata_dir, paste0("Merged_", comparison, ".rds")))
  writeLines(c("SPARKS synthetic demo data.",
               "Generated by sparks_demo_results(); not real biology."),
             file.path(dir, "DEMO_DATA.txt"))
  invisible(dir)
}
