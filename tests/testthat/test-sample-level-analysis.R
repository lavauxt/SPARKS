test_that("pseudobulk aggregation sums raw counts by group and sample", {
  counts <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:6))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$cell_group <- rep(c("A", "B"), each = 3)
  object$sample <- rep(c("s1", "s2", "s3"), times = 2)
  object$condition <- rep(c("Control", "Treatment", "Control"), times = 2)

  aggregated <- SPARKS:::.aggregate_group_counts(
    object, "cell_group", "A", "sample", "condition"
  )

  expect_identical(colnames(aggregated$counts), aggregated$sample_meta$sample)
  expect_equal(as.numeric(aggregated$counts["g1", ]), c(1, 3, 5))
  expect_equal(aggregated$sample_meta$condition,
               c("Control", "Treatment", "Control"))

  object$cell_group <- c("A", "A", "B", "B", "B", "B")
  only_group_a <- SPARKS:::.aggregate_group_counts(
    object, "cell_group", "A", "sample", "condition"
  )
  expect_setequal(colnames(only_group_a$counts), c("s1", "s2"))
  expect_equal(ncol(only_group_a$counts), 2L)
})

test_that("edgeR pseudobulk fit checks condition and replicate design", {
  counts <- matrix(
    c(12, 21, 15, 19, 55, 60, 58, 63,
      8, 13, 10, 16, 35, 41, 38, 45),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("s", 1:4))
  )
  sample_meta <- data.frame(
    sample = colnames(counts),
    condition = c("A", "A", "B", "B")
  )
  if (requireNamespace("edgeR", quietly = TRUE)) {
    fit <- SPARKS:::.fit_pseudobulk_edger(counts, sample_meta, min_replicates = 2L)
    expect_identical(fit$status, "ok")
    expect_true(all(c("gene", "logFC", "p_val_adj", "rank_stat") %in% names(fit$result)))
    expect_identical(unique(fit$result$condition_1), "A")
    expect_identical(unique(fit$result$condition_2), "B")
  } else {
    fit <- SPARKS:::.fit_pseudobulk_edger(counts, sample_meta, min_replicates = 2L)
    expect_match(fit$status, "edgeR is not installed")
  }

  insufficient <- SPARKS:::.fit_pseudobulk_edger(
    counts[, 1:3], sample_meta[1:3, ], min_replicates = 2L
  )
  if (requireNamespace("edgeR", quietly = TRUE)) {
    expect_match(insufficient$status, "at least 2 samples per condition")
  } else {
    expect_match(insufficient$status, "edgeR is not installed")
  }
})

test_that("sample-level analysis writes Hallmark GSEA and pseudobulk DEG separately", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("fgsea")

  set.seed(2491)
  counts <- matrix(rpois(400, lambda = 25), nrow = 50, ncol = 8)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$cell_group <- "GroupA"
  object$sample <- rep(paste0("s", 1:4), each = 2)
  object$condition <- rep(c("A", "B", "A", "B"), each = 2)

  out <- tempfile("sample-level-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  pathways <- list(test_set = paste0("g", 1:15))
  collection <- "Test"

  gsea_files <- SPARKS:::.run_sample_level_analysis(
    object, "cell_group", file.path(out, "Pseudobulk", "DEG"),
    file.path(out, "GSEA"), "Main", "Mouse",
    run_pseudobulk_deg = TRUE, run_gsea = TRUE, pathways = pathways,
    collection = collection,
    gsea_min_size = 5L, gsea_max_size = 20L
  )

  expect_length(list.files(file.path(out, "Pseudobulk", "DEG"),
                           pattern = "^PseudobulkDEG_"), 1L)
  expect_length(list.files(file.path(out, "GSEA"),
                           pattern = "^GSEA_Test_"), 1L)
  expect_length(gsea_files, 1L)
  gsea_table <- utils::read.delim(gsea_files, check.names = FALSE)
  expect_true(all(c("pathway", "group", "NES", "padj") %in% names(gsea_table)))
  expect_identical(utils::head(gsea_table$group, 1L), "GroupA")
})

test_that("GSEA writes a status explanation when replicate design is insufficient", {
  skip_if_not_installed("edgeR")

  counts <- matrix(
    c(10, 12, 14, 16, 20, 21, 19, 22, 11, 13, 15, 17, 18, 20, 23, 24), nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", seq_len(8)))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$group <- "A"
  object$sample <- rep(c("s1", "s2"), each = 4)
  object$condition <- rep(c("WT", "KO"), each = 4)
  out <- tempfile("gsea-status-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  SPARKS:::.run_sample_level_analysis(
    object, "group", file.path(out, "DEG"), file.path(out, "GSEA"),
    "Main", "Mouse", run_gsea = TRUE, pathways = list(test = c("g1", "g2")),
    gsea_min_size = 1L, gsea_max_size = 1L
  )

  status_file <- list.files(out, pattern = "^GSEA_Status_",
                            full.names = TRUE)
  expect_length(status_file, 1L)
  expect_false(dir.exists(file.path(out, "GSEA")))
  expect_match(readLines(status_file), "at least 2 samples per condition")
  snapshot <- SPARKS:::.snapshot_result_files(out)
  expect_true(any(grepl("GSEA_Status_", rownames(snapshot), fixed = TRUE)))
})

test_that("GSEA does not load MSigDB when replicate design is insufficient", {
  skip_if_not_installed("edgeR")

  counts <- matrix(
    c(10, 12, 14, 16, 20, 21, 19, 22, 11, 13, 15, 17, 18, 20, 23, 24), nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", seq_len(8)))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$group <- "A"
  object$sample <- rep(c("s1", "s2"), each = 4)
  object$condition <- rep(c("WT", "KO"), each = 4)
  loader_calls <- 0L

  expect_message(
    SPARKS:::.run_sample_level_analysis(
      object, "group", tempfile("pseudobulk-"), tempfile("gsea-"),
      "Main", "Mouse", run_gsea = TRUE,
      pathway_loader = function(species, collection) {
        loader_calls <<- loader_calls + 1L
        list(test = c("g1", "g2"))
      }
    ),
    "at least 2 samples per condition"
  )
  expect_identical(loader_calls, 0L)
})

test_that("GSEA pathway cache loads each collection once", {
  cache <- new.env(parent = emptyenv())
  loader_calls <- 0L
  loader <- function(species, collection) {
    loader_calls <<- loader_calls + 1L
    list(test = c("g1", "g2"))
  }

  expect_identical(SPARKS:::.get_gsea_pathways("Mouse", "H", cache, loader),
                   list(test = c("g1", "g2")))
  expect_identical(SPARKS:::.get_gsea_pathways("Mouse", "H", cache, loader),
                   list(test = c("g1", "g2")))
  expect_identical(loader_calls, 1L)

  failed_cache <- new.env(parent = emptyenv())
  fail_loader <- function(species, collection) {
    loader_calls <<- loader_calls + 1L
    stop("unavailable")
  }
  expect_message(
    SPARKS:::.get_gsea_pathways("Mouse", "H", failed_cache, fail_loader),
    "Could not load MSigDB collection"
  )
  expect_message(
    pathways <- SPARKS:::.get_gsea_pathways(
      "Mouse", "H", failed_cache, fail_loader
    ),
    NA
  )
  expect_null(pathways)
  expect_identical(loader_calls, 2L)
})

test_that("disabled GSEA does not load a pathway collection", {
  loader_calls <- 0L
  expect_null(SPARKS:::.run_sample_level_analysis(
    seurat_obj = NULL, group_col = "group", pseudobulk_dir = tempfile(),
    gsea_dir = tempfile(), file_prefix = "Main", species = "Mouse",
    run_pseudobulk_deg = FALSE, run_gsea = FALSE,
    pathway_loader = function(species, collection) {
      loader_calls <<- loader_calls + 1L
      list(test = "g1")
    }
  ))
  expect_identical(loader_calls, 0L)
})

test_that("GSEA skips before loading MSigDB when fgsea is unavailable", {
  skip_if_not_installed("edgeR")
  skip_if(requireNamespace("fgsea", quietly = TRUE))

  counts <- matrix(
    c(10, 12, 14, 16, 20, 21, 19, 22, 11, 13, 15, 17, 18, 20, 23, 24), nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", seq_len(8)))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$group <- "A"
  object$sample <- rep(paste0("s", seq_len(4)), each = 2)
  object$condition <- rep(c("WT", "KO"), each = 4)
  loader_calls <- 0L

  SPARKS:::.run_sample_level_analysis(
    object, "group", tempfile("pseudobulk-"), tempfile("gsea-"),
    "Main", "Mouse", run_gsea = TRUE,
    pathway_loader = function(species, collection) {
      loader_calls <<- loader_calls + 1L
      list(test = c("g1", "g2"))
    }
  )
  expect_identical(loader_calls, 0L)
})

test_that("GSEA loads its pathway collection only once per cache", {
  skip_if_not_installed("edgeR")
  skip_if_not_installed("fgsea")

  set.seed(6274)
  counts <- matrix(stats::rpois(800, lambda = 25), nrow = 50, ncol = 16)
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$group <- rep(c("A", "B"), each = 8)
  object$sample <- c(rep(paste0("A", 1:4), each = 2),
                     rep(paste0("B", 1:4), each = 2))
  object$condition <- rep(c(rep("Control", 2), rep("Treatment", 2)), times = 2)
  out <- tempfile("gsea-cache-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  pathway_cache <- new.env(parent = emptyenv())
  loader_calls <- 0L
  pathway_loader <- function(species, collection) {
    loader_calls <<- loader_calls + 1L
    list(test_set = paste0("g", 1:15))
  }

  for (group in c("A", "B")) {
    group_object <- subset(object, cells = colnames(object)[object$group == group])
    group_object$group <- "group"
    SPARKS:::.run_sample_level_analysis(
      group_object, "group", file.path(out, "DEG"), file.path(out, "GSEA"),
      "Main", "Mouse", run_gsea = TRUE, pathway_cache = pathway_cache,
      pathway_loader = pathway_loader, gsea_min_size = 5L, gsea_max_size = 20L
    )
  }

  expect_identical(loader_calls, 1L)
})

test_that("GSEA projection maps pathway NES to cells by group", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:24, nrow = 2,
      dimnames = list(c("gene-a", "gene-b"), paste0("cell", 1:12))
    )
  ))
  object$cell_group <- rep(c("A", "B"), each = 6)
  emb <- matrix(seq_len(24), ncol = 2,
                dimnames = list(colnames(object), c("UMAP_1", "UMAP_2")))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = "RNA"
  )
  gsea_dir <- tempfile("gsea-projection-")
  dir.create(gsea_dir)
  on.exit(unlink(gsea_dir, recursive = TRUE), add = TRUE)
  utils::write.table(
    data.frame(pathway = c("HALLMARK_A", "HALLMARK_A", "HALLMARK_B"),
               group = c("A", "B", "A"), NES = c(1.5, -1.2, 0.5),
               padj = c(0.01, 0.02, 0.2)),
    file.path(gsea_dir, "GSEA_H_Main.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )

  paths <- SPARKS:::.generate_gsea_projection(
    object, gsea_dir, group_by_col = "cell_group", reduction = "umap",
    top_n = 1L
  )

  expect_length(paths, 1L)
  expect_match(basename(paths), "HALLMARK_A")
  expect_true(file.exists(paths))
  expect_length(SPARKS:::.generate_gsea_projection(
    object, gsea_dir, group_by_col = "missing", reduction = "umap"
  ), 0L)

  projection_files <- SPARKS:::.snapshot_result_files(gsea_dir)
  expect_true(any(grepl("GSEA_Projection_", rownames(projection_files), fixed = TRUE)))
})

test_that("GSEA projection skips when embedding is missing", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(1:24, nrow = 2,
                    dimnames = list(c("g1", "g2"), paste0("cell", 1:12)))
  ))
  object$cell_group <- rep(c("A", "B"), each = 6)
  gsea_dir <- tempfile("gsea-no-embedding-")
  dir.create(gsea_dir)
  on.exit(unlink(gsea_dir, recursive = TRUE), add = TRUE)
  utils::write.table(
    data.frame(pathway = "HALLMARK_A", group = "A", NES = 1.5, padj = 0.01),
    file.path(gsea_dir, "GSEA_H_Main.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )

  expect_message(
    paths <- SPARKS:::.generate_gsea_projection(
      object, gsea_dir, group_by_col = "cell_group", reduction = "umap"
    ),
    "Reduction 'umap' is not available"
  )
  expect_length(paths, 0L)
})

test_that("FeaturePlots and violins use RNA genes and split by condition", {
  counts <- matrix(
    c(1:12, 12:1), nrow = 2,
    dimnames = list(c("rna-only", "shared"), paste0("cell", seq_len(12)))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  object[["SCT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(counts["shared", , drop = FALSE], sparse = TRUE),
    data = Matrix::Matrix(counts["shared", , drop = FALSE], sparse = TRUE)
  )
  object$condition <- rep(c("WT", "KO"), each = 6)
  object$group <- rep(c("A", "B"), each = 6)
  emb <- matrix(seq_len(24), ncol = 2,
                dimnames = list(colnames(object), c("UMAP_1", "UMAP_2")))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = "RNA"
  )
  Seurat::DefaultAssay(object) <- "SCT"

  feature_dir <- tempfile("feature-split-")
  violin_dir <- tempfile("violin-split-")
  on.exit(unlink(c(feature_dir, violin_dir), recursive = TRUE), add = TRUE)

  SPARKS::generate_feature_plots(
    object, genes = c("rna-only", "sct-absent"), out_dir = feature_dir,
    prefix = "Main", reduction = "umap", assay = "RNA", split_by = "condition",
    keep_scale = "all"
  )
  SPARKS::generate_violin_plots(
    object, genes = "rna-only", out_dir = violin_dir, prefix = "Main",
    group_by_col = "group", split_by = "condition", assay = "RNA"
  )

  expect_true(file.exists(file.path(feature_dir, "Main_FeaturePlot_rna-only.png")))
  expect_false(file.exists(file.path(feature_dir, "Main_FeaturePlot_sct-absent.png")))
  expect_true(file.exists(file.path(violin_dir, "VlnPlot_Main_rna-only.png")))
})

test_that("split FeaturePlots retain a shared visible colour scale", {
  counts <- matrix(
    c(1:12, 12:1), nrow = 2,
    dimnames = list(c("gene-a", "gene-b"), paste0("cell", seq_len(12)))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object <- Seurat::NormalizeData(object, assay = "RNA", verbose = FALSE)
  object$condition <- factor(rep(c("WT", "KO"), each = 6), levels = c("WT", "KO"))
  emb <- matrix(seq_len(24), ncol = 2,
                dimnames = list(colnames(object), c("UMAP_1", "UMAP_2")))
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = emb, key = "UMAP_", assay = "RNA"
  )

  plot <- SPARKS:::.make_feature_plot(
    object, "gene-a", "umap", "RNA", "condition", "all"
  )
  expect_s3_class(plot, "patchwork")
  expect_s3_class(plot$patches$plots[[1]]$scales$get_scales("colour"),
                  "ScaleContinuous")
  expect_length(plot$patches$plots, 2L)
})

test_that("analysis output folders distinguish analysis types", {
  root <- tempfile("analysis-folders-")
  dirs <- SPARKS:::.make_analysis_dirs(root)

  created_dirs <- dirs[names(dirs) != "GSEA"]
  expect_true(all(vapply(created_dirs, dir.exists, logical(1))))
  expect_false(dir.exists(dirs$GSEA))
  expect_false(identical(dirs$Correlation, dirs$DEG))
  expect_false(identical(dirs$Pseudobulk, dirs$DEG))
})

test_that("escape heatmaps handle a single pathway and one grouping level", {
  skip_if_not_installed("escape")

  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:24, nrow = 2,
      dimnames = list(c("pathway-a", "pathway-b"), paste0("cell", 1:12))
    )
  ))
  suppressWarnings(object[["escape.UCell"]] <- SeuratObject::CreateAssay5Object(
    data = Matrix::Matrix(
      matrix(1:12, nrow = 1,
             dimnames = list("pathway-a", paste0("cell", 1:12))),
      sparse = TRUE
    )
  ))
  object$one_group <- "A"
  out <- tempfile("escape-plots-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)

  path <- SPARKS::generate_escape_plots(
    object, method = "UCell", group_col = "one_group", out_dir = out,
    prefix = "test"
  )
  expect_true(file.exists(path))
  expect_true(file.exists(file.path(
    out, "Escape", "Escape_Scores_UCell_test_one_group.tsv"
  )))
})

test_that("Hallmark escape scoring filters to detected genes for UCell and AUCell", {
  skip_if_not_installed("escape")
  skip_if_not_installed("msigdb")

  pathways <- SPARKS:::.hallmark_gene_sets("Human", "H")
  genes <- unique(c(pathways[[1]], paste0("TEST", 1:270)))
  set.seed(6571)
  counts <- matrix(stats::rpois(length(genes) * 2, lambda = 4),
                   nrow = length(genes),
                   dimnames = list(genes, c("WT_cell", "KO_cell")))
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$condition <- c("WT", "KO")

  for (method in c("UCell", "AUCell")) {
    scored <- SPARKS::run_escape_enrichment(
      object, species = "Human", library = "H", method = method,
      min_size = 5L, gene_sets = pathways
    )
    expect_true(paste0("escape.", method) %in% names(scored@assays))
    scores <- SeuratObject::LayerData(
      scored, assay = paste0("escape.", method), layer = "data"
    )
    expect_identical(colnames(scores), colnames(object))
    expected_pathways <- pathways[lengths(lapply(pathways, intersect, genes)) >= 5L]
    expect_equal(nrow(scores), length(expected_pathways))
  }
})

test_that("escape gene sets are cached across scoring methods", {
  skip_if_not_installed("escape")
  skip_if_not_installed("msigdb")
  cache <- new.env(parent = emptyenv())

  sets <- SPARKS:::.get_escape_gene_sets("Human", "H", cache)
  expect_type(sets, "list")
  expect_identical(
    SPARKS:::.get_escape_gene_sets("Human", "H", cache),
    sets
  )
  expect_true(exists("human::H", envir = cache, inherits = FALSE))
})

test_that("report file tracking excludes stale result tables", {
  root <- tempfile("result-snapshot-")
  dir.create(file.path(root, "DEG"), recursive = TRUE)
  stale <- file.path(root, "DEG", "AllMarkers_old.txt")
  writeLines("old", stale)
  snapshot <- SPARKS:::.snapshot_result_files(root)

  writeLines("old but changed", stale)

  expect_identical(
    basename(SPARKS:::.changed_result_files(snapshot, root)),
    basename(stale)
  )
  unlink(root, recursive = TRUE)
})
