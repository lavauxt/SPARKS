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

  SPARKS:::.run_sample_level_analysis(
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
})

test_that("analysis output folders distinguish analysis types", {
  root <- tempfile("analysis-folders-")
  dirs <- SPARKS:::.make_analysis_dirs(root)

  expect_true(all(vapply(dirs, dir.exists, logical(1))))
  expect_false(identical(dirs$Correlation, dirs$DEG))
  expect_false(identical(dirs$Pseudobulk, dirs$DEG))
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
