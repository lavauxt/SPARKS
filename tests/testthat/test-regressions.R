test_that("missing override config paths fail instead of falling back", {
  config <- tempfile(fileext = ".yml")
  on.exit(unlink(config), add = TRUE)
  yaml::write_yaml(list(), config)

  expect_error(
    SPARKS::load_pipeline_config(config, tempfile("missing-override-")),
    "Override config not found"
  )
})

test_that("sample metadata is validated and normalized regardless of source", {
  expect_error(
    SPARKS:::.validate_sample_table(data.frame(folder_id = "s1")),
    "missing columns: protocol, comparison_group"
  )
  expect_error(
    SPARKS:::.validate_sample_table(data.frame()),
    "empty"
  )

  normalized <- SPARKS:::.validate_sample_table(data.frame(
    folder_id = 1,
    protocol = 2,
    comparison_group = 3,
    quant_format = NA
  ))
  expect_type(normalized$folder_id, "character")
  expect_type(normalized$protocol, "character")
  expect_type(normalized$comparison_group, "character")
  expect_type(normalized$quant_format, "character")
})

test_that("cell-level DEG skips designs with more than two conditions", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(1:6, nrow = 2, dimnames = list(c("g1", "g2"),
                                                   paste0("cell", 1:3)))
  ))
  object$condition <- c("A", "B", "C")
  object$group <- "same"

  expect_message(
    result <- SPARKS::run_deg_analysis(
      object, logfc_threshold = 0.25, min_pct = 0.1, group_by_col = "group"
    ),
    "exactly 2 conditions are required"
  )
  expect_null(result)
})

test_that("condition violin plots include cell points", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:16,
      nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:8))
    )
  ))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object$condition <- rep(c("A", "B"), each = 4)
  object$group <- rep(c("x", "y"), each = 4)

  plot <- SPARKS:::.make_condition_violin(
    object, feature = "g1", group_by_col = "group",
    split_by = "condition", assay = "RNA"
  )

  expect_true(any(vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomPoint")
  }, logical(1))))
})

test_that("split condition DotPlots use a discrete legend", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:24,
      nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:12))
    )
  ))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object$condition <- rep(c("A", "B"), each = 6)
  object$group <- rep(c("x", "y", "z"), each = 4)
  out_dir <- tempfile("signature-plots-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_no_warning(SPARKS::generate_gene_signature_plots(
    object, genes = c("g1", "g2"), out_dir = out_dir, prefix = "test",
    group_by_col = "group", signature_name = "test", condition_col = "condition"
  ))
})

test_that("small proportion tables use a simulated chi-squared p-value", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:8,
      nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:4))
    )
  ))
  object$condition <- c("A", "A", "B", "B")
  object$group <- c("x", "x", "x", "y")
  out_dir <- tempfile("proportion-test-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  expect_no_warning(SPARKS::run_proportion_analysis(
    object, grouping_col = "group", output_dir = out_dir, file_prefix = "test"
  ))
  expect_true(file.exists(file.path(out_dir, "ChiSq_test_group.txt")))
})
