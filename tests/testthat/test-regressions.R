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
