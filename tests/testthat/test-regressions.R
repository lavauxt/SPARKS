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

test_that("report file links resolve relative to the report and encode spaces", {
  root <- tempfile("report-links-")
  report_dir <- file.path(root, "group", "QC")
  result_dir <- file.path(root, "group", "tables with spaces")
  dir.create(report_dir, recursive = TRUE)
  dir.create(result_dir)
  result <- file.path(result_dir, "table one.tsv")
  writeLines("a\tb", result)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  href <- SPARKS:::.report_file_href(result, report_dir)
  expect_identical(href, "../tables%20with%20spaces/table%20one.tsv")
  expect_true(file.exists(file.path(report_dir, utils::URLdecode(href))))
  expect_error(SPARKS:::.report_file_href(file.path(root, "missing.tsv"), report_dir))
  expect_error(SPARKS:::.report_file_href(result, file.path(root, "missing-dir")))
})

test_that("report download section exposes only valid local files", {
  skip_if_not_installed("htmltools")
  root <- tempfile("report-downloads-")
  report_dir <- file.path(root, "group")
  dir.create(report_dir, recursive = TRUE)
  report <- file.path(report_dir, "report.html")
  result <- file.path(report_dir, "table.tsv")
  writeLines("<html><body><h1>Report</h1></body></html>", report)
  writeLines("a\tb", result)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  SPARKS:::.append_report_download_section(report, c(result, file.path(root, "missing.tsv")),
                                           "Downloads")
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(html, '<a href="table.tsv"')
  expect_equal(length(gregexpr("<section id=\"sparks-downloads\"", html,
                               fixed = TRUE)[[1L]]), 1L)
  expect_lt(regexpr("id=\"sparks-downloads\"", html, fixed = TRUE)[[1L]],
            regexpr("<h1>", html, fixed = TRUE)[[1L]])
  expect_match(html, "position:sticky")
  expect_match(html, "sparks-download-search")
  expect_match(html, "<details")
  expect_true(file.exists(file.path(report_dir, "table.tsv")))
  SPARKS:::.append_report_download_section(report, result, "Downloads")
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_equal(length(gregexpr("<section id=\"sparks-downloads\"", html,
                               fixed = TRUE)[[1L]]), 1L)
})

test_that("current-run file tracking includes escape results", {
  root <- tempfile("escape-report-files-")
  dir.create(root)
  files <- c("Escape_Heatmap_pathway.png", "Escape_Scores_UCell.tsv")
  for (name in files) writeLines("result", file.path(root, name))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  expect_setequal(basename(rownames(SPARKS:::.snapshot_result_files(root))), files)
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

test_that("cell-level DEG passes explicit cell identities to Seurat", {
  counts <- matrix(
    c(rep(c(1, 2, 3, 4), 6), rep(c(4, 3, 2, 1), 6)),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("cell", 1:12))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object$condition <- rep(c("A", "B"), each = 6)
  object$group <- "all"
  object[["SCT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    data = Matrix::Matrix(log1p(counts), sparse = TRUE)
  )

  result <- SPARKS::run_deg_analysis(
    object, logfc_threshold = 0, min_pct = 0, group_by_col = "group",
    min_cells_per_group = 3L
  )

  expect_true(is.data.frame(result))
  expect_true(all(result$cluster == "all"))
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

test_that("signature DotPlot dimensions stay within raster device limits", {
  expect_equal(
    SPARKS:::.signature_dotplot_dimensions(n_genes = 2, n_combo = 4),
    c(width = 14, height = 3.4)
  )
  dimensions <- SPARKS:::.signature_dotplot_dimensions(n_genes = 5, n_combo = 200)
  expect_true(all(dimensions <= 48))
})

test_that("scProportionTest skips without independent sample replication", {
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      1:8, nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:4))
    )
  ))
  object$condition <- c("A", "A", "B", "B")
  object$sample <- c("s1", "s1", "s2", "s2")
  object$group <- c("x", "x", "y", "y")

  expect_message(
    SPARKS::run_scproportion_test(
      object, "group", tempfile("scproportion-"), "test"
    ),
    "at least 2 independent samples per condition"
  )
})
