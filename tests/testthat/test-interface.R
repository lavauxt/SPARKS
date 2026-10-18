test_that("result catalog classifies files and excludes targets internals", {
  results <- tempfile("sparks-results-")
  dir.create(results)
  dir.create(file.path(results, "_targets"))
  writeLines("a,b\n1,2", file.path(results, "table.csv"))
  writeLines("private", file.path(results, "_targets", "private.rds"))
  writeBin(as.raw(c(1, 2, 3)), file.path(results, "plot.png"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  catalog <- SPARKS::sparks_results(results)
  expect_setequal(catalog$name, c("table.csv", "plot.png"))
  expect_identical(catalog$type[match("table.csv", catalog$name)], "Table")
  expect_identical(catalog$type[match("plot.png", catalog$name)], "Plot")
  expect_true(all(file.exists(catalog$path)))
})

test_that("report catalog lists only HTML reports", {
  results <- tempfile("sparks-report-catalog-")
  dir.create(file.path(results, "group"), recursive = TRUE)
  writeLines("<html></html>", file.path(results, "group", "QC_report.html"))
  writeLines("<html></html>", file.path(results, "group", "Results_report.html"))
  writeLines("table", file.path(results, "group", "counts.tsv"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  reports <- SPARKS:::.report_catalog(results)
  expect_equal(nrow(reports), 2L)
  expect_true(all(reports$type == "Report"))
  expect_setequal(basename(reports$path),
                  c("QC_report.html", "Results_report.html"))
})

test_that("Reports tab refreshes and selects available HTML reports", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  skip_if_not_installed("base64enc")
  skip_if_not_installed("htmltools")

  app <- SPARKS::sparks_app()
  results <- tempfile("sparks-report-browser-")
  dir.create(results)
  on.exit(unlink(results, recursive = TRUE), add = TRUE)
  writeLines("<html><body>report</body></html>",
             file.path(results, "Results_report_group.html"))

  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(report_dir = results, report_refresh = 1,
                      selected_report = file.path(results,
                                                  "Results_report_group.html"))
    session$flushReact()
    expect_match(output$report_status, "1 report")
    expect_identical(
      selected_report(),
      normalizePath(file.path(results, "Results_report_group.html"),
                    winslash = "/")
    )
  })
})

test_that("Shiny navigation uses task-focused top-level tabs", {
  expect_identical(
    SPARKS:::.sparks_nav_labels,
    c(run = "Run", explore = "Explore", qc = "QC & annotation",
      results = "Results", reports = "Reports")
  )
})

test_that("QC and results templates provide summaries and report navigation", {
  qc <- paste(readLines(testthat::test_path("..", "..", "inst", "rmd",
                                            "qc_report.Rmd"), warn = FALSE),
              collapse = "\n")
  results <- paste(readLines(testthat::test_path("..", "..", "inst", "rmd",
                                                 "results_report.Rmd"),
                              warn = FALSE), collapse = "\n")

  expect_match(qc, "Report guide")
  expect_match(qc, "QC metrics available")
  expect_match(results, "Report guide")
  expect_match(results, "SingleR annotations")
  expect_lt(regexpr("SPARKS_DOWNLOADS", results, fixed = TRUE)[[1L]],
            regexpr("Session Info", results, fixed = TRUE)[[1L]])
  expect_gt(regexpr("Pipeline diagnostics", results, fixed = TRUE)[[1L]],
            regexpr("Session Info", results, fixed = TRUE)[[1L]])
  expect_lt(regexpr("SPARKS_DOWNLOADS", qc, fixed = TRUE)[[1L]],
            regexpr("Session info", qc, fixed = TRUE)[[1L]])
  expect_gt(regexpr("Pipeline diagnostics", qc, fixed = TRUE)[[1L]],
            regexpr("Session info", qc, fixed = TRUE)[[1L]])
})

test_that("QC and results reports render with summary and section order", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("Seurat")

  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      c(1, 2, 3, 4, 5, 6), nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:3))
    )
  ))
  object$orig.ident <- c("s1", "s1", "s2")
  object$condition <- c("A", "A", "B")
  object$seurat_clusters <- c("0", "1", "1")
  object$singleR_labels_fine <- c("T cell", "B cell", "B cell")
  object$nFeature_RNA <- c(2, 2, 2)
  object$nCount_RNA <- c(3, 7, 11)
  object$percent.mt <- c(1, 2, 3)
  object[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(c(0, 0, 1, 1, 2, 2), ncol = 2,
                        dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_", assay = "RNA"
  )
  output_dir <- tempfile("sparks-report-smoke-")
  dir.create(output_dir)
  marker_dir <- file.path(output_dir, "fixture", "cluster", "DEG")
  dir.create(marker_dir, recursive = TRUE)
  utils::write.table(data.frame(
    gene = c("g1", "g2", "g1"),
    cluster = c("0", "0", "1"),
    avg_log2FC = c(1.2, -0.4, 0.7),
    p_val_adj = c(0.01, 0.2, 0.03)
  ), file.path(marker_dir, "AllMarkers_fixture_by_seurat_clusters.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE)
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)
  qc_template <- testthat::test_path("..", "..", "inst", "rmd", "qc_report.Rmd")
  results_template <- testthat::test_path("..", "..", "inst", "rmd",
                                          "results_report.Rmd")

  qc_path <- SPARKS::generate_qc_report(
    object, "fixture", output_dir, rmd_template = qc_template
  )
  results_path <- SPARKS::generate_results_report(
    object, "fixture", output_dir, rmd_template = results_template
  )
  qc_html <- paste(readLines(qc_path, warn = FALSE), collapse = "\n")
  results_html <- paste(readLines(results_path, warn = FALSE), collapse = "\n")

  expect_match(qc_html, "Report guide")
  expect_match(qc_html, "Cells after filtering")
  expect_match(qc_html, "Per-sample QC summary")
  expect_match(qc_html, "Recorded pipeline stage counts")
  expect_match(results_html, "Report guide")
  expect_match(results_html, "SingleR annotations")
  expect_match(results_html, "Marker results at a glance")
  expect_lt(regexpr("sparks-downloads", results_html, fixed = TRUE)[[1L]],
            regexpr("Session Info", results_html, fixed = TRUE)[[1L]])
  expect_gt(regexpr("Pipeline diagnostics", results_html, fixed = TRUE)[[1L]],
            regexpr("Session Info", results_html, fixed = TRUE)[[1L]])
})

test_that("targets helper group calculations agree with independent mode", {
  table <- data.frame(folder_id = c("a", "b"),
                      comparison_group = c("g", "g"))
  cfg <- list(pipeline = list(independent_samples = TRUE))
  result <- SPARKS:::.targets_sample_table(cfg, table)
  expect_identical(result$comparison_group, c("a", "b"))
})

test_that("run monitor formats elapsed time and target progress", {
  started <- as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
  expect_identical(SPARKS:::.format_run_elapsed(started,
                                                started + 3661), "01:01:01")
  expect_identical(SPARKS:::.format_run_elapsed(NULL), "Not started")
  progress <- data.frame(skipped = 1L, dispatched = 2L, completed = 3L,
                         errored = 0L, canceled = 0L)
  expect_match(SPARKS:::.format_target_progress(progress), "completed: 3")
  expect_match(SPARKS:::.format_target_progress(progress), "dispatched: 2")
  expect_identical(SPARKS:::.targets_run_progress(tempfile()),
                   "Targets have not initialized yet.")
  expect_identical(SPARKS:::.format_run_elapsed(0, 3661), "01:01:01")
})

test_that("run log appends timestamped status events", {
  log <- tempfile()
  on.exit(unlink(log), add = TRUE)
  SPARKS:::.log_run_event(log, "STARTING", "unit test")
  contents <- readLines(log)
  expect_length(contents, 1L)
  expect_match(contents, "STARTING unit test")
  expect_match(contents, "[0-9]{4}-[0-9]{2}-[0-9]{2}")
})

test_that("UI configuration validation confirms samples and target graph", {
  skip_if_not_installed("targets")
  fixture <- tempfile("ui-validate-")
  dir.create(fixture)
  base <- file.path(fixture, "base.yaml")
  override <- file.path(fixture, "override.yaml")
  template <- testthat::test_path("..", "..", "inst",
                                  "config_template_mouse.yaml")
  writeLines(readLines(template), base)
  writeLines(c("pipeline:", "  sample_table:", "    - folder_id: s1",
               "      protocol: control", "      comparison_group: g1"),
             override)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  validated <- SPARKS:::.validate_ui_configuration(
    base, override, file.path(fixture, "results"), "sample_00001"
  )
  expect_equal(nrow(validated$sample_table), 1L)
  expect_true("group_g1" %in% validated$targets)
  expect_error(SPARKS:::.validate_ui_configuration(
    base, override, "", NULL
  ), "non-empty")
})

test_that("UI override merges common controls without duplicate YAML keys", {
  fixture <- tempfile("ui-override-")
  dir.create(fixture)
  override <- file.path(fixture, "override.yaml")
  writeLines(c("pipeline:", "  species_target: Mouse", "qc:",
               "  min_features: 100", "  max_mt_percent: 20",
               "processing:", "  cluster_resolution: 0.8",
               "deg:", "  run: false"), override)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  input <- list(yaml_text = "", override_config = override,
                species = "Human", min_features = 250, max_mt = 5,
                resolution = 0.3, run_deg = TRUE,
                run_pseudobulk = TRUE, sample_file = NULL)

  path <- SPARKS:::.write_ui_targets_override(input, fixture)
  result <- yaml::yaml.load_file(path)
  expect_identical(result$pipeline$species_target, "Human")
  expect_equal(result$qc$min_features, 250)
  expect_true(result$deg$run)
})

test_that("run action validates config before starting targets", {
  skip_if_not_installed("shiny")
  app <- SPARKS::sparks_app()
  fixture <- tempfile("shiny-run-validation-")
  dir.create(fixture)
  base <- file.path(fixture, "base.yaml")
  override <- file.path(fixture, "override.yaml")
  template <- testthat::test_path("..", "..", "inst",
                                  "config_template_mouse.yaml")
  writeLines(readLines(template), base)
  writeLines(c("pipeline:", "  sample_table:", "    - folder_id: s1",
               "      protocol: control", "      comparison_group: g1"),
             override)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  shiny::testServer(app$serverFuncSource(), {
    session$setInputs(
      base_config = base, override_config = override,
      results_dir = file.path(fixture, "results"), yaml_text = "",
      sample_file = NULL, species = "Mouse", min_features = 200,
      max_mt = 10, resolution = 0.5, run_deg = TRUE,
      run_pseudobulk = FALSE, target_names = "not_a_target", run = 1
    )
    session$flushReact()
    expect_match(output$status, "Target selection is invalid")
  })
})

test_that("results overview summarizes group, type, and size", {
  catalog <- data.frame(name = c("group_a/QC/qc.png", "group_a/counts.tsv",
                                 "group_b/report.html"),
                        path = c("a", "b", "c"),
                        type = c("Plot", "Table", "Report"),
                        size = c(1024, 2048, 3072))
  summary <- SPARKS:::.summarise_results(catalog)
  expect_setequal(summary$group, c("group_a", "group_b"))
  expect_setequal(summary$type, c("Plot", "Table", "Report"))
  expect_equal(SPARKS:::.result_overview_metrics(catalog)$files, 3L)
  expect_equal(SPARKS:::.result_overview_metrics(catalog)$groups, 2L)
})

test_that("merged object discovery only includes RData group objects", {
  results <- tempfile("merged-object-index-")
  group <- file.path(results, "g1", "RData")
  dir.create(group, recursive = TRUE)
  saveRDS(list(), file.path(group, "Merged_g1.rds"))
  saveRDS(list(), file.path(group, "other.rds"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  objects <- SPARKS:::.discover_merged_objects(results)
  expect_identical(objects$group, "g1")
  expect_identical(basename(objects$path), "Merged_g1.rds")
})

test_that("UMAP selection maps lasso keys back to cell rows", {
  cells <- data.frame(cell = c("c1", "c2", "c3"), x = 1:3)
  selected <- data.frame(key = c("c1", "c3"))
  expect_identical(SPARKS:::.selected_cell_rows(cells, selected), cells[c(1, 3), ])
  expect_equal(nrow(SPARKS:::.selected_cell_rows(NULL, selected)), 0L)
})

test_that("overview explorer handles metadata and continuous coloring", {
  expect_identical(SPARKS:::.explorer_feature_kind("gene:Cd3e", "celltype"),
                   "gene")
  expect_identical(SPARKS:::.explorer_feature_kind("celltype", "celltype"),
                   "metadata")
  expect_equal(SPARKS:::.overview_color_values(c(1, 2)), c(1, 2))
  expect_s3_class(SPARKS:::.overview_color_values(c("a", "b")), "factor")
})

test_that("explorer supports chosen reductions and gene co-expression", {
  skip_if_not_installed("Seurat")
  object <- suppressWarnings(Seurat::CreateSeuratObject(
    counts = matrix(
      c(1, 2, 3, 4, 2, 4, 6, 8), nrow = 2,
      dimnames = list(c("g1", "g2"), paste0("cell", 1:4))
    )
  ))
  object$seurat_clusters <- c("0", "0", "1", "1")
  object[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(1:8, ncol = 2,
                        dimnames = list(colnames(object), c("UMAP_1", "UMAP_2"))),
    key = "UMAP_", assay = "RNA"
  )
  object[["pca"]] <- Seurat::CreateDimReducObject(
    embeddings = matrix(8:1, ncol = 2,
                        dimnames = list(colnames(object), c("PC_1", "PC_2"))),
    key = "PC_", assay = "RNA"
  )

  pca_data <- SPARKS:::.explorer_plot_data(object, "seurat_clusters",
                                            reduction = "pca")
  coexpression <- SPARKS:::.gene_coexpression_frame(object, "g1", "g2")
  expect_equal(as.numeric(pca_data$x),
               as.numeric(Seurat::Embeddings(object, "pca")[, 1L]))
  expect_equal(nrow(coexpression), 4L)
  expect_true(all(coexpression$expression_y > coexpression$expression_x))
  expect_error(SPARKS:::.gene_coexpression_frame(object, "g1", "g1"),
               "two different genes")
})

test_that("overview composition summaries count metadata levels", {
  counts <- data.frame(group = c("A", "A", "B"),
                       value = c("T", "T", "B"))
  summary <- dplyr::count(counts, group, value, name = "cells") |>
    dplyr::rename(level = value)
  expect_equal(sum(summary$cells), 3L)
  expect_setequal(summary$level, c("T", "B"))
})

test_that("result catalog rejects missing result directory", {
  expect_error(SPARKS::sparks_results(file.path(tempdir(), "missing-results")),
               "does not exist")
})

test_that("consolidated report embeds tabs, previews, and downloadable files", {
  skip_if_not_installed("htmltools")
  skip_if_not_installed("base64enc")
  results <- tempfile("sparks-results-")
  dir.create(results)
  writeLines("group,value\na,1\nb,2", file.path(results, "counts.tsv"))
  writeBin(as.raw(c(137, 80, 78, 71)), file.path(results, "plot.png"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  report <- SPARKS::generate_all_results_report(results)
  expect_true(file.exists(report))
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(html, "data-tab")
  expect_match(html, "counts.tsv")
  expect_match(html, "data:image/png;base64")
  expect_match(html, "data:text/csv;base64")
  expect_match(html, "<th>group</th>")
  expect_match(html, "<th>value</th>")
  expect_match(html, "Download original")
  expect_match(html, "Download preview as CSV")
  expect_match(html, "<!DOCTYPE html>")
  expect_equal(length(gregexpr("<html", html, fixed = TRUE)[[1]]), 1L)
  expect_equal(length(gregexpr("</html>", html, fixed = TRUE)[[1]]), 1L)
  expect_match(html, "if(tab)tab.classList.add('active')", fixed = TRUE)
  expect_false("SPARKS_all_results.html" %in% SPARKS::sparks_results(results)$name)
})

test_that("consolidated report renders an empty result directory", {
  skip_if_not_installed("htmltools")
  skip_if_not_installed("base64enc")
  results <- tempfile("empty-sparks-results-")
  dir.create(results)
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  report <- SPARKS::generate_all_results_report(results)
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(html, "No result files were found")
  expect_match(html, "SPARKS — all results")
  expect_match(html, "<!DOCTYPE html>")
  expect_false(grepl("<button data-tab", html, fixed = TRUE))
})

test_that("consolidated report rejects a directory output path", {
  skip_if_not_installed("htmltools")
  skip_if_not_installed("base64enc")
  results <- tempfile("report-output-")
  dir.create(results)
  out_dir <- file.path(results, "destination")
  dir.create(out_dir)
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  expect_error(SPARKS::generate_all_results_report(results, out_dir),
               "points to a directory")
})

test_that("sample table source resolves paths relative to override config", {
  config_dir <- tempfile("sparks-config-")
  dir.create(config_dir)
  sample_path <- file.path(config_dir, "samples.tsv")
  writeLines("folder_id\tprotocol\tcomparison_group\na\tx\tg", sample_path)
  on.exit(unlink(config_dir, recursive = TRUE), add = TRUE)
  cfg <- list(pipeline = list(sample_table = "samples.tsv",
                              config_dir = config_dir))

  expect_identical(SPARKS:::.sample_table_source(cfg, "base.yaml"),
                   normalizePath(sample_path, winslash = "/"))
})

test_that("sample target returns an explicit skip for unavailable data", {
  cfg_path <- testthat::test_path("..", "..", "inst",
                                  "config_template_mouse.yaml")
  cfg <- SPARKS::load_pipeline_config(cfg_path)
  cfg$pipeline$data_dir <- tempfile("no-such-input-")
  sample_table <- data.frame(folder_id = "missing", protocol = "control",
                             comparison_group = "group")
  output_dir <- tempfile("sparks-target-")
  on.exit(unlink(output_dir, recursive = TRUE), add = TRUE)

  sample_dir <- file.path(cfg$pipeline$data_dir, "missing")
  dir.create(sample_dir, recursive = TRUE)
  writeLines("fixture", file.path(sample_dir, "raw.txt"))
  on.exit(unlink(cfg$pipeline$data_dir, recursive = TRUE), add = TRUE)
  result <- SPARKS:::.process_sample_target(
    cfg, sample_table, 1L, output_dir,
    SPARKS:::.sample_input_files(SPARKS:::.sample_stage_config(cfg),
                                 sample_table, 1L)
  )
  expect_true(file.exists(result))
  expect_null(readRDS(result))
})

test_that("sample file targets declare generated QC and matrix sidecars", {
  results <- tempfile("sample-output-files-")
  rds <- file.path(results, "_targets", "samples", "sample.rds")
  qc <- file.path(results, "Sample_QC", "sample", "qc.png")
  matrix <- file.path(results, "RawMatrix", "sample", "counts.tsv")
  dir.create(dirname(rds), recursive = TRUE)
  dir.create(dirname(qc), recursive = TRUE)
  dir.create(dirname(matrix), recursive = TRUE)
  writeBin(as.raw(1), rds)
  writeBin(as.raw(2), qc)
  writeLines("gene\tcount", matrix)
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  outputs <- SPARKS:::.sample_target_outputs(rds, results, "sample")
  expect_setequal(basename(outputs), c("sample.rds", "qc.png", "counts.tsv"))
})

test_that("results manifest records catalogued output files", {
  results <- tempfile("sparks-results-")
  dir.create(results)
  writeLines("x\n1", file.path(results, "table.tsv"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  report_path <- file.path(results, "report.html")
  writeLines("report", report_path)
  manifest <- SPARKS:::.write_results_manifest(results, report_path)
  expect_true(file.exists(manifest))
  expect_match(paste(readLines(manifest), collapse = "\n"), "table.tsv")
})

test_that("targets runner reports missing optional dependency or invalid config", {
  if (!requireNamespace("targets", quietly = TRUE)) {
    expect_error(SPARKS::run_sparks_targets("missing.yaml"), "requires package 'targets'")
  } else {
    expect_error(SPARKS::run_sparks_targets("missing.yaml"), "Base config not found")
  }
})

test_that("target selection validates explicit target names", {
  expect_identical(
    SPARKS:::.validate_targets_selection(c("sample_00001", "group_g1"),
                                         c("sample_00001", "group_g1")),
    invisible(c("sample_00001", "group_g1"))
  )
  expect_error(SPARKS:::.validate_targets_selection("unknown", "sample_00001"),
               "unique names")
  expect_error(SPARKS:::.validate_targets_selection(c("x", "x"), "x"),
               "unique names")
  expect_error(SPARKS:::.validate_targets_selection(character(), "x"),
               "unique names")
})

test_that("tabular results can be downloaded as CSV or TSV", {
  source <- tempfile(fileext = ".tsv")
  csv <- tempfile(fileext = ".csv")
  tsv <- tempfile(fileext = ".tsv")
  writeLines("sample\tcount\ns1\t2", source)
  on.exit(unlink(c(source, csv, tsv)), add = TRUE)

  SPARKS:::.write_result_download(source, csv, "CSV")
  SPARKS:::.write_result_download(source, tsv, "TSV")
  expect_match(paste(readLines(csv), collapse = "\n"), "sample.*count")
  expect_match(paste(readLines(tsv), collapse = "\n"), "sample\\tcount")
})

test_that("result download helper preserves original non-table files", {
  source <- tempfile(fileext = ".png")
  dest <- tempfile(fileext = ".png")
  writeBin(as.raw(c(1, 2, 3)), source)
  on.exit(unlink(c(source, dest)), add = TRUE)

  SPARKS:::.write_result_download(source, dest)
  expect_identical(readBin(dest, "raw", n = 3L), as.raw(c(1, 2, 3)))
})

test_that("targets runner releases its lock after target selection errors", {
  skip_if_not_installed("targets")
  fixture <- tempfile("targets-lock-")
  dir.create(fixture)
  config <- file.path(fixture, "config.yaml")
  override <- file.path(fixture, "override.yaml")
  template <- testthat::test_path("..", "..", "inst",
                                  "config_template_mouse.yaml")
  writeLines(readLines(template), config)
  writeLines(c("pipeline:", "  sample_table:", "    - folder_id: s1",
               "      protocol: control", "      comparison_group: g1"),
             override)
  out <- file.path(fixture, "results")
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  expect_error(SPARKS::run_sparks_targets(config, override, results_dir = out,
                                          target_names = "not_a_target"),
               "unique names")
  expect_false(dir.exists(file.path(out, "_targets", "sparks-run.lock")))
})

test_that("targets runner refuses to share an active store", {
  skip_if_not_installed("targets")
  fixture <- tempfile("targets-lock-contention-")
  dir.create(fixture)
  config <- file.path(fixture, "config.yaml")
  template <- testthat::test_path("..", "..", "inst",
                                  "config_template_mouse.yaml")
  writeLines(readLines(template), config)
  out <- file.path(fixture, "results")
  lock <- file.path(out, "_targets", "sparks-run.lock")
  dir.create(lock, recursive = TRUE)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  expect_error(SPARKS::run_sparks_targets(config, results_dir = out),
               "already holds this store lock")
  expect_true(dir.exists(lock))
})

test_that("sample stage config only includes sample-preprocessing settings", {
  cfg <- list(
    pipeline = list(data_dir = "raw", results_dir = "results"),
    input = list(format = "10x"),
    species = list(gene_removal_pattern = "^mt-", mt_pattern = "^mt-",
                   genes_to_remove = "gene", ref_primary = "reference"),
    qc = list(min_features = 200),
    processing = list(cluster_resolution = 0.7),
    deg = list(run = FALSE)
  )
  stage <- SPARKS:::.sample_stage_config(cfg)
  expect_named(stage, c("pipeline", "input", "species", "qc", "plot",
                        "parallel"))
  expect_identical(stage$pipeline$results_dir, "results")
  expect_null(stage$species$ref_primary)
})

test_that("analysis stage config isolates shared config and run output path", {
  cfg <- list(pipeline = list(results_dir = "old"), qc = list(min_features = 1))
  analysis <- SPARKS:::.analysis_stage_config(cfg, "new")
  expect_s3_class(analysis, "sparks_analysis_config")
  expect_identical(analysis$cfg, cfg)
  expect_identical(analysis$results_dir, "new")
})

test_that("independent sample grouping follows pipeline semantics", {
  table <- data.frame(folder_id = c("a", "b"),
                      comparison_group = c("g", "g"))
  cfg <- list(pipeline = list(independent_samples = TRUE))
  transformed <- SPARKS:::.targets_sample_table(cfg, table)
  expect_identical(transformed$comparison_group, c("a", "b"))
})

test_that("cached sample objects are scoped to their comparison group", {
  cached <- list(1L, 2L, 3L)
  expect_identical(SPARKS:::.load_cached_samples(cached, c(1L, 3L)),
                   list(1L, 3L))
})

test_that("targets config prevents independent-sample regrouping twice", {
  cfg <- list(pipeline = list(results_dir = "old",
                              independent_samples = TRUE))
  grouped <- SPARKS:::.group_analysis_config(cfg, "new")
  expect_identical(grouped$pipeline$results_dir, "new")
  expect_false(grouped$pipeline$independent_samples)
})

test_that("group row selection matches comparison group", {
  table <- data.frame(comparison_group = c("a", "b", "a"))
  expect_identical(SPARKS:::.targets_sample_rows(table, "a"), c(1L, 3L))
})

test_that("targets plan separates sample, group, report, and manifest stages", {
  skip_if_not_installed("targets")
  fixture <- tempfile("targets-config-")
  dir.create(fixture)
  base <- file.path(fixture, "base.yaml")
  override <- file.path(fixture, "override.yaml")
  writeLines(readLines(testthat::test_path("..", "..", "inst",
                                           "config_template_mouse.yaml")), base)
  writeLines(c("pipeline:", "  sample_table:", "    - folder_id: s1",
               "      protocol: control", "      comparison_group: g1",
               "    - folder_id: s2", "      protocol: treated",
               "      comparison_group: g1"), override)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  plan <- SPARKS:::.build_sparks_targets_plan(base, override,
                                              file.path(fixture, "results"))
  names <- vapply(plan, `[[`, character(1), "name")
  expect_true(all(c("sample_00001", "sample_00002", "group_g1",
                    "consolidated_report", "results_manifest") %in% names))
  expect_true(all(c("sample_input_00001", "sample_input_00002") %in% names))
  expect_true(all(c("analysis_config", "reports_g1") %in% names))
  script <- file.path(fixture, "_targets.R")
  targets::tar_script(plan, script = script, ask = FALSE)
  expect_no_error(targets::tar_validate(callr_function = NULL, script = script,
                                       store = file.path(fixture, "store")))
})

test_that("targets plan rejects report target name collisions", {
  skip_if_not_installed("targets")
  fixture <- tempfile("targets-collision-")
  dir.create(fixture)
  base <- file.path(fixture, "base.yaml")
  override <- file.path(fixture, "override.yaml")
  writeLines(readLines(testthat::test_path("..", "..", "inst",
                                           "config_template_mouse.yaml")), base)
  writeLines(c("pipeline:", "  sample_table:", "    - folder_id: s1",
               "      protocol: a", "      comparison_group: g-1",
               "    - folder_id: s2", "      protocol: b",
               "      comparison_group: g_1"), override)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)

  expect_error(SPARKS:::.build_sparks_targets_plan(base, override,
                                                   file.path(fixture, "results")),
               "collide")
})

test_that("targets dependency failure is actionable", {
  if (requireNamespace("targets", quietly = TRUE)) skip("targets is installed")
  cfg <- testthat::test_path("..", "..", "inst",
                             "config_template_mouse.yaml")
  expect_error(SPARKS:::.build_sparks_targets_plan(cfg, NULL, tempdir()),
               "Package 'targets' is required")
})

test_that("sample targets are row-specific and group target settings are isolated", {
  fixture <- tempfile("sample-target-")
  raw <- file.path(fixture, "s1")
  dir.create(raw, recursive = TRUE)
  writeLines("test", file.path(raw, "matrix.mtx"))
  table <- data.frame(folder_id = "s1", protocol = "control",
                      comparison_group = "g1")
  sample_cfg <- list(pipeline = list(data_dir = fixture,
                                    results_dir = fixture),
                     input = list(format = "auto", save_raw_matrix = FALSE),
                     species = list(gene_removal_pattern = "", mt_pattern = "^mt-",
                                    genes_to_remove = character()),
                     qc = list(min_features = 1, max_features = 100,
                               max_counts = 1000, max_mt_percent = 100,
                               min_cells = 1),
                     plot = list(), parallel = list(enable = FALSE))
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  sample_file <- SPARKS:::.process_sample_target(
    sample_cfg, table, 1L, fixture,
    SPARKS:::.sample_input_files(sample_cfg, table, 1L)
  )
  expect_true(file.exists(sample_file))
  expect_named(SPARKS:::.targets_sample_rows(
    data.frame(comparison_group = c("g1", "g2")), "g1"), NULL)
})

test_that("target sample path helper filters failed sample branches", {
  existing <- tempfile(fileext = ".rds")
  saveRDS(NULL, existing)
  on.exit(unlink(existing), add = TRUE)

  expect_identical(SPARKS:::.target_sample_paths(c(existing, NA_character_)),
                   existing)
})

test_that("target sample path helper retains one null sample branch as a path", {
  existing <- tempfile(fileext = ".rds")
  saveRDS(NULL, existing)
  on.exit(unlink(existing), add = TRUE)

  expect_identical(SPARKS:::.target_sample_paths(list(existing)), existing)
})

test_that("group artifact target only includes its own outputs", {
  results <- tempfile("sparks-results-")
  dir.create(file.path(results, "group_a"), recursive = TRUE)
  dir.create(file.path(results, "group_b"), recursive = TRUE)
  writeLines("a", file.path(results, "group_a", "a.tsv"))
  writeLines("b", file.path(results, "group_b", "b.tsv"))
  on.exit(unlink(results, recursive = TRUE), add = TRUE)

  paths <- SPARKS:::.group_result_files(results, "group_a")
  expect_identical(basename(paths), "a.tsv")
})

test_that("sample file inventory is tracked for each sample input", {
  root <- tempfile("sample-inputs-")
  sample_dir <- file.path(root, "s1")
  dir.create(sample_dir, recursive = TRUE)
  writeLines("matrix", file.path(sample_dir, "matrix.mtx"))
  writeLines("genes", file.path(sample_dir, "features.tsv"))
  cfg <- list(pipeline = list(data_dir = root))
  table <- data.frame(folder_id = "s1")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  paths <- SPARKS:::.sample_input_files(cfg, table, 1L)
  expect_setequal(basename(paths), c("matrix.mtx", "features.tsv"))
})

test_that("Shiny app launcher returns an app object when optional UI packages exist", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")

  app <- SPARKS::sparks_app()
  expect_s3_class(app, "shiny.appobj")
  response <- app$httpHandler(list(
    REQUEST_METHOD = "GET", PATH_INFO = "/", QUERY_STRING = "",
    SCRIPT_NAME = ""
  ))
  expect_s3_class(response, "httpResponse")
  expect_equal(response$status, 200L)
})
