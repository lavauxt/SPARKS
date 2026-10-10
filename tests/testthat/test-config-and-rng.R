test_that("config merge replaces YAML sequences and preserves null", {
  base <- list(
    subsets = list(list(name = "template")),
    singler = list(labels = list(list(name = "main"))),
    groupings_main = "clusters"
  )
  override <- list(
    subsets = list(),
    singler = list(labels = list()),
    groupings_main = NULL
  )

  merged <- SPARKS:::.merge_config(base, override)

  expect_length(merged$subsets, 0L)
  expect_length(merged$singler$labels, 0L)
  expect_null(merged$groupings_main)
})

test_that("seed helper restores the caller RNG state", {
  set.seed(8137)
  expected_seed <- get(".Random.seed", envir = .GlobalEnv)
  restore_rng <- SPARKS:::.set_seed_preserving_state(2861L)
  seeded_draw <- runif(1)
  restore_rng()

  expect_identical(get(".Random.seed", envir = .GlobalEnv), expected_seed)
  expect_type(seeded_draw, "double")
})

test_that("seed helper produces the same stream for the same seed", {
  draw_once <- function() {
    restore_rng <- SPARKS:::.set_seed_preserving_state(73L)
    on.exit(restore_rng(), add = TRUE)
    runif(3)
  }

  expect_identical(draw_once(), draw_once())
})

test_that("seed helper restores an uninitialized RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) set.seed(8137), add = TRUE)

  restore_rng <- SPARKS:::.set_seed_preserving_state(2861L)
  runif(1)
  restore_rng()

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("seed helper rejects invalid values", {
  expect_error(SPARKS:::.validate_seed(NA_integer_), "seed must be")
  expect_error(SPARKS:::.validate_seed(-1L), "seed must be")
  expect_error(SPARKS:::.validate_seed(1.5), "seed must be")
  expect_error(SPARKS:::.validate_seed(c(1L, 2L)), "seed must be")
})

test_that("report template lookup honors override and search precedence", {
  config_dir <- tempfile("sparks-config-")
  working_dir <- tempfile("sparks-working-")
  dir.create(config_dir)
  dir.create(working_dir)
  on.exit(unlink(c(config_dir, working_dir), recursive = TRUE), add = TRUE)

  filename <- "qc_report.Rmd"
  config_template <- file.path(config_dir, filename)
  working_template <- file.path(working_dir, filename)
  writeLines("config", config_template)
  writeLines("working", working_template)

  expect_identical(
    SPARKS:::.find_report_template(NULL, config_dir, filename),
    config_template
  )
  expect_identical(
    SPARKS:::.find_report_template("custom.Rmd", config_dir, filename),
    "custom.Rmd"
  )

  unlink(config_template)
  old_wd <- setwd(working_dir)
  on.exit(setwd(old_wd), add = TRUE)
  expect_identical(
    SPARKS:::.find_report_template(NULL, config_dir, filename),
    normalizePath(working_template, winslash = "/")
  )
})

test_that("report template lookup finds package templates without system.file", {
  config_dir <- tempfile("sparks-config-")
  working_dir <- tempfile("sparks-working-")
  dir.create(config_dir)
  dir.create(working_dir)
  on.exit(unlink(c(config_dir, working_dir), recursive = TRUE), add = TRUE)

  old_wd <- setwd(working_dir)
  on.exit(setwd(old_wd), add = TRUE)

  template <- SPARKS:::.find_report_template(NULL, config_dir, "qc_report.Rmd")

  expect_true(file.exists(template))
  expect_identical(basename(template), "qc_report.Rmd")
})

test_that("pipeline log handlers preserve messages and warnings", {
  log_file <- tempfile()
  log_con <- file(log_file, open = "wt")
  on.exit({
    close(log_con)
    unlink(log_file)
  }, add = TRUE)
  logged <- character()
  handler <- function(type, condition) {
    logged <<- c(logged, paste(type, conditionMessage(condition)))
    writeLines(tail(logged, 1L), log_con)
    flush(log_con)
  }
  expect_warning(
    withCallingHandlers(
      {
        message("test message")
        warning("test warning", call. = FALSE)
      },
      message = function(m) handler("MESSAGE", m),
      warning = function(w) handler("WARNING", w)
    ),
    "test warning"
  )
  expect_true(any(grepl("MESSAGE test message", logged, fixed = TRUE)))
  expect_true(any(grepl("WARNING test warning", logged, fixed = TRUE)))
  expect_true(any(grepl("test warning", readLines(log_file), fixed = TRUE)))
})

test_that("raw matrix export rejects unsupported formats", {
  expect_error(
    SPARKS::export_raw_matrix(NULL, tempfile(), format = "csv"),
    "format.*mtx.*h5"
  )
  expect_error(
    SPARKS::export_raw_matrix(NULL, tempfile(), format = c("mtx", "h5")),
    "format.*mtx.*h5"
  )
})

test_that("only the known Seurat aggregate advisory is muffled", {
  known <- "As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis."
  other <- "another informative message"

  expect_message(
    SPARKS:::.with_seurat_aggregate_advisory_muffled(message(known)),
    NA
  )
  expect_message(
    SPARKS:::.with_seurat_aggregate_advisory_muffled(message(other)),
    "another informative message"
  )
})

test_that("QC report looks up metadata columns, not cell row names", {
  report <- readLines(testthat::test_path("..", "..", "inst", "rmd",
                                          "qc_report.Rmd"), warn = FALSE)
  expect_true(any(grepl("colnames\\(obj\\[\\[\\]\\]\\)", report)))
  expect_false(any(grepl("qc_features.*rownames\\(obj", report)))
})

test_that("results report covers current-run analyses and diagnostics", {
  report <- readLines(testthat::test_path("..", "..", "inst", "rmd",
                                          "results_report.Rmd"), warn = FALSE)
  expect_true(any(grepl("run_files", report, fixed = TRUE)))
  expect_true(any(grepl("correlation-results", report, fixed = TRUE)))
  expect_true(any(grepl("Pipeline diagnostics", report, fixed = TRUE)))
  expect_true(any(grepl("Per-cell pathway scoring", report, fixed = TRUE)))
})

test_that("sample-level and Hallmark analysis config has deliberate defaults", {
  cfg <- list(
    pipeline = list(), qc = list(), processing = list(), deg = list(),
    plot = list(), species = list(), singler = list(), groupings_main = NULL
  )
  file <- tempfile(fileext = ".yml")
  on.exit(unlink(file), add = TRUE)
  yaml::write_yaml(cfg, file)

  loaded <- SPARKS::load_pipeline_config(file)
  expect_false(loaded$pseudobulk$run)
  expect_true(loaded$pseudobulk$save_counts)
  expect_equal(loaded$pseudobulk$min_replicates, 2L)
  expect_true(loaded$gsea$run)
  expect_identical(loaded$gsea$collection, "H")
})
