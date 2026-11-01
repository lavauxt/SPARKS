# Extract the embedded JSON data block from a built report.
report_data <- function(html_path) {
  html <- paste(readLines(html_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  m <- regmatches(html, regexec('<script id="sparks-data" type="application/json">(.*?)</script>', html, perl = TRUE))[[1L]]
  jsonlite::fromJSON(gsub("<\\/", "</", m[[2L]], fixed = TRUE), simplifyVector = FALSE)
}

build_demo_report <- function() {
  skip_if_no_report_stack()
  dir <- demo_copy(parent.frame())
  out <- file.path(dir, "SPARKS_all_results.html")
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE))
  list(dir = dir, path = out, html = paste(readLines(out, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
}

test_that("sparks_report writes one self-contained HTML file", {
  r <- build_demo_report()
  expect_true(file.exists(r$path))
  expect_gt(file.size(r$path), 1e6)   # plotly.js is inlined
  # no external scripts, stylesheets or fonts: it must work offline
  expect_false(grepl("<script[^>]+src=", r$html))
  expect_false(grepl("<link[^>]+href=", r$html))
  expect_false(grepl("https?://[^\"' ]+\\.(js|css|woff2?)", r$html))
  # plotly.js ships exactly once (its bundle sets window.PLOTLYENV)
  expect_equal(lengths(regmatches(r$html, gregexpr("window.PLOTLYENV", r$html, fixed = TRUE))), 1L)
  expect_equal(lengths(regmatches(r$html, gregexpr("<script", r$html, fixed = TRUE))), 3L)
  expect_match(r$html, "<title>SPARKS report - ", fixed = FALSE)
  expect_match(r$html, "Synthetic demo data", fixed = TRUE)
})

test_that("the report data block is valid JSON with every section", {
  r <- build_demo_report()
  d <- report_data(r$path)
  expect_named(d, c("meta", "theme_map", "kind_labels", "summary", "datasets", "tables",
                    "figures", "images", "files", "run"), ignore.order = TRUE)
  expect_true(d$meta$demo)
  expect_equal(d$summary$cells, 240)
  expect_equal(d$summary$samples, 4)
  expect_gte(length(d$datasets), 1L)
  ds <- d$datasets[[1L]]
  expect_length(ds$cells, 240L)
  expect_length(ds$x, 240L)
  expect_true(all(c("cell_type", "orig.ident", "condition") %in% names(ds$meta)))
  expect_equal(ds$meta$cell_type$type, "cat")
  expect_length(ds$meta$cell_type$codes, 240L)
  expect_true("CD3E" %in% names(ds$genes))
  # sparse gene entries are arrays even when they hold one value
  expect_true(is.list(ds$genes$CD3E$i))
  kinds <- vapply(d$tables, function(t) t$kind, character(1))
  expect_true(all(c("markers_all", "gsea", "proportions", "deg_counts") %in% kinds))
  tabs <- vapply(d$figures, function(f) f$tab, character(1))
  expect_true(all(c("overview", "qc", "de", "enrichment", "composition") %in% tabs))
  # every figure spec was spliced in as an object, not left as a placeholder string
  expect_false(any(vapply(d$figures, function(f) is.character(f$spec), logical(1))))
  expect_false(grepl("@@SPARKS_FIG", r$html, fixed = TRUE))
})

test_that("tables are column-oriented, rounded and capped", {
  skip_if_no_report_stack()
  dir <- demo_copy()
  out <- file.path(dir, "r.html")
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE, max_table_rows = 10L))
  d <- report_data(out)
  mk <- Filter(function(t) t$kind == "markers_all", d$tables)[[1L]]
  expect_lte(mk$n, 10L)
  expect_true(mk$truncated)
  expect_gt(mk$n_total, 10L)
  expect_true(all(c("gene", "cluster", "avg_log2FC", "p_val_adj") %in% unlist(mk$cols)))
  expect_length(mk$data$gene, mk$n)
  # a capped marker table keeps every cluster represented
  expect_length(unique(unlist(mk$data$cluster)), 3L)
})

test_that("a re-run of the report reflects interactive analyses", {
  skip_if_no_report_stack()
  skip_if_not_installed("targets")
  dir <- demo_copy()
  b <- sparks_read_results(dir)
  obj <- readRDS(b$objects$path[[1L]])
  sparks_interactive_submit(dir, sparks_spec("de", cells_a = colnames(obj)[obj$cell_type == "B cells"], label = "B vs rest"), b$objects$path[[1L]])
  sparks_interactive_run(dir, b$objects$path[[1L]])
  out <- file.path(dir, "r.html")
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE))
  d <- report_data(out)
  ex <- Filter(function(t) t$kind == "explorer", d$tables)
  expect_length(ex, 1L)
  expect_match(ex[[1L]]$label, "B vs rest")
  # the report indexes the file tree but never its own interactive store
  expect_false(any(grepl("_targets", vapply(d$files, function(f) f$rel, character(1)))))
})

test_that("images are embedded within the budget and linked beyond it", {
  skip_if_no_report_stack()
  dir <- demo_copy()
  png_dir <- file.path(dir, "Demo", "QC")
  grDevices::png(file.path(png_dir, "VlnPlot_a.png"), width = 200, height = 200); graphics::plot(1); grDevices::dev.off()
  grDevices::png(file.path(png_dir, "VlnPlot_b.png"), width = 200, height = 200); graphics::plot(2); grDevices::dev.off()
  out <- file.path(dir, "r.html")
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE))
  imgs <- report_data(out)$images
  expect_length(imgs, 2L)
  expect_true(all(vapply(imgs, function(i) startsWith(i$src, "data:image/png;base64,"), logical(1))))
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE, embed_images = FALSE))
  imgs <- report_data(out)$images
  expect_true(all(vapply(imgs, function(i) is.null(i$src), logical(1))))
  expect_setequal(vapply(imgs, function(i) i$rel, character(1)),
                  c("Demo/QC/VlnPlot_a.png", "Demo/QC/VlnPlot_b.png"))
})

test_that("JSON stays inert inside the script element", {
  json <- .report_json(list(a = "</script><!-- x"), character())
  expect_false(grepl("</script", json, fixed = TRUE))
  expect_false(grepl("<!--", json, fixed = TRUE))
})

test_that("sparks_report needs its optional packages and says so", {
  testthat::local_mocked_bindings(.report_missing_packages = function() c("plotly", "jsonlite"))
  expect_error(sparks_report(tempdir()), "plotly, jsonlite")
})

test_that("the pipeline-level report falls back to the static page without the interactive stack", {
  skip_if_not_installed("htmltools")
  skip_if_not_installed("base64enc")
  dir <- tempfile("fallback_"); dir.create(dir)
  writeLines("a\tb\n1\t2", file.path(dir, "t.txt"))
  testthat::local_mocked_bindings(.report_missing_packages = function() "plotly")
  out <- file.path(dir, "SPARKS_all_results.html")
  res <- suppressMessages(.write_unified_report(dir, output_file = out))
  expect_true(file.exists(res))
  expect_false(grepl("sparks-data", paste(readLines(res, warn = FALSE), collapse = "")))
})

test_that("generate_all_results_report now writes the unified report", {
  skip_if_no_report_stack()
  dir <- demo_copy()
  out <- suppressMessages(generate_all_results_report(dir))
  expect_equal(basename(out), "SPARKS_all_results.html")
  expect_match(paste(readLines(out, warn = FALSE), collapse = ""), "sparks-data")
})

test_that("an empty results directory still yields a (mostly empty) report", {
  skip_if_no_report_stack()
  dir <- tempfile("empty_"); dir.create(dir)
  out <- file.path(dir, "r.html")
  suppressMessages(sparks_report(dir, output_file = out, quiet = TRUE))
  d <- report_data(out)
  expect_length(d$datasets, 0L)
  expect_length(d$tables, 0L)
})

test_that("the report: section of config_used.yaml sets report options; bad values are ignored", {
  skip_if_not_installed("yaml")
  dir <- tempfile("cfg_"); dir.create(dir)
  expect_equal(.report_options(dir), list())            # no config at all
  yaml::write_yaml(list(report = list(max_cells = 5000, max_genes = "many", embed_images = FALSE,
                                      max_table_rows = -1, image_budget_mb = 10, title = "My run",
                                      typo_key = 1)),
                   file.path(dir, "config_used.yaml"))
  opt <- .report_options(dir)
  expect_equal(opt, list(max_cells = 5000L, image_budget_mb = 10, embed_images = FALSE, title = "My run"))
  yaml::write_yaml(list(report = NULL, other = 1), file.path(dir, "config_used.yaml"))
  expect_equal(.report_options(dir), list())
  writeLines("{ not: [yaml", file.path(dir, "config_used.yaml"))
  expect_equal(.report_options(dir), list())
})

test_that("the pipeline-level report honours the configured options", {
  skip_if_no_report_stack()
  skip_if_not_installed("yaml")
  dir <- demo_copy()
  yaml::write_yaml(list(report = list(max_cells = 100, title = "Configured title")),
                   file.path(dir, "config_used.yaml"))
  out <- suppressMessages(.write_unified_report(dir))
  d <- report_data(out)
  expect_equal(d$meta$title, "Configured title")
  expect_lte(length(d$datasets[[1L]]$cells), 100L)
  expect_equal(d$datasets[[1L]]$n_total, 240)
})
