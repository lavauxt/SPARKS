test_that("the report's JavaScript parses", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "node is not available")
  js <- system.file("report", "sparks_report.js", package = "SPARKS")
  skip_if(!nzchar(js), "report assets not found")
  # `node --check` is a syntax check only: it never runs the script.
  out <- suppressWarnings(system2(node, c("--check", shQuote(js)), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
})

test_that("the report's CSS and JS ship with the package", {
  for (f in c("sparks_report.js", "sparks_report.css")) {
    expect_true(nzchar(system.file("report", f, package = "SPARKS")), info = f)
  }
})

test_that("the JS only references data fields the R side writes", {
  js <- paste(readLines(system.file("report", "sparks_report.js", package = "SPARKS"), warn = FALSE), collapse = "\n")
  # top-level DATA.* members used by the client
  used <- unique(sub("^DATA\\.", "", regmatches(js, gregexpr("DATA\\.[a-z_]+", js))[[1L]]))
  written <- c("meta", "theme_map", "kind_labels", "summary", "datasets", "tables", "figures",
               "images", "files", "run")
  expect_true(all(used %in% written), info = paste(setdiff(used, written), collapse = ", "))
})
