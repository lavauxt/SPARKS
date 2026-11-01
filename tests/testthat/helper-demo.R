# Shared fixtures. Building the synthetic results takes a few seconds, so it is
# done once per test session and copied for tests that write into it.

.fixture <- new.env(parent = emptyenv())

skip_if_no_seurat <- function() {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("SeuratObject")
  testthat::skip_if_not_installed("uwot")
}

skip_if_no_report_stack <- function() {
  skip_if_no_seurat()
  for (pkg in c("plotly", "htmltools", "jsonlite", "base64enc")) {
    testthat::skip_if_not_installed(pkg)
  }
}

#' A pristine synthetic results directory (read-only by convention)
demo_results <- function() {
  skip_if_no_seurat()
  if (is.null(.fixture$dir) || !dir.exists(.fixture$dir)) {
    dir <- tempfile("sparks_demo_")
    suppressWarnings(suppressMessages(
      sparks_demo_results(dir, n_cells = 240L, n_genes = 120L, seed = 7L)
    ))
    .fixture$dir <- dir
  }
  .fixture$dir
}

#' A throw-away copy of the demo results, removed when the calling test ends
demo_copy <- function(env = parent.frame()) {
  src <- demo_results()
  dest <- tempfile("sparks_copy_")
  dir.create(dest)
  file.copy(list.files(src, full.names = TRUE, all.files = FALSE), dest, recursive = TRUE)
  withr::defer(unlink(dest, recursive = TRUE), envir = env)
  dest
}

demo_object <- function() {
  if (is.null(.fixture$obj)) {
    b <- sparks_read_results(demo_results())
    .fixture$obj <- readRDS(b$objects$path[[1L]])
  }
  .fixture$obj
}

demo_object_path <- function() sparks_read_results(demo_results())$objects$path[[1L]]

# ORA does not read the object, but sparks_run_spec() validates that it is a Seurat object.
demo_object_or_stub <- function() {
  skip_if_no_seurat()
  demo_object()
}
