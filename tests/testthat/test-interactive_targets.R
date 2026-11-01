skip_if_no_targets <- function() {
  skip_if_no_seurat()
  skip_if_not_installed("targets")
  skip_if_not_installed("digest")
}

prep_store <- function(env = parent.frame()) {
  dir <- demo_copy(env)
  list(dir = dir, object = sparks_read_results(dir)$objects$path[[1L]])
}

test_that("an analysis is queued, run as a target, and read back", {
  skip_if_no_targets()
  s <- prep_store()
  obj <- readRDS(s$object)
  tc <- colnames(obj)[obj$cell_type == "T cells"]
  id <- sparks_interactive_submit(s$dir, sparks_spec("de", cells_a = tc, label = "T vs rest"), s$object)
  expect_true(file.exists(file.path(s$dir, "_targets_interactive", "specs", paste0(id, ".rds"))))
  st <- sparks_interactive_status(s$dir)
  expect_equal(st$status, "pending")

  st <- sparks_interactive_run(s$dir, s$object)
  expect_equal(st$id, id)
  expect_equal(st$status, "built")
  expect_equal(st$label, "T vs rest")
  res <- sparks_interactive_result(s$dir, id)
  expect_s3_class(res, "sparks_analysis")
  expect_true("CD3E" %in% res$table$gene)
})

test_that("asking for the same analysis again is a cache hit", {
  skip_if_no_targets()
  s <- prep_store()
  obj <- readRDS(s$object)
  tc <- colnames(obj)[obj$cell_type == "T cells"]
  spec <- sparks_spec("de", cells_a = tc)
  id1 <- sparks_interactive_submit(s$dir, spec, s$object)
  first <- sparks_interactive_run(s$dir, s$object)
  id2 <- sparks_interactive_submit(s$dir, spec, s$object)
  second <- sparks_interactive_run(s$dir, s$object)
  expect_equal(id1, id2)
  expect_equal(nrow(second), 1L)
  expect_equal(second$seconds, first$seconds)   # not re-run: targets kept the recorded run time
})

test_that("a failing analysis is recorded without blocking the others", {
  skip_if_no_targets()
  s <- prep_store()
  obj <- readRDS(s$object)
  tc <- colnames(obj)[obj$cell_type == "T cells"]
  bad <- sparks_interactive_submit(s$dir, sparks_spec("signature", genes = c("nope1", "nope2")), s$object)
  good <- sparks_interactive_submit(s$dir, sparks_spec("de", cells_a = tc), s$object)
  st <- suppressMessages(sparks_interactive_run(s$dir, s$object))
  expect_equal(st$status[st$id == good], "built")
  expect_equal(st$status[st$id == bad], "error")
  expect_match(st$error[st$id == bad], "at least 2 genes")
})

test_that("the same selection on two objects gets two ids", {
  skip_if_no_targets()
  s <- prep_store()
  other <- tempfile(fileext = ".rds"); file.copy(s$object, other)
  spec <- sparks_spec("de", cells_a = c("a", "b", "c"))
  expect_false(identical(sparks_interactive_submit(s$dir, spec, s$object),
                         sparks_interactive_submit(s$dir, spec, other)))
})

test_that("analyses can be deleted", {
  skip_if_no_targets()
  s <- prep_store()
  obj <- readRDS(s$object)
  id <- sparks_interactive_submit(s$dir, sparks_spec("de", cells_a = colnames(obj)[obj$cell_type == "T cells"]), s$object)
  sparks_interactive_run(s$dir, s$object)
  expect_equal(nrow(sparks_interactive_status(s$dir)), 1L)
  sparks_interactive_delete(s$dir, id)
  expect_equal(nrow(sparks_interactive_status(s$dir)), 0L)
})

test_that("an empty store reports an empty history", {
  skip_if_no_targets()
  dir <- demo_copy()
  expect_equal(nrow(sparks_interactive_status(dir)), 0L)
  expect_named(sparks_interactive_status(dir), c("id", "type", "label", "created", "status", "error", "seconds", "object"))
})

test_that("interactive analyses are separate from the pipeline's own _targets store", {
  skip_if_no_targets()
  s <- prep_store()
  sparks_interactive_submit(s$dir, sparks_spec("de", cells_a = c("a", "b", "c")), s$object)
  expect_equal(dirname(.interactive_paths(s$dir)$dir), s$dir)
  expect_false(identical(.interactive_paths(s$dir)$store, file.path(s$dir, "_targets")))
})
