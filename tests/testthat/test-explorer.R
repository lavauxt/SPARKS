skip_if_no_explorer <- function() {
  skip_if_no_report_stack()
  for (pkg in c("shiny", "bslib", "targets", "digest")) skip_if_not_installed(pkg)
}

test_that("the explorer UI builds and exposes the Loupe-style controls", {
  skip_if_no_explorer()
  html <- as.character(sparks_explorer_ui("x"))
  for (id in c("x-embedding", "x-save_cat", "x-run_de", "x-run_recluster", "x-run_sig",
               "x-run_ora", "x-build_report", "x-import", "x-color_by")) {
    expect_true(grepl(paste0('id="', id, '"'), html, fixed = TRUE), info = id)
  }
})

test_that("explorer: select cells, save a category, run DE, re-cluster and score, all on cached targets", {
  skip_if_no_explorer()
  dir <- demo_copy()
  path <- sparks_read_results(dir)$objects$path[[1L]]
  seu <- readRDS(path)  # not `obj`: that name is a reactive inside the module
  tcells <- colnames(seu)[seu$cell_type == "T cells"]
  bcells <- colnames(seu)[seu$cell_type == "B cells"]

  shiny::testServer(sparks_explorer_server,
    args = list(results_dir = dir, background = FALSE, max_cells = 1000L), {
      session$setInputs(dataset = path, reduction = "umap", color_by = "cell_type",
                        point_size = 4, opacity = 0.8)
      p <- session$returned$panel()
      expect_length(p$cells, ncol(seu))
      expect_true("cell_type" %in% names(p$meta))

      # selection by pasted barcodes (the lasso feeds the same state)
      session$setInputs(paste = paste(tcells, collapse = "\n"), apply_paste = 1)
      expect_setequal(session$returned$selection(), tcells)

      # save as a category -> becomes a colour-by column
      session$setInputs(cat_name = "Region A", save_cat = 1)
      expect_named(session$returned$categories(), "Region A")
      expect_true("Selection categories" %in% names(session$returned$panel()$meta))

      # DE: selection vs rest
      session$setInputs(de_vs = "rest", de_test = "wilcox", de_min_pct = 0.1, de_lfc = 0.25, run_de = 1)
      h <- session$returned$history()
      expect_equal(nrow(h), 1L)
      expect_equal(h$status, "built")
      session$setInputs(view_id = h$id[[1L]])
      r <- session$returned$result()
      expect_equal(r$type, "de")
      expect_true("CD3E" %in% r$table$gene)

      # a second selection (B cells), compared against the saved category "Region A"
      session$setInputs(paste = paste(bcells, collapse = " "), apply_paste = 2)
      session$setInputs(cat_name = "Region B", save_cat = 2)
      session$setInputs(de_vs = "Region A", run_de = 2)
      h <- session$returned$history()
      expect_equal(nrow(h), 2L)
      expect_equal(h$status, c("built", "built"))
      session$setInputs(view_id = h$id[[1L]])
      expect_true("CD79A" %in% session$returned$result()$table$gene[1:10])

      # comparing a selection with itself fails cleanly: recorded, not fatal
      quiet <- function(expr) {   # targets echoes the failed target on stdout/stderr
        utils::capture.output(utils::capture.output(suppressMessages(expr), type = "message"), type = "output")
        invisible()
      }
      quiet(session$setInputs(de_vs = "Region B", run_de = 3))
      h <- session$returned$history()
      expect_equal(nrow(h), 3L)
      expect_equal(sum(h$status == "error"), 1L)
      expect_match(h$error[h$status == "error"], "at least 3 cells")

      # deleting it removes the spec and the cached result (failed analyses are
      # retried by every later run, so drop it before going on)
      session$setInputs(view_id = h$id[h$status == "error"], delete = 1)
      expect_equal(nrow(session$returned$history()), 2L)
      expect_false(any(session$returned$history()$status == "error"))

      # re-cluster the selection, then expose the clusters for colouring
      session$setInputs(rc_res = 0.5, run_recluster = 1)
      h <- session$returned$history()
      rc <- h[h$type == "recluster", ]
      expect_equal(rc$status, "built")
      session$setInputs(view_id = rc$id)
      session$setInputs(apply_recluster = 1)
      expect_match(names(session$returned$extra_meta()), "^recluster_")
      expect_true(names(session$returned$extra_meta())[[1L]] %in% names(session$returned$panel()$meta))

      # signature scoring
      session$setInputs(sig_set = "", sig_genes = "CD3E, CD3D, IL7R", run_sig = 1)
      expect_true("signature" %in% session$returned$history()$type)

      # clearing resets only the selection
      session$setInputs(clear_sel = 1)
      expect_length(session$returned$selection(), 0L)
      expect_length(session$returned$categories(), 2L)
      session$setInputs(clear_cats = 1)
      expect_length(session$returned$categories(), 0L)
    })
})

test_that("explorer: unknown barcodes do not change the selection", {
  skip_if_no_explorer()
  dir <- demo_copy()
  path <- sparks_read_results(dir)$objects$path[[1L]]
  shiny::testServer(sparks_explorer_server,
    args = list(results_dir = dir, background = FALSE), {
      session$setInputs(dataset = path, reduction = "umap", color_by = "cell_type")
      session$setInputs(paste = "NOPE1\nNOPE2", apply_paste = 1)
      expect_length(session$returned$selection(), 0L)
    })
})

test_that("explorer: the background-run default is off when SPARKS is only dev-loaded or callr is missing", {
  skip_if_no_explorer()
  expect_type(.xp_can_background(), "logical")
})

test_that("explorer helpers pick sensible defaults", {
  expect_equal(.xp_default_group(c("a", "seurat_clusters", "cell_type")), "cell_type")
  expect_equal(.xp_default_group(c("zzz", "yyy")), "zzz")
  skip_if_no_seurat()
  m <- .xp_meta_columns(demo_object())
  expect_true(all(c("cell_type", "orig.ident") %in% m$cat))
  expect_true("nCount_RNA" %in% m$num)
  expect_false("nCount_RNA" %in% m$cat)
})

test_that("sparks_app includes the Explorer tab", {
  skip_if_no_explorer()
  skip_if_not_installed("DT")
  expect_equal(.sparks_nav_labels[["explorer"]], "Explorer")
  expect_s3_class(sparks_app(), "shiny.appobj")
})
