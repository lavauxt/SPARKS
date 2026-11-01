test_that(".classify_result_file recognises the pipeline's file prefixes", {
  expect_equal(.classify_result_file("AllMarkers_All_by_cell_type.txt"), "markers_all")
  expect_equal(.classify_result_file("DEG_All_T_cells_by_cell_type.txt"), "deg_sig")
  expect_equal(.classify_result_file("DEG_Counts_All_by_cell_type.txt"), "deg_counts")
  expect_equal(.classify_result_file("GSEA_Status_x.txt"), "gsea_status")
  expect_equal(.classify_result_file("GSEA_H_x.tsv"), "gsea")
  expect_equal(.classify_result_file("PseudobulkDEG_x_by_y.txt"), "pseudobulk_deg")
  expect_equal(.classify_result_file("PseudobulkCounts_x.txt"), "pseudobulk_counts")
  expect_equal(.classify_result_file("before_filtering_S1_cell_counts.txt"), "cell_counts")
  expect_equal(.classify_result_file("UMAP_x_y.png"), "plot_umap")
  expect_equal(.classify_result_file("UMAP_notes.txt"), "other_table")
  expect_equal(.classify_result_file("qc_report.html"), "report")
  expect_equal(.classify_result_file("Merged_A.rds"), "object")
  expect_equal(.classify_result_file("something.bin"), "file")
})

test_that(".parse_result_path recovers comparison / suffix / unit / grouping", {
  p <- .parse_result_path("Cmp/Cmp/All_cells/cell_type/DEG/AllMarkers_x.txt")
  expect_equal(c(p$comparison, p$suffix, p$unit, p$grouping, p$category),
               c("Cmp", "Cmp", "All_cells", "cell_type", "DEG"))
  q <- .parse_result_path("Cmp/QC/VlnPlot_x.png")
  expect_equal(c(q$comparison, q$category), c("Cmp", "QC"))
  expect_true(is.na(q$unit))
  top <- .parse_result_path("pipeline.log")
  expect_true(is.na(top$comparison))
})

test_that("sparks_read_results errors clearly on a missing directory", {
  expect_error(sparks_read_results(file.path(tempdir(), "no-such-dir")), "does not exist")
})

test_that("the bundle indexes a results directory and skips pipeline internals", {
  dir <- demo_copy()
  dir.create(file.path(dir, "_targets", "objects"), recursive = TRUE)
  writeLines("x", file.path(dir, "_targets", "objects", "a"))
  writeLines("report", file.path(dir, "SPARKS_all_results.html"))
  b <- sparks_read_results(dir)
  expect_s3_class(b, "sparks_bundle")
  expect_false(any(grepl("_targets", b$catalog$rel)))
  expect_false("SPARKS_all_results.html" %in% b$catalog$name)
  expect_true(all(c("markers_all", "deg_sig", "gsea", "proportions", "cell_counts") %in% b$catalog$kind))
  expect_equal(nrow(b$objects), 1L)
  expect_true(all(c("comparison", "unit", "grouping") %in% names(b$units)))
  expect_output(print(b), "sparks_bundle")
})

test_that("sparks_table stacks per-group files and keeps provenance columns", {
  b <- sparks_read_results(demo_results())
  sig <- sparks_table(b, "deg_sig", grouping = "cell_type")
  expect_gt(nrow(sig), 0L)
  expect_true(all(c(".comparison", ".unit", ".grouping", ".file", "gene", "avg_log2FC") %in% names(sig)))
  expect_setequal(unique(sig$cluster), c("T cells", "B cells", "Monocytes"))
  expect_equal(nrow(sparks_table(b, "no_such_kind")), 0L)
  one <- sparks_table(b, "markers_all", grouping = "cell_type", max_rows = 5L)
  expect_lte(nrow(one), 5L)
})

test_that(".rbind_fill pads missing columns", {
  out <- .rbind_fill(list(data.frame(a = 1, b = 2), data.frame(a = 3, c = "x")))
  expect_equal(names(out), c("a", "b", "c"))
  expect_true(is.na(out$b[[2]]) && is.na(out$c[[1]]))
})

test_that("QC text files are parsed", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("Total cells: 30", "\nCells by orig.ident:", "S1: 10", "S2: 20", "",
               "Cells per cluster by orig.ident:"), f)
  p <- .parse_cell_counts(f)
  expect_equal(p$total, 30L)
  expect_equal(p$by_sample$sample, c("S1", "S2"))
  expect_equal(p$by_sample$cells, c(10L, 20L))

  g <- tempfile(fileext = ".txt")
  # what print(table(x)) writes
  writeLines(c("", "doublet singlet ", "    12     488 "), g)
  d <- .parse_doublet_stats(g)
  expect_equal(d$class, c("doublet", "singlet"))
  expect_equal(d$cells, c(12L, 488L))

  b <- sparks_read_results(demo_results())
  expect_true(all(c("before_filtering", "after_filtering") %in% b$qc$stages$stage))
  expect_setequal(unique(b$qc$doublets$class), c("doublet", "singlet"))
})

test_that("sparks_panel extracts a compact, decodable view of a Seurat object", {
  obj <- demo_object()
  p <- sparks_panel(obj, max_cells = 100L, max_genes = 10L)
  expect_s3_class(p, "sparks_panel")
  expect_length(p$cells, 100L)
  expect_equal(p$n_total, ncol(obj))
  expect_length(p$x, 100L)
  expect_lte(length(p$genes), 10L)
  ct <- .panel_meta_values(p, "cell_type")
  expect_s3_class(ct, "factor")
  expect_equal(as.character(ct), as.character(obj$cell_type[p$cells]))
  g <- names(p$genes)[[1L]]
  dense <- .panel_gene_values(p, g)
  expect_equal(unname(dense), unname(as.numeric(SeuratObject::LayerData(obj, layer = "data")[g, p$cells])),
               tolerance = 0.01)
  # down-sampling is reproducible
  expect_equal(sparks_panel(obj, max_cells = 50L, seed = 3L)$cells,
               sparks_panel(obj, max_cells = 50L, seed = 3L)$cells)
  expect_output(print(p), "sparks_panel")
  expect_error(sparks_panel(list()), "Seurat object")
})

test_that("panels can be extended with genes and metadata", {
  obj <- demo_object()
  p <- sparks_panel(obj, max_genes = 0L)
  expect_length(p$genes, 0L)
  p <- sparks_panel_add_genes(p, obj, c("CD3E", "NOT_A_GENE", "LYZ"))
  expect_setequal(names(p$genes), c("CD3E", "LYZ"))
  vals <- stats::setNames(rep(c("a", "b"), length.out = length(p$cells)), p$cells)
  p <- sparks_panel_add_meta(p, "custom", vals)
  expect_equal(p$meta$custom$levels, c("a", "b"))
  p <- sparks_panel_add_meta(p, "score", stats::setNames(seq_along(p$cells) / 7, p$cells))
  expect_equal(p$meta$score$type, "num")
  expect_error(sparks_panel_add_meta(p, "bad", 1:3), "match the panel")
})
