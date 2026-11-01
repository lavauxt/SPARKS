test_that("specs are content-addressed and order-insensitive", {
  a <- sparks_spec("de", cells_a = c("c3", "c1", "c2"))
  b <- sparks_spec("de", cells_a = c("c1", "c2", "c3", "c1"))
  expect_equal(a$id, b$id)
  expect_false(identical(a$id, sparks_spec("de", cells_a = c("c1", "c2"))$id))
  expect_false(identical(a$id, sparks_spec("de", cells_a = c("c1", "c2", "c3"), test = "t")$id))
  # the label and creation time are presentation, not identity
  expect_equal(a$id, sparks_spec("de", cells_a = c("c1", "c2", "c3"), label = "renamed")$id)
  expect_equal(sparks_spec("de", cells_a = "c1", label = "mine")$label, "mine")
  expect_match(sparks_spec("recluster", cells = letters)$label, "Re-cluster 26 cells")
  expect_output(print(a), "sparks_spec")
  expect_error(sparks_spec("nope"), "should be one of")
})

test_that("differential expression recovers the planted markers", {
  skip_if_no_seurat()
  obj <- demo_object()
  tc <- colnames(obj)[obj$cell_type == "T cells"]
  res <- sparks_run_spec(sparks_spec("de", cells_a = tc), obj)
  expect_s3_class(res, "sparks_analysis")
  expect_equal(res$type, "de")
  expect_equal(res$n_a, length(tc))
  expect_equal(res$n_b, ncol(obj) - length(tc))
  expect_true(all(c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj") %in% names(res$table)))
  top <- res$table$gene[res$table$avg_log2FC > 0 & res$table$p_val_adj < 0.05]
  expect_true(all(c("CD3E", "CD3D") %in% top))
  expect_true("CD79A" %in% res$table$gene[res$table$avg_log2FC < 0])
  expect_false(is.unsorted(res$table$p_val_adj))
})

test_that("DE against an explicit group excludes shared cells and validates sizes", {
  skip_if_no_seurat()
  obj <- demo_object()
  tc <- colnames(obj)[obj$cell_type == "T cells"]
  bc <- colnames(obj)[obj$cell_type == "B cells"]
  res <- sparks_run_spec(sparks_spec("de", cells_a = tc, cells_b = bc, min_pct = 0.05), obj)
  expect_equal(res$n_b, length(bc))
  overlap <- sparks_run_spec(sparks_spec("de", cells_a = tc, cells_b = c(bc, tc[1:3])), obj)
  expect_equal(overlap$n_b, length(bc))
  expect_error(sparks_run_spec(sparks_spec("de", cells_a = tc[1:2]), obj), "at least 3 cells")
  expect_error(sparks_run_spec(sparks_spec("de", cells_a = tc, cells_b = bc[1:2]), obj), "at least 3 cells")
  expect_error(sparks_run_spec(sparks_spec("de", cells_a = c("nope1", "nope2", "nope3")), obj), "at least 3 cells")
  expect_error(sparks_run_spec(sparks_spec("de", cells_a = tc), "not an object"), "No such file|cannot open|does not exist")
})

test_that("an object whose data layer was never normalised is normalised before DE", {
  skip_if_no_seurat()
  obj <- demo_object()
  raw <- SeuratObject::CreateSeuratObject(SeuratObject::LayerData(obj, layer = "counts"))
  raw$cell_type <- obj$cell_type[colnames(raw)]
  fixed <- .ensure_log_data(raw, "RNA")
  d <- SeuratObject::LayerData(fixed, layer = "data")
  k <- SeuratObject::LayerData(fixed, layer = "counts")
  expect_false(isTRUE(all.equal(as.numeric(d[1:20, 1:20]), as.numeric(k[1:20, 1:20]))))
  tc <- colnames(raw)[raw$cell_type == "T cells"]
  expect_true("CD3E" %in% sparks_run_spec(sparks_spec("de", cells_a = tc), raw)$table$gene)
})

test_that("re-clustering a subset returns clusters, an embedding and markers", {
  skip_if_no_seurat()
  obj <- demo_object()
  cells <- colnames(obj)[obj$cell_type %in% c("T cells", "B cells")]
  res <- sparks_run_spec(sparks_spec("recluster", cells = cells, resolution = 0.5, npcs = 8L, nfeatures = 60L), obj)
  expect_equal(sort(names(res$clusters)), sort(cells))
  expect_equal(dim(res$embedding), c(length(cells), 2L))
  expect_gte(res$n_clusters, 2L)
  # the two planted populations should end up in different clusters
  tab <- table(obj$cell_type[names(res$clusters)], res$clusters)
  expect_true(all(apply(tab[c("T cells", "B cells"), , drop = FALSE], 1L, function(r) max(r) / sum(r)) > 0.9))
  expect_s3_class(res$markers, "data.frame")
  expect_error(sparks_run_spec(sparks_spec("recluster", cells = cells[1:10]), obj), "at least 30 cells")
})

test_that("signature scoring separates the cell type the genes mark", {
  skip_if_no_seurat()
  obj <- demo_object()
  genes <- c("CD3E", "CD3D", "IL7R", "CD8A")
  res <- sparks_run_spec(sparks_spec("signature", genes = c(genes, "NOT_A_GENE"), name = "T sig"), obj)
  expect_equal(res$missing, "NOT_A_GENE")
  expect_setequal(res$genes, genes)
  expect_equal(names(res$scores), colnames(obj))
  by <- tapply(res$scores, obj$cell_type, mean)
  expect_equal(names(which.max(by)), "T cells")
  mean_res <- sparks_run_spec(sparks_spec("signature", genes = genes, method = "mean"), obj)
  expect_equal(mean_res$method, "mean")
  expect_gt(stats::cor(res$scores, mean_res$scores), 0.9)
  expect_error(sparks_run_spec(sparks_spec("signature", genes = c("A", "B")), obj), "at least 2 genes")
})

test_that("over-representation uses the hypergeometric tail", {
  sets <- list(setA = paste0("g", 1:10), setB = paste0("g", 11:20))
  universe <- paste0("g", 1:100)
  res <- sparks_run_spec(sparks_spec("ora", genes = paste0("g", 1:5), universe = universe, gene_sets = sets), object = demo_object_or_stub())
  tab <- res$table
  expect_equal(tab$pathway[[1L]], "setA")
  expect_equal(tab$overlap[tab$pathway == "setA"], 5L)
  expect_equal(tab$pval[tab$pathway == "setA"], stats::phyper(4, 10, 90, 5, lower.tail = FALSE))
  expect_equal(tab$padj, stats::p.adjust(tab$pval, "BH"))
  expect_error(sparks_run_spec(sparks_spec("ora", genes = "x"), object = demo_object_or_stub()), "gene set")
  expect_error(sparks_run_spec(sparks_spec("ora", genes = "zzz", universe = universe, gene_sets = sets),
                               object = demo_object_or_stub()), "universe")
})

test_that("GMT files and config signatures become gene sets", {
  gmt <- tempfile(fileext = ".gmt")
  writeLines(c("SET_A\tdesc\tG1\tG2\tG3", "SET_B\tdesc\tG4\tG5", ""), gmt)
  sets <- sparks_read_gmt(gmt)
  expect_equal(names(sets), c("SET_A", "SET_B"))
  expect_equal(sets$SET_A, c("G1", "G2", "G3"))
  cfg <- list(gene_signatures = list(list(name = "Sig1", genes = list("A", "B")),
                                     list(name = "Empty", genes = list())))
  expect_equal(sparks_config_gene_sets(cfg), list(Sig1 = c("A", "B")))
  expect_equal(sparks_config_gene_sets(list()), list())
})

test_that("barcodes can be read from the report's CSV, a header-less list or a TSV", {
  csv <- tempfile(fileext = ".csv")
  writeLines(c("barcode,cell_type", "C1,T", "C2,B", "C1,T"), csv)
  expect_equal(sparks_read_barcodes(csv), c("C1", "C2"))
  txt <- tempfile(fileext = ".txt")
  writeLines(c("C9", "", "C8"), txt)
  expect_equal(sparks_read_barcodes(txt), c("C9", "C8"))
  tsv <- tempfile(fileext = ".tsv")
  writeLines(c("sample\tCell\tx", "s\tC5\t1"), tsv)
  expect_equal(sparks_read_barcodes(tsv), "C5")
  empty <- tempfile(); writeLines(character(), empty)
  expect_equal(sparks_read_barcodes(empty), character())
})

test_that("categories become a factor aligned to cells; the last category wins", {
  cells <- paste0("c", 1:6)
  f <- sparks_categories_to_factor(list(A = c("c1", "c2"), B = c("c2", "c3")), cells)
  expect_equal(as.character(f), c("A", "B", "B", "Unassigned", "Unassigned", "Unassigned"))
  expect_equal(levels(f), c("A", "B", "Unassigned"))
  all_assigned <- sparks_categories_to_factor(list(A = cells), cells)
  expect_equal(levels(all_assigned), "A")
})
