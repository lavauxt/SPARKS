test_that("sparks_colors keeps the validated palette in fixed order", {
  light <- sparks_colors(8L)
  expect_equal(light[1:3], c("#2a78d6", "#eb6834", "#1baf7a"))
  expect_equal(sparks_colors(8L, "dark")[1L], "#3987e5")
  # fixed order: the first n colours do not depend on how many are requested
  expect_equal(sparks_colors(3L), sparks_colors(12L)[1:3])
  many <- sparks_colors(20L)
  expect_length(many, 20L)
  expect_true(all(grepl("^#[0-9a-fA-F]{6}$", many)))
  expect_equal(length(unique(many)), 20L)
  expect_length(sparks_colors(0L), 0L)
})

test_that("the theme map covers every palette colour and swaps light for dark", {
  m <- .sparks_theme_map()
  expect_equal(unname(m[["#2a78d6"]]), "#3987e5")
  expect_true(all(c("#fcfcfb", "#0b0b0b") %in% names(m)))
  expect_equal(unname(m[["#fcfcfb"]]), "#1a1a19")
})

test_that("figure builders return NULL on empty input instead of failing", {
  skip_if_not_installed("plotly")
  expect_null(sparks_fig_qc_stages(data.frame()))
  expect_null(sparks_fig_doublets(NULL))
  expect_null(sparks_fig_deg_counts(data.frame(x = 1)))
  expect_null(sparks_fig_volcano(data.frame()))
  expect_null(sparks_fig_composition(data.frame(a = 1)))
  expect_null(sparks_fig_gsea(data.frame()))
  expect_null(sparks_fig_heatmap(matrix(numeric(), 0, 0)))
  expect_null(sparks_fig_cells_per_sample(data.frame()))
})

test_that("results figures are plotly objects with a serialisable spec", {
  skip_if_no_report_stack()
  b <- sparks_read_results(demo_results())
  figs <- list(
    sparks_fig_qc_stages(b$qc$stages),
    sparks_fig_doublets(b$qc$doublets),
    sparks_fig_deg_counts(sparks_table(b, "deg_counts", grouping = "cell_type")),
    sparks_fig_volcano(sparks_table(b, "markers_all", grouping = "cell_type")),
    sparks_fig_gsea(sparks_table(b, "gsea")),
    sparks_fig_cells_per_sample(data.frame(sample = c("a", "b"), cells = c(3L, 4L))),
    sparks_fig_heatmap(matrix(stats::rnorm(6), 2, 3, dimnames = list(c("a", "b"), c("x", "y", "z"))), diverging = TRUE)
  )
  for (f in figs) {
    expect_s3_class(f, "plotly")
    spec <- jsonlite::fromJSON(.figure_spec_json(f), simplifyVector = FALSE)
    expect_true(all(c("data", "layout") %in% names(spec)))
    expect_gt(length(spec$data), 0L)
  }
  # a single-bar chart must still serialise x / y as arrays (plotly.js rejects scalars)
  one <- .figure_spec_json(sparks_fig_deg_counts(data.frame(Group = "a", DEG_Count = 3L)))
  expect_match(one, "\"x\":\\[3\\]")
})

test_that("the volcano classifies genes with the requested thresholds", {
  skip_if_not_installed("plotly")
  deg <- data.frame(gene = c("up", "down", "ns", "weak"), avg_log2FC = c(2, -2, 3, 0.1),
                    p_val_adj = c(1e-5, 1e-5, 0.5, 1e-5))
  spec <- jsonlite::fromJSON(.figure_spec_json(sparks_fig_volcano(deg, lfc_cut = 0.25, p_cut = 0.05)),
                             simplifyVector = FALSE)
  names_ <- vapply(spec$data, function(t) t$name %||% "", character(1))
  expect_true(any(grepl("^Up \\(1\\)", names_)))
  expect_true(any(grepl("^Down \\(1\\)", names_)))
  expect_true(any(grepl("^Not significant \\(2\\)", names_)))
  # alternative column names are picked up (pseudobulk tables use logFC)
  pb <- data.frame(gene = "g", logFC = 1, padj = 0.001)
  expect_s3_class(sparks_fig_volcano(pb), "plotly")
})

test_that("composition and GSEA figures stack the right numbers", {
  skip_if_not_installed("plotly")
  prop <- data.frame(Var1 = rep(c("T", "B"), 2), Var2 = rep(c("s1", "s2"), each = 2),
                     Freq = c(0.75, 0.25, 0.5, 0.5), Count = c(3L, 1L, 2L, 2L))
  spec <- jsonlite::fromJSON(.figure_spec_json(sparks_fig_composition(prop)), simplifyVector = FALSE)
  expect_equal(spec$layout$barmode, "stack")
  tot <- Reduce(`+`, lapply(spec$data, function(t) unlist(t$y)))
  expect_equal(tot, c(100, 100))
})

test_that("the dot plot summary matches a direct computation", {
  skip_if_no_report_stack()
  obj <- demo_object()
  p <- sparks_panel(obj, max_genes = 0L)
  p <- sparks_panel_add_genes(p, obj, c("CD3E", "LYZ"))
  d <- sparks_dot_data(p, c("CD3E", "LYZ", "missing"), "cell_type")
  expect_setequal(unique(d$gene), c("CD3E", "LYZ"))
  expect_equal(nrow(d), 2L * nlevels(.panel_meta_values(p, "cell_type")))
  v <- .panel_gene_values(p, "CD3E")
  g <- .panel_meta_values(p, "cell_type")
  expect_equal(d$mean_expr[d$gene == "CD3E" & d$group == "T cells"], mean(v[g == "T cells"]))
  expect_equal(d$pct_expr[d$gene == "CD3E" & d$group == "T cells"], 100 * mean(v[g == "T cells"] > 0))
  # the T-cell marker is highest in T cells
  expect_equal(d$group[d$gene == "CD3E"][which.max(d$mean_expr[d$gene == "CD3E"])], "T cells")
  expect_equal(nrow(sparks_dot_data(p, "missing", "cell_type")), 0L)
})

test_that("embedding and violin figures accept panel genes or external values", {
  skip_if_no_report_stack()
  obj <- demo_object()
  p <- sparks_panel(obj, max_genes = 0L)
  expect_s3_class(sparks_fig_embedding(p, "cell_type"), "plotly")
  expect_s3_class(sparks_fig_embedding(p, "nCount_RNA"), "plotly")
  vals <- stats::runif(length(p$cells))
  e <- sparks_fig_embedding(p, "ignored", values = vals, value_label = "my score")
  spec <- jsonlite::fromJSON(.figure_spec_json(e), simplifyVector = FALSE)
  expect_equal(spec$data[[1]]$marker$colorbar$title$text, "my score")
  sel <- sparks_fig_embedding(p, "cell_type", selected = p$cells[1:5])
  expect_true(any(grepl("^Selection \\(5\\)", vapply(plotly::plotly_build(sel)$x$data, function(t) t$name %||% "", character(1)))))
  expect_s3_class(sparks_fig_violin(p, "x", "cell_type", values = vals), "plotly")
  expect_error(sparks_fig_violin(p, "not_in_panel", "cell_type"), "not in the panel")
})
