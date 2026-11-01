# Re-analysis primitives for the interactive explorer.
#
# An *analysis spec* is a small, plain list describing one on-the-fly analysis
# (differential expression between two cell sets, re-clustering a selection,
# scoring a gene signature, over-representation of a gene list). Specs are
# hashable, so the explorer can run them as cached targets: asking for the
# same analysis twice costs nothing, and the history survives restarts.
#
# Everything here is a pure function of (spec, Seurat object): no Shiny, no
# targets, no global state -- which is what makes it unit-testable.

# ---- Specs -------------------------------------------------------------------------

#' Describe an on-the-fly analysis
#'
#' @param type One of \code{"de"} (differential expression between two cell
#'   sets), \code{"recluster"} (re-cluster a subset), \code{"signature"} (score
#'   a gene set per cell) or \code{"ora"} (over-representation of a gene list).
#' @param ... Type-specific parameters:
#'   \describe{
#'     \item{de}{\code{cells_a}, \code{cells_b} (\code{NULL} = all other cells),
#'       \code{test}, \code{min_pct}, \code{logfc_threshold}, \code{assay}}
#'     \item{recluster}{\code{cells}, \code{resolution}, \code{nfeatures},
#'       \code{npcs}, \code{markers}, \code{seed}}
#'     \item{signature}{\code{genes}, \code{method} (\code{"module"} or
#'       \code{"mean"}), \code{name}}
#'     \item{ora}{\code{genes}, \code{universe}, \code{gene_sets} (named list)}
#'   }
#' @param label Human-readable label shown in the analysis history.
#' @return A \code{sparks_spec}: a list with an \code{id} (content hash).
#' @export
#' @examples
#' spec <- sparks_spec("de", cells_a = c("c1", "c2", "c3"), label = "Cluster 1 vs rest")
#' spec$id
sparks_spec <- function(type = c("de", "recluster", "signature", "ora"), ..., label = NULL) {
  type <- match.arg(type)
  params <- list(...)
  defaults <- switch(
    type,
    de = list(cells_a = character(), cells_b = NULL, test = "wilcox", min_pct = 0.1,
              logfc_threshold = 0.25, assay = NULL),
    recluster = list(cells = character(), resolution = 0.5, nfeatures = 2000L, npcs = 30L,
                     markers = TRUE, seed = 1L),
    signature = list(genes = character(), method = "module", name = "Signature"),
    ora = list(genes = character(), universe = NULL, gene_sets = list())
  )
  params <- utils::modifyList(defaults, params[!vapply(params, is.null, logical(1))])
  for (nm in names(defaults)) if (is.null(params[[nm]]) && is.null(defaults[[nm]])) params[nm] <- list(NULL)
  spec <- c(list(type = type), params)
  for (k in c("cells_a", "cells_b", "cells", "genes")) {
    if (!is.null(spec[[k]])) spec[[k]] <- sort(unique(as.character(spec[[k]])))
  }
  spec$id <- .spec_hash(spec)
  spec$label <- label %||% .default_spec_label(spec)
  spec$created <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  structure(spec, class = c("sparks_spec", "list"))
}

.spec_hash <- function(spec) {
  core <- spec[setdiff(names(spec), c("id", "label", "created"))]
  digest::digest(core, algo = "xxhash64")
}

.default_spec_label <- function(spec) {
  switch(
    spec$type,
    de = sprintf("DE: %d cells vs %s", length(spec$cells_a),
                 if (is.null(spec$cells_b)) "rest" else paste(length(spec$cells_b), "cells")),
    recluster = sprintf("Re-cluster %d cells (res %.2g)", length(spec$cells), spec$resolution),
    signature = sprintf("Score: %s (%d genes)", spec$name, length(spec$genes)),
    ora = sprintf("ORA: %d genes", length(spec$genes))
  )
}

#' @export
print.sparks_spec <- function(x, ...) {
  cat("<sparks_spec ", x$id, "> ", x$label, "\n", sep = "")
  invisible(x)
}

#' Run an analysis spec on a Seurat object
#'
#' @param spec A [sparks_spec()].
#' @param object A Seurat object (or a path to an \code{.rds} file).
#' @return A list of class \code{sparks_analysis} with \code{type}, \code{spec}
#'   and type-specific results (see [sparks_spec()]).
#' @export
sparks_run_spec <- function(spec, object) {
  if (is.character(object) && length(object) == 1L) object <- .load_object_cached(object)
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  res <- switch(
    spec$type,
    de = .an_de(spec, object),
    recluster = .an_recluster(spec, object),
    signature = .an_signature(spec, object),
    ora = .an_ora(spec),
    stop("Unknown analysis type: ", spec$type, call. = FALSE)
  )
  res$type <- spec$type
  res$spec <- unclass(spec)
  structure(res, class = c("sparks_analysis", "list"))
}

.object_cache <- new.env(parent = emptyenv())

# Loading a merged object is the slow step; keep the last one per path in this session.
.load_object_cached <- function(path) {
  info <- file.info(path)
  key <- paste(normalizePath(path, winslash = "/", mustWork = TRUE), info$mtime, info$size)
  if (!is.null(.object_cache$key) && identical(.object_cache$key, key)) return(.object_cache$obj)
  obj <- readRDS(path)
  .object_cache$key <- key
  .object_cache$obj <- obj
  obj
}

# ---- Helpers -----------------------------------------------------------------------

# Choose the assay used for ad-hoc work: RNA when present (LogNormalize-able), else default.
.analysis_assay <- function(object, assay = NULL) {
  assays <- SeuratObject::Assays(object)
  if (!is.null(assay) && assay %in% assays) return(assay)
  if ("RNA" %in% assays) "RNA" else SeuratObject::DefaultAssay(object)
}

# Make sure the "data" layer holds log-normalised values (pipeline objects can carry SCT only).
.ensure_log_data <- function(object, assay) {
  layers <- tryCatch(SeuratObject::Layers(object[[assay]]), error = function(e) character())
  has_data <- "data" %in% layers
  needs <- !has_data
  if (has_data && "counts" %in% layers) {
    d <- SeuratObject::LayerData(object, assay = assay, layer = "data")
    k <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
    # Identical to counts => never normalised.
    probe <- seq_len(min(500L, nrow(d)))
    needs <- isTRUE(all.equal(as.numeric(d[probe, seq_len(min(50L, ncol(d)))]),
                              as.numeric(k[probe, seq_len(min(50L, ncol(k)))])))
  }
  if (needs) {
    object <- suppressWarnings(Seurat::NormalizeData(object, assay = assay, verbose = FALSE))
  }
  object
}

.check_cells <- function(cells, object, what, min_cells = 3L) {
  found <- intersect(cells, colnames(object))
  if (length(found) < min_cells) {
    stop(sprintf("%s needs at least %d cells present in the object (%d of %d given were found).",
                 what, min_cells, length(found), length(cells)), call. = FALSE)
  }
  found
}

# ---- Differential expression ------------------------------------------------------------

.an_de <- function(spec, object) {
  assay <- .analysis_assay(object, spec$assay)
  a <- .check_cells(spec$cells_a, object, "Group A")
  b <- if (is.null(spec$cells_b)) setdiff(colnames(object), a)
       else setdiff(.check_cells(spec$cells_b, object, "Group B"), a)
  if (length(b) < 3L) stop("Group B needs at least 3 cells that are not in group A.", call. = FALSE)
  object <- .ensure_log_data(object, assay)
  SeuratObject::DefaultAssay(object) <- assay
  # FindMarkers(cells.1 =, cells.2 =) only exists on the matrix method; for a Seurat
  # object the supported route is two idents on a throw-away metadata column.
  sub <- subset(object, cells = c(a, b))
  sub$sparks_group <- ifelse(colnames(sub) %in% a, "A", "B")
  SeuratObject::Idents(sub) <- "sparks_group"
  res <- suppressWarnings(suppressMessages(Seurat::FindMarkers(
    sub, ident.1 = "A", ident.2 = "B", assay = assay, test.use = spec$test,
    min.pct = spec$min_pct, logfc.threshold = spec$logfc_threshold,
    only.pos = FALSE, verbose = FALSE)))
  if (!is.data.frame(res) || !nrow(res)) {
    table <- data.frame(gene = character(), avg_log2FC = numeric(), pct.1 = numeric(),
                        pct.2 = numeric(), p_val = numeric(), p_val_adj = numeric())
  } else {
    res$gene <- rownames(res)
    fc <- intersect(c("avg_log2FC", "avg_logFC"), names(res))[1L]
    if (!identical(fc, "avg_log2FC")) names(res)[names(res) == fc] <- "avg_log2FC"
    table <- res[order(res$p_val_adj, -abs(res$avg_log2FC)),
                 c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"), drop = FALSE]
    rownames(table) <- NULL
  }
  list(table = table, n_a = length(a), n_b = length(b), assay = assay)
}

# ---- Re-clustering ------------------------------------------------------------------------

.an_recluster <- function(spec, object) {
  cells <- .check_cells(spec$cells, object, "Re-clustering", min_cells = 30L)
  assay <- .analysis_assay(object)
  sub <- subset(object, cells = cells)
  SeuratObject::DefaultAssay(sub) <- assay
  sub <- suppressWarnings(suppressMessages({
    sub <- Seurat::NormalizeData(sub, assay = assay, verbose = FALSE)
    nf <- min(as.integer(spec$nfeatures), nrow(sub) - 1L)
    sub <- Seurat::FindVariableFeatures(sub, assay = assay, nfeatures = nf, verbose = FALSE)
    sub <- Seurat::ScaleData(sub, assay = assay, features = SeuratObject::VariableFeatures(sub),
                             verbose = FALSE)
    npcs <- max(2L, min(as.integer(spec$npcs), ncol(sub) - 1L, length(SeuratObject::VariableFeatures(sub)) - 1L))
    sub <- Seurat::RunPCA(sub, assay = assay, npcs = npcs, verbose = FALSE, seed.use = spec$seed)
    dims <- seq_len(min(npcs, 20L))
    sub <- Seurat::FindNeighbors(sub, dims = dims, verbose = FALSE)
    sub <- Seurat::FindClusters(sub, resolution = spec$resolution, random.seed = spec$seed,
                                verbose = FALSE)
    sub
  }))
  emb <- tryCatch(
    suppressWarnings(suppressMessages(
      Seurat::Embeddings(Seurat::RunUMAP(sub, dims = seq_len(min(npcs, 20L)), seed.use = spec$seed,
                                         n.neighbors = min(30L, ncol(sub) - 1L), verbose = FALSE),
                         "umap"))),
    error = function(e) Seurat::Embeddings(sub, "pca")[, 1:2, drop = FALSE])
  clusters <- stats::setNames(as.character(sub$seurat_clusters), colnames(sub))
  markers <- NULL
  if (isTRUE(spec$markers) && length(unique(clusters)) > 1L) {
    Seurat::Idents(sub) <- "seurat_clusters"
    markers <- tryCatch(suppressWarnings(suppressMessages(Seurat::FindAllMarkers(
      sub, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
      max.cells.per.ident = 300L, verbose = FALSE))), error = function(e) NULL)
    if (is.data.frame(markers) && nrow(markers)) {
      markers <- markers[, intersect(c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"), names(markers))]
      rownames(markers) <- NULL
    } else markers <- NULL
  }
  list(clusters = clusters, embedding = emb[, 1:2, drop = FALSE], markers = markers,
       n_clusters = length(unique(clusters)), resolution = spec$resolution)
}

# ---- Signature scoring ---------------------------------------------------------------------

.an_signature <- function(spec, object) {
  assay <- .analysis_assay(object)
  genes <- intersect(spec$genes, rownames(object))
  if (length(genes) < 2L) {
    stop(sprintf("Signature scoring needs at least 2 genes present in the object (%d of %d found).",
                 length(genes), length(spec$genes)), call. = FALSE)
  }
  object <- .ensure_log_data(object, assay)
  SeuratObject::DefaultAssay(object) <- assay
  method <- spec$method
  scores <- NULL
  if (identical(method, "module")) {
    scores <- tryCatch({
      nm <- "SPARKSsig"
      o <- suppressWarnings(suppressMessages(Seurat::AddModuleScore(
        object, features = list(genes), name = nm, assay = assay,
        ctrl = max(5L, min(50L, floor(nrow(object) / 10L))), seed = 1L, nbin = min(24L, max(5L, floor(nrow(object) / 20L))))))
      stats::setNames(as.numeric(o[[paste0(nm, "1")]][, 1L]), colnames(o))
    }, error = function(e) NULL)
    if (is.null(scores)) method <- "mean"
  }
  if (is.null(scores)) {
    m <- as.matrix(SeuratObject::LayerData(object, assay = assay, layer = "data")[genes, , drop = FALSE])
    z <- t(scale(t(m)))
    z[is.na(z)] <- 0
    scores <- stats::setNames(colMeans(z), colnames(m))
  }
  list(scores = scores, genes = genes, missing = setdiff(spec$genes, genes),
       method = method, name = spec$name)
}

# ---- Over-representation ------------------------------------------------------------------------

.an_ora <- function(spec) {
  sets <- spec$gene_sets
  if (!length(sets)) stop("Provide at least one gene set (named list).", call. = FALSE)
  genes <- unique(spec$genes)
  universe <- unique(spec$universe %||% unlist(sets, use.names = FALSE))
  genes <- intersect(genes, universe)
  if (!length(genes)) stop("None of the genes are in the universe.", call. = FALSE)
  rows <- lapply(names(sets), function(nm) {
    s <- intersect(sets[[nm]], universe)
    k <- length(intersect(genes, s))
    p <- if (length(s)) stats::phyper(k - 1L, length(s), length(universe) - length(s),
                                      length(genes), lower.tail = FALSE) else NA_real_
    data.frame(pathway = nm, overlap = k, set_size = length(s), n_genes = length(genes),
               fold_enrichment = if (length(s)) (k / length(genes)) / (length(s) / length(universe)) else NA_real_,
               pval = p, overlap_genes = paste(intersect(genes, s), collapse = ","),
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab$padj <- stats::p.adjust(tab$pval, "BH")
  tab <- tab[order(tab$pval), , drop = FALSE]
  rownames(tab) <- NULL
  list(table = tab, n_genes = length(genes), universe = length(universe))
}

#' Read a GMT file into a named list of gene sets
#' @param path Path to a \code{.gmt} file.
#' @return Named list of character vectors.
#' @export
sparks_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  sets <- lapply(parts, function(p) unique(p[-(1:2)][nzchar(p[-(1:2)])]))
  names(sets) <- vapply(parts, `[[`, character(1), 1L)
  sets
}

#' Gene sets defined in the pipeline configuration
#' @param cfg Pipeline config list (e.g. \code{bundle$config}).
#' @return Named list of character vectors (from \code{gene_signatures}).
#' @export
sparks_config_gene_sets <- function(cfg) {
  sigs <- cfg$gene_signatures
  if (!length(sigs)) return(list())
  out <- lapply(sigs, function(s) unlist(s$genes, use.names = FALSE))
  names(out) <- vapply(sigs, function(s) as.character(s$name %||% "signature"), character(1))
  out[vapply(out, length, integer(1)) > 0L]
}

# ---- Categories (named cell sets) ----------------------------------------------------------------

#' Turn named cell sets into a metadata vector
#' @param categories Named list: label -> barcodes.
#' @param cells All cell barcodes (defines the output order).
#' @param unassigned Label for cells in no set.
#' @return Factor aligned to \code{cells}. Cells in several sets take the last.
#' @export
sparks_categories_to_factor <- function(categories, cells, unassigned = "Unassigned") {
  out <- rep(unassigned, length(cells))
  names(out) <- cells
  for (nm in names(categories)) out[names(out) %in% categories[[nm]]] <- nm
  lv <- unique(c(names(categories), unassigned))
  factor(unname(out), levels = lv[lv %in% out])
}

#' Read cell barcodes from a CSV / TXT file
#'
#' Accepts a file with a \code{barcode} (or \code{cell}) column, or a plain
#' one-column list -- including the CSV downloaded from the HTML report.
#' @param path File path.
#' @return Character vector of barcodes.
#' @export
sparks_read_barcodes <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  if (!length(lines)) return(character())
  first <- strsplit(lines[[1L]], "[,\t]")[[1L]]
  header <- tolower(gsub('"', "", first))
  col <- match(TRUE, header %in% c("barcode", "cell", "cells", "cell_id", "barcodes"))
  if (!is.na(col)) {
    vals <- vapply(strsplit(lines[-1L], "[,\t]"), function(p) if (length(p) >= col) p[[col]] else NA_character_, character(1))
  } else {
    vals <- vapply(strsplit(lines, "[,\t]"), `[[`, character(1), 1L)
  }
  unique(gsub('^"|"$', "", trimws(vals[!is.na(vals)])))
}
