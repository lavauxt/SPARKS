#' Detect the quantification format of a sample folder
#'
#' Looks for the on-disk layout of each supported format, in order:
#' AnnData (.h5ad) > Alevin > Cell Ranger MEX > Cell Ranger HDF5 > STARsolo >
#' a bare MatrixMarket triplet.
#' @param sample_dir Character. Path to one sample's folder.
#' @return Character format code, or NA if nothing recognizable was found.
#' @keywords internal
.detect_quant_format <- function(sample_dir) {
  if (length(list.files(sample_dir, pattern = "\\.h5ad$", full.names = TRUE)) > 0L) {
    return("h5ad")
  }
  if (file.exists(file.path(sample_dir, "alevin", "quants_mat.gz"))) {
    return("alevin")
  }
  if (dir.exists(file.path(sample_dir, "outs", "filtered_feature_bc_matrix"))) {
    return("10x")
  }
  if (length(list.files(file.path(sample_dir, "outs"), pattern = "\\.h5$", full.names = TRUE)) > 0L ||
      length(list.files(sample_dir, pattern = "\\.h5$", full.names = TRUE)) > 0L) {
    return("10x_h5")
  }
  if (dir.exists(file.path(sample_dir, "Solo.out"))) {
    return("starsolo")
  }
  if (any(file.exists(file.path(sample_dir, c("matrix.mtx", "matrix.mtx.gz"))))) {
    return("mtx")
  }
  NA_character_
}

#' Read raw counts from an AnnData .h5ad file
#'
#' Tries \code{anndataR} first (pure R, no Python needed); falls back to
#' \code{zellkonverter} (Bioconductor, uses a managed Python env) if that
#' isn't installed. AnnData stores raw (pre-normalisation) counts in
#' different places depending on the pipeline that wrote the file, so
#' \code{raw_slot} controls where to look.
#' @param path Character. Path to the .h5ad file.
#' @param raw_slot Character. "auto" tries layers$counts, then raw.X, then X,
#'   in that order; or force one of "layers:counts", "raw", "X" explicitly.
#' @return Sparse matrix, genes x cells.
#' @keywords internal
.load_h5ad_counts <- function(path, raw_slot = "auto") {
  if (requireNamespace("anndataR", quietly = TRUE)) {
    adata <- anndataR::read_h5ad(path)
    search_order <- if (identical(raw_slot, "auto")) c("layers:counts", "raw", "X") else raw_slot
    mat <- NULL
    used <- NULL
    for (src in search_order) {
      mat <- switch(src,
        "layers:counts" = if (!is.null(adata$layers) && "counts" %in% names(adata$layers))
                             adata$layers[["counts"]] else NULL,
        "raw"           = if (!is.null(adata$raw)) adata$raw$X else NULL,
        "X"             = adata$X,
        NULL
      )
      if (!is.null(mat)) { used <- src; break }
    }
    if (is.null(mat)) {
      stop("Could not find a counts matrix in ", path,
           " (checked layers$counts, raw.X, X). Set input$h5ad_raw_slot explicitly.")
    }
    message("      -> using AnnData slot '", used, "' as raw counts")

    # AnnData is cells x genes (obs x var); Seurat wants genes x cells.
    mat <- Matrix::t(methods::as(mat, "CsparseMatrix"))
    rownames(mat) <- adata$var_names
    colnames(mat) <- adata$obs_names
    return(mat)
  }

  if (requireNamespace("zellkonverter", quietly = TRUE)) {
    sce <- zellkonverter::readH5AD(path)
    assay_name <- if (identical(raw_slot, "auto")) {
      an <- SummarizedExperiment::assayNames(sce)
      if ("counts" %in% an) "counts" else an[1L]
    } else raw_slot
    message("      -> using zellkonverter assay '", assay_name, "' as raw counts")
    return(methods::as(SummarizedExperiment::assay(sce, assay_name), "CsparseMatrix"))
  }

  stop(
    "Reading .h5ad files requires either the 'anndataR' or 'zellkonverter' package.\n",
    "  Install one of:\n",
    "    remotes::install_github('scverse/anndataR')\n",
    "    BiocManager::install('zellkonverter')"
  )
}

#' Read a generic MatrixMarket triplet that doesn't follow 10x file naming
#' @keywords internal
.load_generic_mtx <- function(sample_dir, matrix_file = "matrix.mtx",
                               barcodes_file = "barcodes.tsv",
                               features_file = "features.tsv") {
  .find <- function(name) {
    hits <- file.path(sample_dir, c(name, paste0(name, ".gz")))
    hits <- hits[file.exists(hits)]
    if (length(hits) == 0L) stop("Could not find ", name, " (or .gz) in ", sample_dir)
    hits[[1L]]
  }
  mtx      <- methods::as(Matrix::readMM(.find(matrix_file)), "CsparseMatrix")
  barcodes <- readLines(.find(barcodes_file))
  features <- utils::read.delim(.find(features_file), header = FALSE, stringsAsFactors = FALSE)

  if (nrow(mtx) == length(barcodes) && ncol(mtx) == nrow(features)) {
    mtx <- Matrix::t(mtx)  # was written cells x genes; transpose to genes x cells
  }
  rownames(mtx) <- features[[1L]]
  colnames(mtx) <- barcodes
  mtx
}

#' Warn if a "raw counts" matrix looks like it's actually already normalized
#' @keywords internal
.warn_if_noninteger_counts <- function(cnts, label = "") {
  vals <- if (methods::is(cnts, "sparseMatrix")) cnts@x else as.numeric(as.matrix(cnts))
  if (length(vals) == 0L) return(invisible(NULL))
  samp <- vals[sample.int(length(vals), size = min(length(vals), 20000L))]
  frac <- mean(abs(samp - round(samp)) > 1e-6)
  if (frac > 0.01) {
    warning(sprintf(
      "%s: %.1f%% of a sample of matrix values are non-integer. This usually means the matrix is already normalised/log-transformed rather than raw counts -- double check input$h5ad_raw_slot (for AnnData input) or the detected format.",
      label, 100 * frac))
  }
  invisible(NULL)
}

#' Export a raw count matrix as a standard, tool-agnostic scRNA-seq matrix
#'
#' Regardless of whether the sample was loaded from Alevin, Cell Ranger,
#' STARsolo, AnnData, or a generic MTX triplet, this writes a canonical
#' 10x-style artifact (readable by \code{Seurat::Read10X()}, Scanpy's
#' \code{read_10x_mtx()}, etc.) so every run has a portable raw-counts
#' snapshot, independent of the input format.
#' @param cnts Matrix or sparse Matrix, genes x cells.
#' @param out_dir Character. Directory to write into.
#' @param format Character. "mtx" (MatrixMarket triplet, default) or "h5"
#'   (10x-style HDF5). Both need the Bioconductor 'DropletUtils' package.
#' @return Invisible out_dir, or NULL if DropletUtils isn't installed.
#' @export
export_raw_matrix <- function(cnts, out_dir, format = "mtx") {
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    message("   [SKIP] export_raw_matrix: 'DropletUtils' not installed. ",
            "Install with: BiocManager::install('DropletUtils')")
    return(invisible(NULL))
  }
  make_dir(out_dir)
  cnts <- methods::as(cnts, "CsparseMatrix")

  # write10xCounts() needs a *directory* for the MTX triplet but a *file* path
  # for HDF5. Passing the (already created) directory with type = "HDF5" made it
  # delete the directory and leave an extension-less file named after the
  # sample in its place, instead of <out_dir>/raw_feature_bc_matrix.h5.
  is_h5    <- identical(format, "h5")
  out_path <- if (is_h5) file.path(out_dir, "raw_feature_bc_matrix.h5") else out_dir

  ok <- safe_run({
    DropletUtils::write10xCounts(
      path        = out_path,
      x           = cnts,
      barcodes    = colnames(cnts),
      gene.id     = rownames(cnts),
      gene.symbol = rownames(cnts),
      type        = if (is_h5) "HDF5" else "sparse",
      version     = "3",
      overwrite   = TRUE
    )
    TRUE
  }, label = paste0("export_raw_matrix -> ", out_dir), fallback = FALSE)

  if (isTRUE(ok)) {
    message("   Raw matrix (", nrow(cnts), " genes x ", ncol(cnts),
            " cells) saved to: ", out_path, " [", format, "]")
  }
  invisible(out_dir)
}

#' Load count matrix from Alevin, Cell Ranger (MEX or HDF5), STARsolo,
#' AnnData (.h5ad), or a generic MatrixMarket triplet
#' @param folder_id Character
#' @param data_path Character. Root data directory
#' @param format Character or NULL. "auto" (default when NULL) inspects the
#'   folder; or force one of "alevin", "10x", "10x_h5", "starsolo", "h5ad", "mtx".
#' @param cfg Named list or NULL. Full pipeline config, used for
#'   \code{input$format} / \code{input$starsolo_feature} / \code{input$h5ad_raw_slot}.
#' @return Raw count matrix (genes x cells), tagged with a "quant_format" attribute
#' @keywords internal
.load_counts <- function(folder_id, data_path, format = NULL, cfg = NULL) {
  sample_dir  <- file.path(data_path, folder_id)
  alevin_path <- file.path(sample_dir, "alevin", "quants_mat.gz")
  tenx_path   <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix")

  format <- format %||% cfg$input$format %||% "auto"
  if (identical(format, "auto")) {
    format <- .detect_quant_format(sample_dir)
    if (is.na(format)) {
      stop("No valid input found in: ", sample_dir,
           "\n  Tried Alevin:      ", alevin_path,
           "\n  Tried Cell Ranger: ", tenx_path,
           "\n  Tried STARsolo:    ", file.path(sample_dir, "Solo.out"),
           "\n  Tried .h5ad / .h5 / bare matrix.mtx directly under: ", sample_dir,
           "\n  Set `input$format` in the config, or a `quant_format` column ",
           "in the sample table, to specify the format explicitly.")
    }
  }

  cnts <- switch(format,
    "alevin" = {
      message("   Input: Alevin  ", folder_id)
      tximport::tximport(files = alevin_path, type = "alevin")$counts
    },
    "10x" = {
      message("   Input: Cell Ranger (MEX)  ", folder_id)
      Seurat::Read10X(data.dir = tenx_path)
    },
    "10x_h5" = {
      h5_candidates <- c(
        list.files(file.path(sample_dir, "outs"), pattern = "\\.h5$", full.names = TRUE),
        list.files(sample_dir, pattern = "\\.h5$", full.names = TRUE)
      )
      if (length(h5_candidates) == 0L) stop("No .h5 file found for ", folder_id, " in ", sample_dir)
      message("   Input: Cell Ranger (HDF5)  ", folder_id, " <- ", basename(h5_candidates[1L]))
      Seurat::Read10X_h5(h5_candidates[1L])
    },
    "starsolo" = {
      feature    <- cfg$input$starsolo_feature %||% "Gene"
      candidates <- file.path(sample_dir, "Solo.out", feature, c("filtered", "raw"))
      d          <- candidates[dir.exists(candidates)]
      if (length(d) == 0L) stop("No Solo.out/", feature, "/{filtered,raw} found for ", folder_id)
      message("   Input: STARsolo (", feature, ")  ", folder_id, " <- ", d[1L])
      Seurat::Read10X(data.dir = d[1L], gene.column = 2L)
    },
    "h5ad" = {
      h5ad_files <- list.files(sample_dir, pattern = "\\.h5ad$", full.names = TRUE)
      if (length(h5ad_files) == 0L) stop("No .h5ad file found for ", folder_id, " in ", sample_dir)
      message("   Input: AnnData (.h5ad)  ", folder_id, " <- ", basename(h5ad_files[1L]))
      .load_h5ad_counts(h5ad_files[1L], raw_slot = cfg$input$h5ad_raw_slot %||% "auto")
    },
    "mtx" = {
      message("   Input: generic MatrixMarket triplet  ", folder_id)
      .load_generic_mtx(sample_dir)
    },
    stop("Unknown quant_format '", format, "' for sample '", folder_id,
         "'. Supported: auto, alevin, 10x, 10x_h5, starsolo, h5ad, mtx.")
  )

  if (is.list(cnts)) cnts <- cnts[["Gene Expression"]] %||% cnts[[1L]]
  attr(cnts, "quant_format") <- format
  cnts
}

#' Process a single sample: load, filter, doublet removal
#' @param folder_id Character
#' @param protocol Character. Condition label (e.g. "WT")
#' @param file_prefix Character
#' @param qc_dir Character
#' @param data_path Character
#' @param gene_removal_pattern Character regex
#' @param mt_pattern Character regex
#' @param genes_to_remove Character vector
#' @param min_features Integer. Cells need more than this many detected genes
#' @param max_features Integer. Cells need fewer than this many detected genes
#'   (config default 6000 when \code{qc$max_features} is not set)
#' @param max_counts Integer
#' @param max_mt_percent Numeric
#' @param min_cells Integer
#' @param cfg Named list. Full pipeline config
#' @param raw_matrix_dir Character or NULL. If given and
#'   \code{cfg$input$save_raw_matrix} is TRUE, the as-loaded (pre-QC) counts
#'   matrix is exported here in standard 10x-style format via
#'   \code{export_raw_matrix()}.
#' @param quant_format Character or NULL. Per-sample override for the input
#'   format (see \code{.load_counts()}); NULL defers to \code{cfg$input$format}.
#' @return Seurat object or NULL
#' @export
process_single_sample <- function(folder_id, protocol, file_prefix, qc_dir,
                                  data_path, gene_removal_pattern, mt_pattern,
                                  genes_to_remove, min_features, max_features,
                                  max_counts, max_mt_percent, min_cells, cfg,
                                  raw_matrix_dir = NULL, quant_format = NULL) {
  
  message("\n--- Processing Sample: ", file_prefix, " ---")
  make_dir(qc_dir)

  cnts <- safe_run(.load_counts(folder_id, data_path, format = quant_format, cfg = cfg),
                   label = paste0("load_counts: ", folder_id))
  
  if (is.null(cnts)) {
    message(" [!] Error: Could not load counts for ", folder_id)
    return(NULL)
  }

  safe_run(.warn_if_noninteger_counts(cnts, folder_id), label = "integer-counts check")

  if (isTRUE(cfg$input$save_raw_matrix) && !is.null(raw_matrix_dir)) {
    export_raw_matrix(cnts, file.path(raw_matrix_dir, folder_id),
                       format = cfg$input$raw_matrix_format %||% "mtx")
  }

  message(" -> Creating Seurat object...")
  so <- Seurat::CreateSeuratObject(
    counts       = cnts,
    project      = file_prefix,
    min.cells    = min_cells,
    min.features = min_features
  )

  so$orig.ident      <- file_prefix
  so$condition       <- protocol
  so$quant_format    <- attr(cnts, "quant_format") %||% quant_format %||% "auto"
  
  so[["percent.mt"]] <- Seurat::PercentageFeatureSet(so, pattern = mt_pattern)

  message(" -> Generating pre-filtering QC plots...")
  generate_qc_plots(so, qc_dir, paste0("preFiltering_", file_prefix), cfg)

  save_cell_counts(so, paste0("before_filtering_", file_prefix), qc_dir)

  so <- subset(so,
    subset = nFeature_RNA > min_features &
             nFeature_RNA < max_features  &
             nCount_RNA   < max_counts   &
             percent.mt   < max_mt_percent)

  if (ncol(so) == 0L) {
    message("   [SKIP] No cells after QC filtering: ", file_prefix)
    return(NULL)
  }

  manual_list <- genes_to_remove %||% c()
  # An empty/NA pattern means "no pattern-based removal": grep("", x) matches
  # *every* gene and would wipe the whole matrix.
  pattern_hits <- if (length(gene_removal_pattern) == 1L &&
                      !is.na(gene_removal_pattern) && nzchar(gene_removal_pattern)) {
    rownames(so)[grep(gene_removal_pattern, rownames(so))]
  } else {
    character(0)
  }
  to_remove    <- unique(c(pattern_hits, manual_list))
  present_to_remove <- intersect(to_remove, rownames(so))
  
  if (length(present_to_remove) > 0) {
    genes_to_keep <- setdiff(rownames(so), present_to_remove)
    so <- subset(so, features = genes_to_keep)
    message(" -> Filtered ", length(present_to_remove), " genes (pattern + manual list) after QC.")
  }

  Seurat::DefaultAssay(so) <- "RNA"
  so <- SeuratObject::JoinLayers(so)

  n_before <- ncol(so)
  so <- safe_run({
    sce <- scDblFinder::scDblFinder(
      Seurat::as.SingleCellExperiment(so),  
      samples = "orig.ident"
    )
    so$scDblFinder.class <- sce$scDblFinder.class
    writeLines(
      utils::capture.output(print(table(so$scDblFinder.class))),
      file.path(qc_dir, paste0("DoubletStats_", file_prefix, ".txt"))
    )
    sub <- subset(so, subset = scDblFinder.class == "singlet")

    Seurat::DefaultAssay(sub) <- "RNA"
    sub <- SeuratObject::JoinLayers(sub)
    sub
  }, label = "scDblFinder", fallback = {
    so$scDblFinder.class <- "singlet"
    so
  })

  n_after <- ncol(so)
  writeLines(
    paste("Removed by scDblFinder:", n_before - n_after),
    file.path(qc_dir, paste0("Doublet_Removed_", file_prefix, ".txt"))
  )

  generate_qc_plots(so, qc_dir, paste0("postFiltering_", file_prefix), cfg)
  save_cell_counts(so,           paste0("after_filtering_", file_prefix), qc_dir)
  message("   Final cells count: ", ncol(so))
  so
}

#' SCTransform + PCA + Integration + UMAP + clustering (Seurat v5 Pro Version)
#' @param seurat_obj Seurat object
#' @param dims_pca Integer vector
#' @param resolution Numeric
#' @param npcs Integer
#' @param vars_to_regress Character vector of metadata columns to regress
#' @param split_by Character. Metadata column to split layers for batch integration and Harmony grouping
#' @param prep_sct_findmarkers Logical. Run Seurat's \code{PrepSCTFindMarkers()}
#'   right after SCTransform so the SCT counts/data of all samples are on a common
#'   depth (needed for between-sample DE). It has to happen here: the Assay5
#'   conversion below discards the per-sample SCT models it needs. No-op with a
#'   single sample.
#' @return Processed Seurat object
#' @export
run_seurat_processing <- function(seurat_obj, 
                                  dims_pca = 1:20,
                                  resolution = 0.5,
                                  npcs       = 50L,
                                  vars_to_regress = "percent.mt",
                                  split_by   = "orig.ident",
                                  prep_sct_findmarkers = TRUE) { 
  
  if (ncol(seurat_obj) < 10L) stop("Too few cells: ", ncol(seurat_obj))

  Seurat::DefaultAssay(seurat_obj) <- "RNA"

  message("   [Integration] Splitting RNA assay by: ", split_by)
  seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj@meta.data[[split_by]])

  message("   [SCTransform] Regressing variables: ", paste(vars_to_regress, collapse = ", "))
  seurat_obj <- Seurat::SCTransform(
    seurat_obj,
    assay           = "RNA",
    new.assay.name  = "SCT",
    vars.to.regress = vars_to_regress,
    vst.flavor      = "v2",
    verbose         = FALSE
  )

  # SCTransform() returns an SCTAssay holding one SCT model per layer (here: per
  # sample). Cross-sample DE / average expression needs the corrected counts and
  # data re-expressed at a common sequencing depth, which PrepSCTFindMarkers()
  # computes *from those models*. The Assay5 coercion below drops them (Assay5
  # has no SCTModel.list) and FindMarkers() does no depth check on an Assay5, so
  # Prep has to run here -- afterwards it can no longer run at all.
  # (run_analysis_unit() used to skip it as "not required" on every run.)
  # No-op for a single model; a failure is reported but not fatal.
  if (isTRUE(prep_sct_findmarkers) && inherits(seurat_obj[["SCT"]], "SCTAssay")) {
    message("   [SCTransform] Running PrepSCTFindMarkers (", length(levels(seurat_obj[["SCT"]])),
            " SCT model(s))...")
    seurat_obj <- safe_run(
      Seurat::PrepSCTFindMarkers(seurat_obj, assay = "SCT", verbose = FALSE),
      label    = "PrepSCTFindMarkers",
      fallback = seurat_obj
    )
  }

  if (!inherits(seurat_obj[["SCT"]], "Assay5")) {
    seurat_obj[["SCT"]] <- methods::as(seurat_obj[["SCT"]], "Assay5")
  }

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  actual_npcs <- min(npcs, ncol(seurat_obj) - 1L)
  seurat_obj  <- Seurat::RunPCA(seurat_obj, npcs = actual_npcs, verbose = FALSE)

  actual_dims <- dims_pca[dims_pca <= actual_npcs]
  if (length(actual_dims) < 2L) {
    message("   [WARNING] dims_pca exceeds available PCs : using 1:", actual_npcs)
    actual_dims <- seq_len(actual_npcs)
  }

  message("   [Integration] Rejoining layers before integration...")
  seurat_obj[["RNA"]] <- SeuratObject::JoinLayers(seurat_obj[["RNA"]])
  seurat_obj[["SCT"]] <- SeuratObject::JoinLayers(seurat_obj[["SCT"]])

  n_batches <- length(unique(stats::na.omit(as.character(seurat_obj@meta.data[[split_by]]))))
  if (n_batches >= 2L) {
    message("   [Integration] Running direct Harmony integration (", n_batches, " batches)...")
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("The 'harmony' package is missing. Please run: install.packages('harmony')")
    }

    seurat_obj <- harmony::RunHarmony(
      object         = seurat_obj,
      group.by.vars  = split_by,
      reduction.use  = "pca",
      reduction.save = "integrated.dr",
      verbose        = FALSE
    )
    cluster_red <- "integrated.dr"
  } else {
    message("   [Integration] Single '", split_by, "' value: skipping Harmony, using PCA.")
    cluster_red <- "pca"
  }

  message("   [Clustering] Running UMAP and FindNeighbors on ", cluster_red, "...")
  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = cluster_red, dims = actual_dims, verbose = FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = cluster_red, dims = actual_dims, verbose = FALSE)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  
  return(seurat_obj)
}

#' Assign cell-type labels using positive/negative marker gene rules
#'
#' Rules are applied in order; later rules overwrite earlier ones on the same cells.
#'
#' @param seurat_obj Seurat object
#' @param new_col_name Character. Metadata column to create
#' @param rules List of lists, each with:
#'   \itemize{
#'     \item \code{label}    Character. Label to assign
#'     \item \code{positive} Character vector. Genes that must be expressed
#'     \item \code{negative} Character vector. Genes that must NOT be expressed
#'   }
#' @param unassigned_label Character
#' @param threshold Numeric between 0 and 1. Min fraction of positive genes expressed to qualify
#' @return Seurat object with new metadata column
#' @export
label_cells_by_markers <- function(seurat_obj, new_col_name, rules,
                                    unassigned_label = "Unassigned",
                                    threshold = 0.1) {

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  expr   <- SeuratObject::LayerData(seurat_obj, assay = "SCT", layer = "data")
  labels <- rep(unassigned_label, ncol(seurat_obj))

  for (rule in rules) {
    lbl      <- rule$label
    pos_req  <- as.character(rule$positive %||% character())
    neg_req  <- as.character(rule$negative %||% character())
    pos_genes <- intersect(pos_req, rownames(expr))
    neg_genes <- intersect(neg_req, rownames(expr))
    if (length(pos_req) > 0L && length(pos_genes) == 0L) next
    if (length(pos_genes) > 0L) {
      pos_score <- Matrix::colMeans(expr[pos_genes, , drop = FALSE] > 0)
    } else {
      pos_score <- rep(1, ncol(seurat_obj))   
    }

    if (length(neg_genes) > 0L) {
      neg_fail <- Matrix::colMeans(expr[neg_genes, , drop = FALSE] > 0) >= threshold
    } else {
      neg_fail <- rep(FALSE, ncol(seurat_obj))
    }

    labels[pos_score >= threshold & !neg_fail] <- lbl
  }

  seurat_obj[[new_col_name]] <- labels
  message("   [label_cells] '", new_col_name, "' distribution:")
  print(table(seurat_obj@meta.data[[new_col_name]]))
  seurat_obj
}

#' Annotate cells with SingleR (primary + secondary references)
#' @param seurat_obj Seurat object
#' @param species_target species for celldex db
#' @param ref_celldex1 celldex reference object 1
#' @param ref_celldex2 celldex reference object 2
#' @param singler_cfg Named list from config
#' @return Seurat object with annotation columns added
#' @export
run_singler_annotation <- function(seurat_obj, species_target, ref_celldex1, ref_celldex2,
                                    singler_cfg) {
  message("--- Starting SingleR annotation ---")   
  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  sr_data <- Seurat::GetAssayData(seurat_obj, assay = "SCT", layer = "data")
  
  ref_primary <- do.call(getExportedValue("celldex", ref_celldex1), args = list())
  ref_secondary <- do.call(getExportedValue("celldex", ref_celldex2), args = list())

  unassigned <- singler_cfg$unassigned_label
  min_cells  <- singler_cfg$min_cells_per_group

  for (lt in singler_cfg$labels) {
    message("   SingleR: ", lt$name)

    sr1 <- safe_run(
      SingleR::SingleR(test   = sr_data, ref = ref_primary,
                       labels = ref_primary[[lt$ref_field]]),
      label = paste0("SingleR primary (", lt$name, ")")
    )
    sr2 <- safe_run(
      SingleR::SingleR(test   = sr_data, ref = ref_secondary,
                       labels = ref_secondary[[lt$ref_field]]),
      label = paste0("SingleR secondary (", lt$name, ")")
    )

    if (is.null(sr1) && is.null(sr2)) {
      message("   [SKIP] SingleR failed for: ", lt$name)
      seurat_obj[[lt$name]] <- unassigned
      next
    }

    labels <- if (!is.null(sr1)) sr1$labels
              else                rep(NA_character_, ncol(seurat_obj))
    if (!is.null(sr2)) {
      na_idx         <- is.na(labels)
      labels[na_idx] <- sr2$labels[na_idx]
    }
    labels[is.na(labels)] <- unassigned
    seurat_obj[[lt$name]] <- labels

    counts  <- table(seurat_obj@meta.data[[lt$name]])
    valid   <- names(counts)[counts >= min_cells]
    low_idx <- !seurat_obj@meta.data[[lt$name]] %in% valid
    seurat_obj@meta.data[[lt$name]][low_idx] <- unassigned
    message("   [", lt$name, "] ", sum(low_idx), " cells reassigned to '",
            unassigned, "' (< ", min_cells, ")")
  }
  seurat_obj
}
