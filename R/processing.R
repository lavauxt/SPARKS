#' Load count matrix from Alevin or Cell Ranger output
#' @param folder_id Character
#' @param data_path Character. Root data directory
#' @return Raw count matrix
#' @keywords internal
.load_counts <- function(folder_id, data_path) {
  alevin_path <- file.path(data_path, folder_id, "alevin", "quants_mat.gz")
  tenx_path   <- file.path(data_path, folder_id, "outs",
                            "filtered_feature_bc_matrix")

  if (file.exists(alevin_path)) {
    message("   Input: Alevin  ", folder_id)
    txi  <- tximport::tximport(files = alevin_path, type = "alevin")
    cnts <- txi$counts
  } else if (dir.exists(tenx_path)) {
    message("   Input: Cell Ranger  ", folder_id)
    cnts <- Seurat::Read10X(data.dir = tenx_path)
  } else {
    stop("No valid input found in: ", file.path(data_path, folder_id),
         "\n  Tried Alevin: ",  alevin_path,
         "\n  Tried 10X:    ",  tenx_path)
  }

  if (is.list(cnts)) cnts <- cnts[["Gene Expression"]]
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
#' @param min_features Integer
#' @param max_counts Integer
#' @param max_mt_percent Numeric
#' @param min_cells Integer
#' @param cfg Named list. Full pipeline config
#' @return Seurat object or NULL
#' @export
process_single_sample <- function(folder_id, protocol, file_prefix, qc_dir,
                                  data_path, gene_removal_pattern, mt_pattern,
                                  genes_to_remove, min_features, max_counts,
                                  max_mt_percent, min_cells, cfg) {
  
  message("\n--- Processing Sample: ", file_prefix, " ---")
  make_dir(qc_dir)

   cnts <- safe_run(.load_counts(folder_id, data_path),
                   label = paste0("load_counts: ", folder_id))
  
  if (is.null(cnts)) {
    message(" [!] Error: Could not load counts for ", folder_id)
    return(NULL)
  }

  manual_list <- genes_to_remove %||% c()
  pattern_hits <- rownames(cnts)[grep(gene_removal_pattern, rownames(cnts))]
  to_remove    <- unique(c(pattern_hits, manual_list))
  present_to_remove <- intersect(to_remove, rownames(cnts))
  
  if (length(present_to_remove) > 0) {
    cnts <- cnts[!rownames(cnts) %in% present_to_remove, , drop = FALSE]
    message(" -> Filtered ", length(present_to_remove), " genes (pattern + manual list).")
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
  so[["percent.mt"]] <- Seurat::PercentageFeatureSet(so, pattern = mt_pattern)

  message(" -> Generating pre-filtering QC plots...")
  generate_qc_plots(so, qc_dir, paste0("preFiltering_", file_prefix), cfg)

  save_cell_counts(so, paste0("before_filtering_", file_prefix), qc_dir)

  # QC filtering
  so <- subset(so,
    subset = nFeature_RNA > min_features &
             nCount_RNA   < max_counts   &
             percent.mt   < max_mt_percent)

  if (ncol(so) == 0L) {
    message("   [SKIP] No cells after QC filtering: ", file_prefix)
    return(NULL)
  }

  Seurat::DefaultAssay(so) <- "RNA"
  so <- SeuratObject::JoinLayers(so)

  # Doublet removal 
  n_before <- ncol(so)
  so <- safe_run({
    sce <- scDblFinder::scDblFinder(
      Seurat::as.SingleCellExperiment(so),  
      samples = "orig.ident"
    )
    so$scDblFinder.class <- sce$scDblFinder.class
    writeLines(
      capture.output(print(table(so$scDblFinder.class))),
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

#' SCTransform + PCA + UMAP + clustering (Seurat v5)
#' @param seurat_obj Seurat object
#' @param dims_pca Integer vector
#' @param resolution Numeric
#' @param npcs Integer
#' @param vars_to_regress Character vector of metadata columns to regress
#' @return Processed Seurat object
#' @export
run_seurat_processing <- function(seurat_obj, 
                                  dims_pca = 1:20,
                                  resolution = 0.5,
                                  npcs       = 50L,
                                  vars_to_regress = "percent.mt") { 
  
  if (ncol(seurat_obj) < 10L) stop("Too few cells: ", ncol(seurat_obj))

  Seurat::DefaultAssay(seurat_obj) <- "RNA"
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj)

  message("   [SCTransform] Regressing variables: ", paste(vars_to_regress, collapse = ", "))

  seurat_obj <- Seurat::SCTransform(
    seurat_obj,
    assay           = "RNA",
    new.assay.name  = "SCT",
    vars.to.regress = vars_to_regress,
    vst.flavor      = "v2",
    verbose         = FALSE
  )

  Seurat::DefaultAssay(seurat_obj) <- "SCT"
  actual_npcs <- min(npcs, ncol(seurat_obj) - 1L)
  seurat_obj  <- Seurat::RunPCA(seurat_obj, npcs = actual_npcs, verbose = FALSE)

  actual_dims <- dims_pca[dims_pca <= actual_npcs]
  if (length(actual_dims) < 2L) {
    message("   [WARNING] dims_pca exceeds available PCs : using 1:", actual_npcs)
    actual_dims <- seq_len(actual_npcs)
  }

  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = "pca", dims = actual_dims, verbose = FALSE)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = "pca", dims = actual_dims, verbose = FALSE)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = resolution, verbose = FALSE)
  
  seurat_obj
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
