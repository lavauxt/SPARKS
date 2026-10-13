test_that("celldex reference load failure falls back to the other reference", {
  counts <- matrix(
    1:24, nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:12))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object[["SCT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    data = Matrix::Matrix(counts, sparse = TRUE)
  )

  reference <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = Matrix::Matrix(matrix(
      c(1, 0, 0, 1), nrow = 2,
      dimnames = list(c("g1", "g2"), c("ref1", "ref2"))
    ), sparse = TRUE)),
    colData = S4Vectors::DataFrame(label = c("TypeA", "TypeB"))
  )
  loader <- function(name) {
    if (identical(name, "broken_reference")) {
      stop("HDF5. File accessibility. Unable to open file.")
    }
    reference
  }
  cfg <- list(
    unassigned_label = "Unassigned",
    min_cells_per_group = 1L,
    labels = list(list(name = "cell_type", ref_field = "label"))
  )

  expect_message(
    annotated <- SPARKS::run_singler_annotation(
      object, "Mouse", "broken_reference", "working_reference", cfg,
      reference_loader = loader
    ),
    "HDF5. File accessibility"
  )
  expect_true("cell_type" %in% colnames(annotated@meta.data))
  expect_true(all(annotated$cell_type %in% c("TypeA", "TypeB", "Unassigned")))
  expect_false(all(annotated$cell_type == "Unassigned"))
})

test_that("corrupt MouseRNAseqData cache triggers a forced refresh", {
  reference <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = Matrix::Matrix(matrix(
      c(1, 0, 0, 1), nrow = 2,
      dimnames = list(c("g1", "g2"), c("ref1", "ref2"))
    ), sparse = TRUE)),
    colData = S4Vectors::DataFrame(label.main = c("TypeA", "TypeB"))
  )
  refresh_called <- FALSE

  loaded <- SPARKS:::.safe_load_celldex_reference(
    "MouseRNAseqData",
    loader = function(name) stop("HDF5. File accessibility. Unable to open file."),
    refresh_loader = function(name) {
      refresh_called <<- TRUE
      reference
    },
    legacy_loader = function(name) stop("legacy path should not be needed")
  )

  expect_true(refresh_called)
  expect_s4_class(loaded, "SummarizedExperiment")
  expect_identical(dim(loaded), dim(reference))
})

test_that("annotation continues with unassigned labels if both references fail", {
  counts <- matrix(
    1:24, nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("cell", 1:12))
  )
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts))
  object[["SCT"]] <- SeuratObject::CreateAssay5Object(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    data = Matrix::Matrix(counts, sparse = TRUE)
  )
  loader <- function(name) stop("HDF5. File accessibility. Unable to open file.")
  cfg <- list(
    unassigned_label = "Unassigned",
    min_cells_per_group = 1L,
    labels = list(list(name = "cell_type", ref_field = "label"))
  )

  expect_message(
    annotated <- SPARKS::run_singler_annotation(
      object, "Mouse", "broken_primary", "broken_secondary", cfg,
      reference_loader = loader
    ),
    "HDF5. File accessibility"
  )
  expect_identical(unique(annotated$cell_type), "Unassigned")
})
