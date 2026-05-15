run_sex_scoring <- function(
    seurat_obj,
    config_block,
    assay = "RNA",
    ctrl = NULL,
    set.ident = FALSE
) {

    message("--- Running Sex Scoring ---")

    # ------------------------------------------------------------------
    # Check assay
    # ------------------------------------------------------------------
    if (!assay %in% names(seurat_obj@assays)) {
        message("   [ERROR] Assay '", assay, "' not found.")
        return(seurat_obj)
    }

    # ------------------------------------------------------------------
    # Get features from assay
    # ------------------------------------------------------------------
    all_genes <- SeuratObject::Features(seurat_obj[[assay]])

    female.features <- intersect(
        config_block$markers$female,
        all_genes
    )

    male.features <- intersect(
        config_block$markers$male,
        all_genes
    )

    # ------------------------------------------------------------------
    # Safety check
    # ------------------------------------------------------------------
    if (length(female.features) == 0 &&
        length(male.features) == 0) {

        message("   [WARNING] No sex markers found.")
        return(seurat_obj)
    }

    message(
        "   -> Found ",
        length(female.features),
        " female and ",
        length(male.features),
        " male markers."
    )

    # ------------------------------------------------------------------
    # Control genes
    # ------------------------------------------------------------------
    features <- list(
        female.features,
        male.features
    )

    if (is.null(ctrl)) {
        ctrl <- min(
            vapply(
                X = features,
                FUN = length,
                FUN.VALUE = numeric(1)
            )
        )
    }

    # ------------------------------------------------------------------
    # Module scoring
    # ------------------------------------------------------------------
    object.sex <- Seurat::AddModuleScore(
        object = seurat_obj,
        features = features,
        name = "Sex",
        ctrl = ctrl,
        assay = assay
    )

    # ------------------------------------------------------------------
    # Retrieve scores
    # ------------------------------------------------------------------
    sex.columns <- grep(
        pattern = "^Sex",
        x = colnames(object.sex[[]]),
        value = TRUE
    )

    sex.scores <- object.sex[[sex.columns]]

    # Rename columns
    colnames(sex.scores) <- c(
        "Female.Score",
        "Male.Score"
    )

    # ------------------------------------------------------------------
    # Sex difference
    # ------------------------------------------------------------------
    sex.scores$Sex.Difference <-
        sex.scores$Female.Score -
        sex.scores$Male.Score

    # ------------------------------------------------------------------
    # Assign sex
    # ------------------------------------------------------------------
    sex.scores$Assigned_Sex <- apply(
        X = sex.scores[, c("Female.Score", "Male.Score")],
        MARGIN = 1,
        FUN = function(scores) {

            if (all(scores < 0)) {
                return("Undetermined")
            }

            if (length(which(scores == max(scores))) > 1) {
                return("Undecided")
            }

            if (scores[1] > scores[2]) {
                return("Female")
            }

            return("Male")
        }
    )

    # ------------------------------------------------------------------
    # Add metadata
    # ------------------------------------------------------------------
    seurat_obj[[colnames(sex.scores)]] <- sex.scores

    # ------------------------------------------------------------------
    # Optional identity assignment
    # ------------------------------------------------------------------
    if (set.ident) {

        seurat_obj[["old.ident"]] <- Idents(seurat_obj)

        Idents(seurat_obj) <- "Assigned_Sex"
    }

    message(
        "   -> Added metadata: ",
        paste(colnames(sex.scores), collapse = ", ")
    )

    return(seurat_obj)
}