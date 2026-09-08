run_sex_scoring <- function(
    seurat_obj,
    config_block,
    assay = "RNA",
    ctrl = NULL,
    set.ident = FALSE
) {
    message("--- Running Sex Scoring ---")

    if (!assay %in% names(seurat_obj@assays)) {
        message("   [ERROR] Assay '", assay, "' not found.")
        return(seurat_obj)
    }

    all_genes <- SeuratObject::Features(seurat_obj[[assay]])

    female.features <- intersect(config_block$markers$female, all_genes)
    male.features <- intersect(config_block$markers$male, all_genes)

    if (length(female.features) == 0 && length(male.features) == 0) {
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

    features <- list(female.features, male.features)

    # BUG FIX: ctrl was min(length(female.features), length(male.features)) --
    # e.g. 2 for XIST/TSIX-only female markers. But `ctrl` in AddModuleScore()
    # is the size of the randomly-sampled, expression-binned *background*
    # gene pool used to correct each module score -- unrelated to how many
    # marker genes define the module itself. ctrl=2 makes that background
    # correction extremely noisy. Default to Seurat's own ctrl default (100),
    # capped by how many genes actually exist to sample from.
    if (is.null(ctrl)) {
        ctrl <- min(100L, length(all_genes) - 1L)
    }

    if (ctrl < 1L) {
        message("   [WARNING] ctrl < 1 (too few genes in assay). Cannot run AddModuleScore.")
        return(seurat_obj)
    }

    object.sex <- Seurat::AddModuleScore(
        object = seurat_obj,
        features = features,
        name = "Sex",
        ctrl = ctrl,
        assay = assay
    )

    # BUG FIX: "^Sex" also matches any pre-existing metadata column literally
    # named "Sex" (e.g. biological sex recorded in a sample sheet, common in
    # real datasets) or anything else starting with "Sex" -- silently pulling
    # in a 3rd column would make `colnames(sex.scores) <-
    # c("Female.Score","Male.Score")` below mis-assign names (or error).
    # AddModuleScore(name="Sex") always creates exactly "Sex1"/"Sex2" for a
    # 2-feature-set call, so match that precisely instead.
    sex.columns <- grep(
        pattern = "^Sex[0-9]+$",
        x = colnames(object.sex[[]]),
        value = TRUE
    )

    if (length(sex.columns) != 2L) {
        message("   [WARNING] Expected 2 AddModuleScore columns, found ",
                length(sex.columns), ". Skipping sex scoring.")
        return(seurat_obj)
    }

    sex.scores <- object.sex[[sex.columns]]

    colnames(sex.scores) <- c("Female.Score", "Male.Score")

    sex.scores$Sex.Difference <-
        sex.scores$Female.Score - sex.scores$Male.Score

    sex.scores$Assigned_Sex <- apply(
        X = sex.scores[, c("Female.Score", "Male.Score")],
        MARGIN = 1,
        FUN = function(scores) {
            if (all(scores < 0)) return("Undetermined")
            if (length(which(scores == max(scores))) > 1) return("Undecided")
            if (scores[1] > scores[2]) return("Female")
            return("Male")
        }
    )

    seurat_obj[[colnames(sex.scores)]] <- sex.scores

    if (set.ident) {
        seurat_obj[["old.ident"]] <- Seurat::Idents(seurat_obj)
        Seurat::Idents(seurat_obj) <- "Assigned_Sex"
    }

    message(
        "   -> Added metadata: ",
        paste(colnames(sex.scores), collapse = ", ")
    )

    return(seurat_obj)
}