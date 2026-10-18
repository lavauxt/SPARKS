#' Create the SPARKS targets pipeline
#'
#' Each comparison group is a separate target, allowing successful groups to
#' be reused after a failure. The pipeline uses the same `sparks()` API.
#' @param base_config_path Path to the base YAML configuration.
#' @param override_config_path Optional override YAML configuration.
#' @param results_dir Directory for targets store and group outputs.
#' @param target_names Optional target names to build; dependencies are built
#'   automatically and unchanged targets remain cached.
#' @return Invisibly, the targets store path.
#' @export
run_sparks_targets <- function(base_config_path, override_config_path = NULL,
                               results_dir = NULL, target_names = NULL) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("The targets feature requires package 'targets'. Install it with ",
         "install.packages('targets').", call. = FALSE)
  }
  if (!file.exists(base_config_path)) {
    stop("Base config not found: ", base_config_path, call. = FALSE)
  }
  if (!is.null(override_config_path) && !file.exists(override_config_path)) {
    stop("Override config not found: ", override_config_path, call. = FALSE)
  }

  cfg <- load_pipeline_config(base_config_path, override_config_path)
  results_dir <- results_dir %||% cfg$pipeline$results_dir
  if (length(results_dir) != 1L || is.na(results_dir) || !nzchar(results_dir)) {
    stop("results_dir must be one non-empty path.", call. = FALSE)
  }
  results_dir <- normalizePath(results_dir, winslash = "/", mustWork = FALSE)
  make_dir(results_dir)
  targets_dir <- file.path(results_dir, "_targets")
  make_dir(targets_dir)
  lock_dir <- file.path(targets_dir, "sparks-run.lock")
  if (!dir.create(lock_dir, showWarnings = FALSE)) {
    stop("A SPARKS targets run already holds this store lock: ", lock_dir,
         ". If no run is active, remove the stale lock directory manually.",
         call. = FALSE)
  }
  on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
  base_config_path <- normalizePath(base_config_path, winslash = "/", mustWork = TRUE)
  if (!is.null(override_config_path)) {
    override_config_path <- normalizePath(override_config_path, winslash = "/",
                                          mustWork = TRUE)
  }

  script_path <- file.path(targets_dir, "_targets.R")
  package_path <- getNamespaceInfo(asNamespace("SPARKS"), "path")
  script <- .targets_plan_script(base_config_path, override_config_path,
                                 results_dir, package_path)
  writeLines(script, script_path, useBytes = TRUE)
  if (!is.null(target_names)) {
    valid_names <- .targets_plan_names(base_config_path, override_config_path,
                                       results_dir)
    .validate_targets_selection(target_names, valid_names)
  }

  old_project <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  library_paths <- unique(c(dirname(package_path), .libPaths()))
  Sys.setenv(R_LIBS_USER = paste(library_paths,
                                 collapse = .Platform$path.sep))
  on.exit({
    if (is.na(old_project)) Sys.unsetenv("R_LIBS_USER") else
      Sys.setenv(R_LIBS_USER = old_project)
  }, add = TRUE)
  targets::tar_make(script = script_path,
                    store = file.path(targets_dir, "store"),
                    names = target_names)
  invisible(file.path(targets_dir, "store"))
}

.write_ui_targets_override <- function(input, results_dir) {
  override <- file.path(results_dir, "shiny_override.yaml")
  make_dir(dirname(override))
  cfg <- if (nzchar(input$yaml_text)) {
    yaml::yaml.load(input$yaml_text)
  } else if (nzchar(input$override_config)) {
    yaml::yaml.load_file(input$override_config) %||% list()
  } else {
    list()
  }
  if (is.null(cfg)) cfg <- list()
  if (!nzchar(input$yaml_text)) {
    cfg$pipeline <- modifyList(cfg$pipeline %||% list(),
                               list(species_target = input$species))
    cfg$qc <- modifyList(cfg$qc %||% list(),
                         list(min_features = input$min_features,
                              max_mt_percent = input$max_mt))
    cfg$processing <- modifyList(cfg$processing %||% list(),
                                 list(cluster_resolution = input$resolution))
    cfg$deg <- modifyList(cfg$deg %||% list(), list(run = input$run_deg))
    cfg$pseudobulk <- modifyList(cfg$pseudobulk %||% list(),
                                 list(run = input$run_pseudobulk))
  }
  if (!is.null(input$sample_file)) {
    sample_dest <- file.path(results_dir, "sample_table.tsv")
    file.copy(input$sample_file$datapath, sample_dest, overwrite = TRUE)
    if (is.null(cfg$pipeline)) cfg$pipeline <- list()
    cfg$pipeline$sample_table <- normalizePath(sample_dest, winslash = "/",
                                               mustWork = TRUE)
  }
  yaml::write_yaml(cfg, override)
  override
}

.targets_plan_script <- function(base_config_path, override_config_path,
                                 results_dir, package_path = NULL) {
  deparse_one <- function(x) paste(deparse(x), collapse = "")
  package_loader <- if (is.null(package_path)) "" else paste0(
    "if (!requireNamespace('SPARKS', quietly = TRUE)) ",
    "pkgload::load_all(", deparse_one(package_path), ", quiet = TRUE)\n"
  )
  paste0(
    package_loader,
    "SPARKS:::.build_sparks_targets_plan(",
    deparse_one(normalizePath(base_config_path, winslash = "/", mustWork = FALSE)),
    ", ",
    deparse_one(if (is.null(override_config_path)) NULL else
      normalizePath(override_config_path, winslash = "/", mustWork = FALSE)),
    ", ", deparse_one(results_dir), ")"
  )
}

.targets_plan_names <- function(base_config_path, override_config_path = NULL,
                                results_dir = tempdir()) {
  plan <- .build_sparks_targets_plan(base_config_path, override_config_path,
                                     results_dir)
  vapply(plan, `[[`, character(1), "name")
}

.validate_targets_selection <- function(target_names, valid_names) {
  if (!is.character(target_names) || !length(target_names) ||
      anyNA(target_names) || any(!nzchar(target_names)) ||
      anyDuplicated(target_names) || any(!target_names %in% valid_names)) {
    stop("Target selection is invalid: target_names must contain unique names from the SPARKS targets plan.",
         call. = FALSE)
  }
  invisible(target_names)
}

.write_result_download <- function(path, file, format = "Original") {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "tsv", "txt") && format %in% c("CSV", "TSV")) {
    data <- utils::read.delim(path,
      sep = if (ext == "csv") "," else "\t", check.names = FALSE,
      stringsAsFactors = FALSE)
    if (format == "CSV") {
      utils::write.csv(data, file, row.names = FALSE, na = "")
    } else {
      utils::write.table(data, file, sep = "\t", quote = FALSE,
                         row.names = FALSE, na = "")
    }
  } else if (!file.copy(path, file, overwrite = TRUE)) {
    stop("Could not copy selected result for download.", call. = FALSE)
  }
  invisible(file)
}

.log_run_event <- function(path, event, detail = "") {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), event, detail)
  cat(line, "\n", file = path, append = TRUE, sep = "")
  invisible(line)
}

.validate_ui_configuration <- function(base_config_path,
                                       override_config_path = NULL,
                                       results_dir, target_names = NULL) {
  if (length(base_config_path) != 1L || !file.exists(base_config_path)) {
    stop("Base config not found: ", base_config_path, call. = FALSE)
  }
  if (!is.null(override_config_path) &&
      (length(override_config_path) != 1L ||
       !file.exists(override_config_path))) {
    stop("Override config not found: ", override_config_path, call. = FALSE)
  }
  if (length(results_dir) != 1L || is.na(results_dir) ||
      !nzchar(trimws(results_dir))) {
    stop("Results directory must be one non-empty path.", call. = FALSE)
  }
  cfg <- load_pipeline_config(base_config_path, override_config_path)
  sample_table <- load_sample_table(cfg)
  plan_names <- .targets_plan_names(base_config_path, override_config_path,
                                    results_dir)
  if (!is.null(target_names)) {
    .validate_targets_selection(target_names, plan_names)
  }
  list(config = cfg, sample_table = sample_table, targets = plan_names)
}

.format_target_progress <- function(summary) {
  if (is.null(summary) || !is.data.frame(summary) || !nrow(summary)) {
    return("Target progress is not available yet.")
  }
  fields <- intersect(c("skipped", "dispatched", "completed", "errored",
                        "canceled"), names(summary))
  paste(vapply(fields, function(field) {
    paste0(field, ": ", sum(summary[[field]], na.rm = TRUE))
  }, character(1)), collapse = " | ")
}

.selected_cell_rows <- function(data, selected, cells = NULL) {
  if (!is.data.frame(data) || !nrow(data) || is.null(selected) ||
      !is.data.frame(selected) || !"key" %in% names(selected)) {
    if (!is.null(cells) && is.data.frame(data) && "cell" %in% names(data)) {
      return(data[data$cell %in% cells, , drop = FALSE])
    }
    if (is.data.frame(data)) return(data[0, , drop = FALSE])
    return(data.frame())
  }
  data[data$cell %in% as.character(selected$key), , drop = FALSE]
}

.overview_color_values <- function(values) {
  if (is.numeric(values)) return(as.numeric(values))
  factor(values, levels = unique(as.character(values)))
}

.format_run_elapsed <- function(started_at, now = Sys.time()) {
  if (is.null(started_at) || length(started_at) != 1L ||
      is.na(started_at) ||
      !(inherits(started_at, "POSIXt") || is.numeric(started_at))) {
    return("Not started")
  }
  seconds <- max(0, as.numeric(difftime(now, started_at, units = "secs")))
  sprintf("%02d:%02d:%02d", as.integer(seconds %/% 3600),
          as.integer((seconds %% 3600) %/% 60),
          as.integer(floor(seconds %% 60)))
}

.summarise_results <- function(catalog) {
  if (!is.data.frame(catalog) || !nrow(catalog)) {
    return(data.frame(group = character(), type = character(), files = integer(),
                      size_mb = numeric(), stringsAsFactors = FALSE))
  }
  group <- sub("[/\\\\].*$", "", catalog$name)
  group[!grepl("[/\\\\]", catalog$name)] <- "General"
  size_mb <- catalog$size / 1024^2
  summary <- stats::aggregate(
    cbind(files = rep(1L, nrow(catalog)), size_mb = size_mb),
    by = list(group = group, type = catalog$type), FUN = sum
  )
  summary[order(summary$group, summary$type), , drop = FALSE]
}

.result_overview_metrics <- function(catalog) {
  list(
    files = if (is.data.frame(catalog)) nrow(catalog) else 0L,
    groups = if (is.data.frame(catalog) && nrow(catalog)) {
      length(unique(sub("[/\\\\].*$", "", catalog$name)))
    } else 0L,
    size_mb = if (is.data.frame(catalog) && nrow(catalog)) {
      sum(catalog$size, na.rm = TRUE) / 1024^2
    } else 0
  )
}

.targets_run_progress <- function(store) {
  if (!requireNamespace("targets", quietly = TRUE) ||
      !dir.exists(file.path(store, "meta"))) {
    return("Targets have not initialized yet.")
  }
  summary <- tryCatch(targets::tar_progress_summary(store = store),
                      error = function(e) NULL)
  .format_target_progress(summary)
}

.explorer_features <- function(object, max_genes = 1000L) {
  metadata <- colnames(object[[]])
  genes <- rownames(object)
  if (length(genes) > max_genes) genes <- genes[seq_len(max_genes)]
  c(metadata, paste0("gene:", genes))
}

.explorer_feature_kind <- function(feature, metadata_names) {
  if (startsWith(feature, "gene:") &&
      !substring(feature, 6L) %in% metadata_names) "gene" else "metadata"
}

.discover_merged_objects <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    return(data.frame(group = character(), path = character(),
                      stringsAsFactors = FALSE))
  }
  paths <- list.files(results_dir, pattern = "^Merged_.*\\.rds$",
                      recursive = TRUE, full.names = TRUE)
  paths <- paths[grepl("[/\\\\]RData[/\\\\]", paths)]
  data.frame(
    group = sub("^Merged_(.*)\\.rds$", "\\1", basename(paths)),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

.explorer_feature_values <- function(object, feature) {
  metadata <- object[[]]
  if (feature %in% colnames(metadata)) {
    values <- metadata[[feature]]
    names(values) <- rownames(metadata)
    return(values)
  }
  if (startsWith(feature, "gene:") && requireNamespace("Seurat", quietly = TRUE)) {
    gene <- substring(feature, 6L)
    values <- tryCatch(
      Seurat::FetchData(object, vars = gene, layer = "data")[[gene]],
      error = function(e) NULL
    )
    if (!is.null(values) && is.null(names(values))) {
      names(values) <- rownames(metadata)
    }
    return(values)
  }
  NULL
}

.explorer_plot_data <- function(object, feature, group_by = NULL,
                                max_cells = 50000L) {
  reductions <- names(object@reductions)
  reduction <- reductions[tolower(reductions) == "umap"][1L]
  if (is.na(reduction) || !length(reduction)) {
    stop("The selected Seurat object has no UMAP reduction.", call. = FALSE)
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (ncol(embedding) < 2L) {
    stop("The selected embedding has fewer than two dimensions.", call. = FALSE)
  }
  metadata <- object[[]]
  cells <- intersect(rownames(embedding), rownames(metadata))
  values <- .explorer_feature_values(object, feature)
  if (is.null(values)) {
    stop("Selected metadata field or gene is unavailable.", call. = FALSE)
  }
  names(values) <- names(values) %||% rownames(metadata)
  cells <- intersect(cells, names(values))
  if (!length(cells)) stop("No matching cells for the selected feature.",
                           call. = FALSE)
  if (length(cells) > max_cells) {
    cells <- cells[seq_len(max_cells)]
  }
  data <- data.frame(
    cell = cells,
    x = embedding[cells, 1L],
    y = embedding[cells, 2L],
    value = .overview_color_values(values[cells]),
    stringsAsFactors = FALSE
  )
  if (!is.null(group_by) && group_by %in% colnames(metadata)) {
    data$group <- metadata[cells, group_by, drop = TRUE]
  } else {
    data$group <- "All cells"
  }
  data
}

.gene_expression_frame <- function(object, gene, group_by = NULL,
                                   max_cells = 50000L) {
  metadata <- object[[]]
  expression <- tryCatch(
    Seurat::FetchData(object, vars = gene, layer = "data")[[gene]],
    error = function(e) Seurat::FetchData(object, vars = gene,
                                           layer = "counts")[[gene]]
  )
  cells <- intersect(names(expression), rownames(metadata))
  if (!length(cells)) {
    stop("No expression values found for the selected gene.", call. = FALSE)
  }
  if (length(cells) > max_cells) cells <- cells[seq_len(max_cells)]
  data <- data.frame(cell = cells, expression = expression[cells],
                     stringsAsFactors = FALSE)
  if (!is.null(group_by) && group_by %in% colnames(metadata)) {
    data$group <- metadata[cells, group_by, drop = TRUE]
  } else {
    data$group <- "All cells"
  }
  data
}

.composition_counts <- function(object, feature, group_by) {
  if (is.null(group_by) || !length(group_by) || identical(group_by, "None")) {
    return(data.frame(group = character(), level = character(), cells = integer(),
                      stringsAsFactors = FALSE))
  }
  values <- .explorer_feature_values(object, feature)
  metadata <- object[[]]
  cells <- intersect(names(values), rownames(metadata))
  frame <- data.frame(group = metadata[cells, group_by, drop = TRUE],
                      value = values[cells], stringsAsFactors = FALSE)
  if (is.numeric(frame$value)) {
    frame$level <- cut(frame$value, breaks = 20)
    dplyr::count(frame, group, level, name = "cells")
  } else {
    dplyr::count(frame, group, value, name = "cells") |>
      dplyr::rename(level = value)
  }
}

.build_sparks_targets_plan <- function(base_config_path, override_config_path,
                                       results_dir) {
  if (!requireNamespace("targets", quietly = TRUE)) {
    stop("Package 'targets' is required to build this workflow.")
  }
  cfg <- load_pipeline_config(base_config_path, override_config_path)
  sample_table <- load_sample_table(cfg)
  sample_table <- .targets_sample_table(cfg, sample_table)
  targets <- list(
    targets::tar_target_raw("results_dir", results_dir),
    targets::tar_target_raw("base_config", base_config_path, format = "file"),
    targets::tar_target_raw("override_config",
      if (is.null(override_config_path)) base_config_path else
        override_config_path, format = "file"),
    targets::tar_target_raw("config", quote(SPARKS::load_pipeline_config(
      base_config, if (identical(override_config, base_config)) NULL else
        override_config))),
    targets::tar_target_raw("analysis_config", quote(
      SPARKS:::.analysis_stage_config(config, results_dir))),
    targets::tar_target_raw("sample_config", quote(
      SPARKS:::.sample_stage_config(config))),
    targets::tar_target_raw("sample_table_file", quote(
      SPARKS:::.sample_table_source(config, base_config)), format = "file"),
    targets::tar_target_raw("sample_table", quote(
      SPARKS:::.load_sample_table_target(config, sample_table_file)))
  )
  sample_names <- character(nrow(sample_table))
  for (i in seq_len(nrow(sample_table))) {
    input_name <- sprintf("sample_input_%05d", i)
    name <- sprintf("sample_%05d", i)
    sample_names[i] <- name
    input_command <- substitute(
      SPARKS:::.sample_input_files(sample_config, sample_table, ROW),
      list(ROW = as.integer(i))
    )
    targets[[length(targets) + 1L]] <- targets::tar_target_raw(
      input_name, input_command, format = "file"
    )
    command <- substitute(
      SPARKS:::.process_sample_target(sample_config, sample_table, ROW,
                                      results_dir, INPUT_FILES),
      list(ROW = as.integer(i), INPUT_FILES = as.name(input_name))
    )
    targets[[length(targets) + 1L]] <- targets::tar_target_raw(
      name, command, format = "file"
    )
  }
  group_names <- unique(sample_table$comparison_group)
  group_keys <- .target_group_key(group_names)
  if (anyDuplicated(group_keys)) {
    stop("Comparison group names collide after safe target-name conversion.",
         call. = FALSE)
  }
  group_targets <- character(length(group_names))
  report_targets <- character(length(group_names))
  for (i in seq_along(group_names)) {
    group <- group_names[i]
    group_key <- group_keys[i]
    name <- paste0("group_", group_key)
    group_targets[i] <- name
    sample_symbols <- lapply(.targets_sample_rows(sample_table, group),
                             function(row) as.name(sample_names[row]))
    command <- substitute(
      SPARKS:::.run_sparks_group_target(
        base_config, analysis_config, sample_table, GROUP,
        SAMPLE_FILES, results_dir
      ),
      list(GROUP = group,
           SAMPLE_FILES = as.call(c(list(as.name("list")), sample_symbols)))
    )
    targets[[length(targets) + 1L]] <- targets::tar_target_raw(
      name, command, format = "file"
    )
  }
  existing_target_names <- vapply(targets, `[[`, character(1), "name")
  for (i in seq_along(group_names)) {
    group <- group_names[i]
    name <- paste0("reports_", group_keys[i])
    previous_reports <- if (i > 1L) report_targets[seq_len(i - 1L)] else
      character()
    if (name %in% c(existing_target_names, previous_reports)) {
      stop("Generated target name collision: ", name, call. = FALSE)
    }
    report_targets[i] <- name
    command <- substitute(
      SPARKS:::.generate_group_reports_target(
        analysis_config, GROUP, GROUP_MARKER
      ),
      list(GROUP = group,
           GROUP_MARKER = as.name(group_targets[i]))
    )
    targets[[length(targets) + 1L]] <- targets::tar_target_raw(
      name, command, format = "file"
    )
  }
  report_command <- substitute(
    { invisible(REPORTS); SPARKS:::.write_consolidated_report(results_dir) },
    list(REPORTS = as.call(c(list(as.name("c")),
                             lapply(report_targets, as.name))))
  )
  targets[[length(targets) + 1L]] <- targets::tar_target_raw(
    "consolidated_report", report_command, format = "file"
  )
  manifest_command <- quote(
    SPARKS:::.write_results_manifest(results_dir, consolidated_report)
  )
  targets[[length(targets) + 1L]] <- targets::tar_target_raw(
    "results_manifest", manifest_command, format = "file"
  )
  targets
}

.sample_stage_config <- function(cfg) {
  list(
    pipeline = cfg$pipeline[c("data_dir", "results_dir")],
    input = cfg$input,
    species = cfg$species[c("gene_removal_pattern", "mt_pattern",
                            "genes_to_remove")],
    qc = cfg$qc,
    plot = cfg$plot,
    parallel = cfg$parallel
  )
}

.analysis_stage_config <- function(cfg, results_dir) {
  structure(list(cfg = cfg, results_dir = results_dir),
            class = "sparks_analysis_config")
}

.group_analysis_config <- function(cfg, results_dir) {
  cfg$pipeline$results_dir <- results_dir
  cfg$pipeline$independent_samples <- FALSE
  cfg
}

.targets_sample_rows <- function(sample_table, group) {
  which(sample_table$comparison_group == group)
}

.targets_sample_table <- function(cfg, sample_table) {
  if (isTRUE(cfg$pipeline$independent_samples)) {
    sample_table$comparison_group <- if (!anyDuplicated(sample_table$folder_id)) {
      sample_table$folder_id
    } else {
      paste(sample_table$comparison_group, sample_table$folder_id, sep = "_")
    }
  }
  sample_table
}

.sample_input_files <- function(cfg, sample_table, row) {
  sample_dir <- file.path(cfg$pipeline$data_dir,
                          sample_table$folder_id[row])
  if (!dir.exists(sample_dir)) {
    stop("Sample input directory not found: ", sample_dir, call. = FALSE)
  }
  paths <- list.files(sample_dir, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE, include.dirs = FALSE)
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  if (!length(paths)) {
    stop("No sample input files found under: ", sample_dir, call. = FALSE)
  }
  normalizePath(paths, winslash = "/", mustWork = TRUE)
}

.sample_table_source <- function(cfg, base_config_path) {
  source <- cfg$pipeline$sample_table
  if (is.character(source) && length(source) == 1L) {
    if (file.exists(source)) return(normalizePath(source, winslash = "/"))
    config_dirs <- c(cfg$pipeline$override_config_dir,
                     cfg$pipeline$config_dir)
    config_dirs <- config_dirs[!vapply(config_dirs, is.null, logical(1))]
    ancestors <- unlist(lapply(config_dirs, function(path) {
      path <- normalizePath(path, winslash = "/", mustWork = FALSE)
      dirs <- character()
      repeat {
        dirs <- c(dirs, path)
        parent <- dirname(path)
        if (identical(parent, path)) break
        path <- parent
      }
      dirs
    }), use.names = FALSE)
    candidates <- file.path(unique(ancestors), source)
    found <- candidates[file.exists(candidates)]
    if (length(found)) return(normalizePath(found[1L], winslash = "/"))
  }
  # Inline YAML sample tables are already tracked through the config file.
  normalizePath(base_config_path, winslash = "/", mustWork = TRUE)
}

.load_sample_table_target <- function(cfg, sample_table_file) {
  invisible(sample_table_file)
  load_sample_table(cfg)
}

.process_sample_target <- function(cfg, sample_table, row, results_dir,
                                   input_files = NULL) {
  invisible(input_files)
  meta <- sample_table[row, , drop = FALSE]
  qf <- if ("quant_format" %in% names(meta) &&
            !is.na(meta$quant_format[1L]) && nzchar(meta$quant_format[1L])) {
    as.character(meta$quant_format[1L])
  } else NULL
  sample_dir <- file.path(results_dir, "_targets", "samples")
  make_dir(sample_dir)
  path <- file.path(sample_dir, paste0(sprintf("%05d", row), ".rds"))
  cfg$pipeline$results_dir <- results_dir
  sample_label <- .safe_filename(meta$folder_id[1L])
  object <- process_single_sample(
    folder_id = meta$folder_id[1L], protocol = meta$protocol[1L],
    file_prefix = paste0(meta$protocol[1L], "_", meta$folder_id[1L]),
    qc_dir = file.path(results_dir, "Sample_QC", sample_label),
    data_path = cfg$pipeline$data_dir,
    gene_removal_pattern = cfg$species$gene_removal_pattern,
    mt_pattern = cfg$species$mt_pattern,
    genes_to_remove = cfg$species$genes_to_remove,
    min_features = cfg$qc$min_features, max_features = cfg$qc$max_features,
    max_counts = cfg$qc$max_counts, max_mt_percent = cfg$qc$max_mt_percent,
    min_cells = cfg$qc$min_cells, cfg = cfg,
    raw_matrix_dir = file.path(results_dir, "RawMatrix", sample_label),
    quant_format = qf
  )
  if (is.null(object)) {
    saveRDS(NULL, path)
    return(.sample_target_outputs(path, results_dir, sample_label))
  }
  saveRDS(object, path)
  .sample_target_outputs(path, results_dir, sample_label)
}

.target_sample_paths <- function(sample_file) {
  paths <- if (is.character(sample_file)) sample_file else
    vapply(sample_file, function(x) {
      if (!length(x)) NA_character_ else as.character(x[[1L]])
    }, character(1))
  paths[!is.na(paths) & file.exists(paths)]
}

.sample_target_outputs <- function(rds_path, results_dir, sample_label) {
  sidecar_dirs <- file.path(results_dir, c("Sample_QC", "RawMatrix"),
                            sample_label)
  sidecars <- unlist(lapply(sidecar_dirs, function(path) {
    if (!dir.exists(path)) return(character())
    files <- list.files(path, recursive = TRUE, full.names = TRUE)
    files[!file.info(files)$isdir]
  }), use.names = FALSE)
  normalizePath(unique(c(rds_path, sidecars)), winslash = "/", mustWork = TRUE)
}

.load_cached_samples <- function(processed_samples, rows) {
  lapply(processed_samples[rows], function(item) {
    if (is.character(item) && length(item) == 1L) {
      if (!file.exists(item)) stop("Processed sample target is missing: ", item,
                                   call. = FALSE)
      readRDS(item)
    } else item
  })
}

.write_results_manifest <- function(results_dir, report_path) {
  catalog <- .result_catalog(results_dir)
  manifest <- file.path(results_dir, "results_manifest.tsv")
  utils::write.table(catalog, manifest, sep = "\t", quote = FALSE,
                     row.names = FALSE, na = "")
  invisible(report_path)
  normalizePath(manifest, winslash = "/", mustWork = TRUE)
}

.run_sparks_group_target <- function(base_config_path, analysis_config,
                                     sample_table, group_key, sample_files,
                                     results_dir) {
  safe_group <- gsub("[^A-Za-z0-9._-]+", "_", group_key)
  group_root <- file.path(results_dir, group_key)
  make_dir(group_root)

  analysis_override <- .group_analysis_config(analysis_config$cfg,
                                              analysis_config$results_dir)
  override_dir <- file.path(results_dir, "_targets", "overrides")
  make_dir(override_dir)
  group_override <- file.path(override_dir, paste0(safe_group, ".yaml"))
  yaml::write_yaml(analysis_override, group_override)

  cached_paths <- .target_sample_paths(sample_files)
  selected <- sample_table$comparison_group == group_key
  if (length(cached_paths) != sum(selected)) {
    stop("Expected one completed sample target for every sample in group.",
         call. = FALSE)
  }
  missing_samples <- vapply(cached_paths, function(path) is.null(readRDS(path)),
                            logical(1))
  if (any(missing_samples)) {
    stop("Sample processing returned no usable data for group '", group_key,
         "' (", sum(missing_samples), " skipped sample(s)).", call. = FALSE)
  }
  group_table <- sample_table[selected, , drop = FALSE]
  sparks(
    base_config_path = base_config_path,
    override_config_path = group_override,
    sample_metadata = group_table,
    processed_samples = as.list(cached_paths),
    render_reports = FALSE
  )
  marker <- file.path(results_dir, "_targets",
                      paste0("group-complete-", safe_group, ".txt"))
  writeLines(group_root, marker)
  completed_files <- c(marker, .group_result_files(results_dir, group_key))
  completed_files <- unique(completed_files[file.exists(completed_files)])
  normalizePath(completed_files, winslash = "/", mustWork = TRUE)
}

.generate_group_reports_target <- function(analysis_config, group_key,
                                           group_marker) {
  invisible(group_marker)
  cfg <- analysis_config$cfg
  cfg$pipeline$results_dir <- analysis_config$results_dir
  dirs <- .setup_group_dirs(cfg$pipeline$results_dir, group_key)
  group_root <- dirs$base
  object_path <- file.path(dirs$rdata, paste0("Merged_", group_key, ".rds"))
  if (!file.exists(object_path)) {
    marker <- file.path(cfg$pipeline$results_dir, "_targets",
                        paste0("reports-skipped-", .safe_filename(group_key), ".txt"))
    writeLines("No processed object was produced.", marker)
    return(normalizePath(marker, winslash = "/", mustWork = TRUE))
  }
  group_files <- list.files(dirs$base, recursive = TRUE, full.names = TRUE)
  object <- readRDS(object_path)
  qc_template <- .find_report_template(cfg$report$rmd_template,
                                      cfg$pipeline$config_dir, "qc_report.Rmd")
  results_template <- .find_report_template(
    cfg$report$results_rmd_template, cfg$pipeline$config_dir,
    "results_report.Rmd"
  )
  reports <- character()
  if (!is.null(qc_template)) {
    reports <- c(reports, generate_qc_report(
      seurat_obj = object, comp_group = group_key, out_dir = dirs$qc,
      author = cfg$report$author %||% "Pipeline User",
      title = cfg$report$title %||% paste("QC Report -", group_key),
      rmd_template = qc_template, cfg = cfg, run_files = group_files,
      log_path = file.path(cfg$pipeline$results_dir, "pipeline.log")
    ))
  }
  if (!is.null(results_template)) {
    groupings <- cfg$groupings_main
    if (is.null(groupings)) {
      groupings <- unique(c(cfg$processing$cluster_col, cfg$singler$label_names))
    }
    reports <- c(reports, generate_results_report(
      seurat_obj = object, comp_group = group_key, out_dir = group_root,
      groupings = intersect(groupings, colnames(object@meta.data)),
      author = cfg$report$author %||% "Pipeline User",
      title = cfg$report$title %||% paste("Results Report -", group_key),
      rmd_template = results_template, cfg = cfg, run_files = group_files,
      log_path = file.path(cfg$pipeline$results_dir, "pipeline.log")
    ))
  }
  marker <- file.path(cfg$pipeline$results_dir, "_targets",
                      paste0("reports-complete-", .safe_filename(group_key), ".txt"))
  writeLines(paste(reports[file.exists(reports)], collapse = "\n"), marker)
  completed_files <- c(marker, reports[file.exists(reports)])
  completed_files <- unique(completed_files[file.exists(completed_files)])
  normalizePath(completed_files, winslash = "/", mustWork = TRUE)
}

.group_result_files <- function(results_dir, group_root) {
  path <- file.path(results_dir, group_root)
  paths <- list.files(path, recursive = TRUE, full.names = TRUE)
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  paths[basename(paths) != "targets_override.yaml"]
}

.result_group_key <- function(group) gsub("[^A-Za-z0-9._-]+", "_", group)

.target_group_key <- function(group) gsub("[^A-Za-z0-9_]+", "_", group)

.result_catalog <- function(results_dir) {
  if (!dir.exists(results_dir)) {
    stop("Results directory does not exist: ", results_dir, call. = FALSE)
  }
  paths <- list.files(results_dir, recursive = TRUE, full.names = TRUE,
                      all.files = FALSE, include.dirs = FALSE)
  paths <- paths[file.info(paths)$isdir %in% FALSE]
  results_root <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  paths <- paths[!grepl("(^|[/\\\\])_targets([/\\\\]|$)", paths)]
  paths <- paths[!grepl("(^|[/\\\\])_targets_override_[^/\\\\]+\\.yaml$", paths)]
  paths <- paths[!basename(paths) %in%
                   c("targets_override.yaml", "shiny_override.yaml")]
  paths <- paths[!basename(paths) %in%
                   c("SPARKS_all_results.html", "SPARKS_large_results.zip",
                     "results_manifest.tsv")]
  if (!length(paths)) {
    return(data.frame(name = character(), path = character(), type = character(),
                      size = numeric(), stringsAsFactors = FALSE))
  }
  ext <- tolower(tools::file_ext(paths))
  type <- ifelse(ext %in% c("csv", "tsv", "txt"), "Table",
                 ifelse(ext %in% c("png", "jpg", "jpeg", "gif", "webp", "svg", "pdf"),
                        "Plot", ifelse(ext %in% c("html", "htm"), "Report", "File")))
  data.frame(
    name = substring(paths, nchar(results_root) + 2L),
    path = paths,
    type = type,
    size = unname(file.info(paths)$size),
    stringsAsFactors = FALSE
  )
}

#' List downloadable SPARKS results
#'
#' @param results_dir Pipeline results directory.
#' @return A data frame with relative names, absolute paths, types, and sizes.
#' @export
sparks_results <- function(results_dir) .result_catalog(results_dir)

.write_consolidated_report <- function(results_dir, dependencies = NULL) {
  invisible(dependencies)
  for (pkg in c("htmltools", "base64enc")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The consolidated report requires package '", pkg,
           "'. Install it with install.packages('", pkg, "').", call. = FALSE)
    }
  }
  catalog <- .result_catalog(results_dir)
  output <- file.path(results_dir, "SPARKS_all_results.html")
  data_uri <- function(path) {
    ext <- tolower(tools::file_ext(path))
    mime <- switch(ext, png = "image/png", jpg = "image/jpeg", jpeg = "image/jpeg",
                   gif = "image/gif", webp = "image/webp",
                   svg = "image/svg+xml", pdf = "application/pdf",
                   csv = "text/csv", tsv = "text/tab-separated-values",
                   txt = "text/plain", html = "text/html", htm = "text/html",
                   "application/octet-stream")
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    paste0("data:", mime, ";base64,", base64enc::base64encode(raw))
  }
  html_table <- function(x) {
    htmltools::tags$table(
      class = "result-table",
      htmltools::tags$thead(htmltools::tags$tr(lapply(names(x), htmltools::tags$th))),
      htmltools::tags$tbody(lapply(seq_len(nrow(x)), function(i) {
        htmltools::tags$tr(lapply(x[i, , drop = FALSE], function(cell) {
          htmltools::tags$td(as.character(cell))
        }))
      }))
    )
  }
  tabs <- if (!nrow(catalog)) {
    list(htmltools::tags$p("No result files were found."))
  } else {
    lapply(seq_len(nrow(catalog)), function(i) {
      row <- catalog[i, ]
      ext <- tolower(tools::file_ext(row$path))
      download_link <- htmltools::tags$a("Download original",
                                        href = data_uri(row$path),
                                        download = basename(row$path))
      content <- if (ext %in% c("png", "jpg", "jpeg", "gif", "webp", "svg")) {
        htmltools::tags$img(src = data_uri(row$path),
                            style = "max-width:100%;height:auto")
      } else if (ext %in% c("csv", "tsv", "txt")) {
        header <- tryCatch(readLines(row$path, n = 1L, warn = FALSE),
                          error = function(e) character())
        tabs <- if (length(header)) lengths(regmatches(header,
          gregexpr("\t", header, fixed = TRUE))) else 0L
        commas <- if (length(header)) lengths(regmatches(header,
          gregexpr(",", header, fixed = TRUE))) else 0L
        sep <- if (ext == "csv" || commas > tabs) "," else "\t"
        tab <- tryCatch(utils::read.delim(row$path, sep = sep,
                                           check.names = FALSE,
                                           stringsAsFactors = FALSE,
                                           nrows = 1000),
                        error = function(e) NULL)
        if (is.null(tab)) {
          htmltools::tags$pre(paste(readLines(row$path, warn = FALSE),
                                    collapse = "\n"))
        } else {
          csv_path <- tempfile(fileext = ".csv")
          utils::write.csv(tab, csv_path, row.names = FALSE, na = "")
          csv_uri <- data_uri(csv_path)
          unlink(csv_path)
          htmltools::tags$div(
            htmltools::tags$input(type = "search", class = "table-filter",
                                  placeholder = "Filter rows"),
            htmltools::tags$a("Download preview as CSV", href = csv_uri,
                              download = paste0(tools::file_path_sans_ext(
                                basename(row$path)), ".csv")),
            html_table(tab)
          )
        }
      } else if (ext %in% c("html", "htm")) {
        htmltools::tags$iframe(srcdoc = paste(readLines(row$path, warn = FALSE),
                                              collapse = "\n"),
                               style = "width:100%;height:600px;border:0")
      } else {
        htmltools::tags$p("Preview unavailable; use the download link.")
      }
      htmltools::tags$section(
        id = paste0("result-", i), class = "result",
        htmltools::tags$h3(row$name),
        download_link,
        content
      )
    })
  }
  head <- htmltools::tagList(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("SPARKS results"),
      htmltools::tags$style(htmltools::HTML(
        "body{font-family:system-ui;margin:2rem}nav{display:flex;gap:.5rem;flex-wrap:wrap}.result{border-top:1px solid #ddd;padding:1rem 0}.table-filter{margin:1rem 0}table{border-collapse:collapse}td,th{border:1px solid #ccc;padding:.3rem}.result-tab{display:none}.result-tab.active{display:block}nav button{cursor:pointer}"
      )),
      htmltools::tags$script(htmltools::HTML(
        "document.addEventListener('click',e=>{if(e.target.matches('[data-tab]')){document.querySelectorAll('.result-tab').forEach(x=>x.classList.remove('active'));let tab=document.getElementById(e.target.dataset.tab);if(tab)tab.classList.add('active')}});document.addEventListener('input',e=>{if(e.target.matches('.table-filter')){let q=e.target.value.toLowerCase(),t=e.target.parentElement.querySelector('table');if(t)[...t.tBodies[0].rows].forEach(r=>r.hidden=!r.innerText.toLowerCase().includes(q))}});"
      ))
  )
  body <- htmltools::tags$body(
      htmltools::tags$h1("SPARKS — all results"),
      htmltools::tags$p("Each catalogued output is embedded below with an original-file download. Tabular previews show up to 1,000 rows; large result collections can make this document large."),
      if (nrow(catalog)) {
        htmltools::tagList(
          htmltools::tags$nav(lapply(seq_len(nrow(catalog)), function(i)
            htmltools::tags$button(catalog$name[i],
                                   `data-tab` = paste0("result-", i)))),
          lapply(seq_len(nrow(catalog)), function(i) {
            tab <- tabs[[i]]
            tab$attribs$class <- if (i == 1L) "result-tab active" else
              "result-tab"
            tab
          })
        )
      } else tabs
    )
  head_html <- htmltools::renderTags(head)$html
  body_html <- htmltools::renderTags(body)$html
  writeLines(c("<!DOCTYPE html>", "<html>", "<head>", head_html,
               "</head>", body_html, "</html>"), output, useBytes = TRUE)
  normalizePath(output, winslash = "/", mustWork = TRUE)
}

#' Generate one portable HTML document for all pipeline results
#'
#' @param results_dir Pipeline results directory.
#' @param output_file Optional output path. Defaults to
#'   \,\code{SPARKS_all_results.html} inside \,\code{results_dir}.
#' @return The generated HTML path.
#' @export
generate_all_results_report <- function(results_dir,
                                        output_file = file.path(results_dir,
                                                               "SPARKS_all_results.html")) {
  path <- .write_consolidated_report(results_dir)
  if (!identical(normalizePath(output_file, winslash = "/", mustWork = FALSE),
                 path)) {
    make_dir(dirname(output_file))
    if (dir.exists(output_file)) {
      stop("output_file points to a directory: ", output_file, call. = FALSE)
    }
    if (!file.copy(path, output_file, overwrite = TRUE)) {
      stop("Unable to copy consolidated report to: ", output_file,
           call. = FALSE)
    }
    path <- normalizePath(output_file, winslash = "/", mustWork = TRUE)
  }
  path
}

#' Launch the SPARKS Shiny interface
#'
#' @param ... Reserved for future app-launch options.
#' @return A Shiny application object.
#' @export
sparks_app <- function(...) {
  for (pkg in c("shiny", "bslib", "DT", "base64enc")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("The SPARKS interface requires package '", pkg,
           "'. Install it with install.packages('", pkg, "').",
           call. = FALSE)
    }
  }
  ui <- bslib::page_navbar(
    title = "SPARKS",
    theme = bslib::bs_theme(version = 5),
    bslib::nav_panel("Configure and run",
      shiny::sidebarLayout(
        shiny::sidebarPanel(
          shiny::textInput("base_config", "Base YAML path"),
          shiny::textInput("override_config", "Override YAML path (optional)"),
          shiny::fileInput("sample_file", "Sample table (optional TSV)"),
          shiny::textInput("results_dir", "Results directory", "results"),
          shiny::selectInput("species", "Species", c("Mouse", "Human")),
          shiny::numericInput("min_features", "QC: minimum features", 200, min = 0),
          shiny::numericInput("max_mt", "QC: maximum mitochondrial percent", 10, min = 0),
          shiny::numericInput("resolution", "Clustering resolution", 0.5, min = 0),
          shiny::checkboxInput("run_deg", "Run differential expression", TRUE),
          shiny::checkboxInput("run_pseudobulk", "Run pseudobulk", FALSE),
          shiny::textAreaInput("yaml_text", "Complete override YAML",
                               rows = 16, placeholder = "Paste YAML here"),
          shiny::actionButton("validate", "Check configuration only"),
          shiny::selectizeInput("target_names", "Targets to run (empty = all)",
                                choices = NULL, multiple = TRUE,
                                options = list(placeholder = "All targets")),
          shiny::actionButton("run", "Run / resume targets", class = "btn-primary"),
          shiny::actionButton("stop", "Stop active run")
        ),
            shiny::mainPanel(
            bslib::layout_columns(
              bslib::value_box("Run state", shiny::uiOutput("run_state_box"),
                               theme = "primary"),
              bslib::value_box("Elapsed", shiny::textOutput("elapsed"),
                               theme = "info"),
              bslib::value_box("Target counts", shiny::textOutput("target_counts"),
                               theme = "success"),
              col_widths = c(4, 4, 4)
            ),
            bslib::card(
              bslib::card_header("Run status and validation"),
              shiny::verbatimTextOutput("status"),
              shiny::verbatimTextOutput("validation")
            ),
            bslib::card(
              bslib::card_header("Run log"),
              shiny::verbatimTextOutput("run_log")
            )
          )
      )
    ),
      bslib::nav_panel("Overview",
        shiny::sidebarLayout(
          shiny::sidebarPanel(
            shiny::textInput("overview_dir", "Results directory", "results"),
            shiny::actionButton("overview_refresh", "Refresh overview"),
            shiny::uiOutput("overview_group_choice"),
            shiny::uiOutput("overview_feature_choice"),
            shiny::uiOutput("overview_split_choice")
          ),
          shiny::mainPanel(
            bslib::layout_columns(
              bslib::value_box("Files", shiny::textOutput("overview_files"),
                               theme = "primary"),
              bslib::value_box("Groups", shiny::textOutput("overview_groups"),
                               theme = "info"),
              bslib::value_box("Total size", shiny::textOutput("overview_size"),
                               theme = "success"),
              col_widths = c(4, 4, 4)
            ),
            bslib::navset_card_tab(
              bslib::nav_panel("Summary", DT::DTOutput("overview_summary")),
              bslib::nav_panel("Dimensional reduction",
                shiny::uiOutput("overview_umap_ui"),
        shiny::uiOutput("overview_download_ui"),
        DT::DTOutput("overview_selected_cells")
              ),
              bslib::nav_panel("Gene expression",
                shiny::uiOutput("overview_expression_ui"),
                shiny::uiOutput("overview_gene_violin_ui"),
                DT::DTOutput("overview_expression_values")
              ),
              bslib::nav_panel("Composition",
                shiny::uiOutput("overview_composition_ui"),
                DT::DTOutput("overview_composition_table")
              ),
              bslib::nav_panel("Files", DT::DTOutput("overview_catalog"))
            )
          )
        )
      ),
    bslib::nav_panel("Results",
      shiny::sidebarLayout(
        shiny::sidebarPanel(shiny::textInput("browse_dir", "Results directory"),
                            shiny::actionButton("refresh", "Refresh"),
                            shiny::uiOutput("download_choice"),
                            shiny::selectInput("download_format", "Download format",
                                               choices = c("Original", "CSV", "TSV"))),
        shiny::mainPanel(DT::DTOutput("catalog"),
                         DT::DTOutput("table_preview"),
                         shiny::uiOutput("plot_preview"),
                         shiny::downloadButton("download", "Download selected result"))
      )
    )
  )

  server <- function(input, output, session) {
    state <- shiny::reactiveValues(process = NULL, report = "Ready",
                                   started_at = NULL, log_path = NULL,
                                   validated = FALSE, last_heartbeat = NULL)
    shiny::observeEvent(input$run, {
      if (!is.null(state$process) && state$process$is_alive()) {
        state$report <- "A targets run is already active."
        return()
      }
      state$report <- "VALIDATING — checking config, sample table and target selection."
      state$validated <- FALSE
      state$started_at <- NULL
      if (is.null(input$results_dir) || !nzchar(input$results_dir)) {
        state$report <- "VALIDATION FAILED: select a results directory."
        return()
      }
      dir.create(input$results_dir, recursive = TRUE, showWarnings = FALSE)
      state$log_path <- file.path(input$results_dir, "shiny_run.log")
      .log_run_event(state$log_path, "VALIDATING",
                     "Checking configuration and sample table")
      tryCatch({
        if (is.null(input$base_config) || !nzchar(input$base_config)) {
          stop("Select a base YAML path.")
        }
        if (is.null(input$results_dir) || !nzchar(input$results_dir)) {
          stop("Select a results directory.")
        }
        dir.create(input$results_dir, recursive = TRUE, showWarnings = FALSE)
        temp_override <- .write_ui_targets_override(input, input$results_dir)
        selected_targets <- if (length(input$target_names)) {
          input$target_names
        } else NULL
        validation <- .validate_ui_configuration(
          input$base_config, temp_override, input$results_dir, selected_targets
        )
        sample_table <- validation$sample_table
        plan_names <- validation$targets
        shiny::updateSelectizeInput(session, "target_names",
                                    choices = plan_names, selected = character(),
                                    server = TRUE)
        state$validated <- TRUE
        state$report <- paste0("VALIDATED — ", nrow(sample_table),
                               " sample(s), ", length(plan_names),
                               " targets ready. Launching analysis.")
        if (!is.null(state$process) && !state$process$is_alive()) {
          state$process <- NULL
        }
        .log_run_event(state$log_path, "VALIDATED",
                       paste(nrow(sample_table), "samples;",
                             length(plan_names), "targets"))
      }, error = function(e) {
        state$validated <- FALSE
        state$report <- paste("VALIDATION FAILED:", conditionMessage(e))
        if (!is.null(input$results_dir) && nzchar(input$results_dir)) {
          dir.create(input$results_dir, recursive = TRUE, showWarnings = FALSE)
          state$log_path <- file.path(input$results_dir, "shiny_run.log")
          .log_run_event(state$log_path, "VALIDATION FAILED",
                         conditionMessage(e))
        }
      })
      if (!isTRUE(state$validated)) return()
      if (!requireNamespace("callr", quietly = TRUE)) {
        state$report <- "Install package 'callr' to run targets in the background."
        if (!is.null(state$log_path)) {
          .log_run_event(state$log_path, "NOT STARTED", "Missing package 'callr'")
        }
        return()
      }
      if (is.null(state$process) || !state$process$is_alive()) {
        override <- file.path(input$results_dir, "shiny_override.yaml")
        selected_targets <- if (length(input$target_names)) input$target_names else NULL
        state$started_at <- Sys.time()
        state$last_heartbeat <- Sys.time()
        .log_run_event(state$log_path, "STARTING",
                       paste("base config:", input$base_config,
                             "| target selection:",
                             if (is.null(selected_targets)) "all" else
                               paste(selected_targets, collapse = ",")))
        package_path <- getNamespaceInfo(asNamespace("SPARKS"), "path")
        state$process <- tryCatch(callr::r_bg(
          function(base, over, out, selected, log_path, package_path,
                   library_paths) {
            .libPaths(unique(c(library_paths, .libPaths())))
            log <- function(event, detail = "") {
              cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), event,
                  detail, "\n", file = log_path, append = TRUE, sep = "")
            }
            log("RUNNING", "Starting targets workflow")
            tryCatch({
              if (!requireNamespace("SPARKS", quietly = TRUE)) {
                if (!file.exists(file.path(package_path, "DESCRIPTION")) ||
                    !requireNamespace("pkgload", quietly = TRUE)) {
                  stop("Package 'SPARKS' is unavailable in the background R session. Install SPARKS or run from a development session with 'pkgload' installed.",
                       call. = FALSE)
                }
                pkgload::load_all(package_path, quiet = TRUE)
              }
            SPARKS::run_sparks_targets(base, over, out, selected)
              log("SUCCEEDED", "Targets workflow completed")
            }, error = function(e) {
              log("FAILED", conditionMessage(e))
              stop(e)
            })
          },
          args = list(input$base_config, override, input$results_dir,
                      selected_targets, state$log_path,
                      getNamespaceInfo(asNamespace("SPARKS"), "path"),
                      .libPaths()),
          stdout = file.path(input$results_dir, "shiny_process.log"),
          stderr = "2>&1"
        ), error = identity)
        if (inherits(state$process, "error")) {
          detail <- conditionMessage(state$process)
          state$process <- NULL
          state$started_at <- NULL
          state$report <- paste("FAILED TO START:", detail)
          .log_run_event(state$log_path, "FAILED TO START", detail)
        } else {
          state$report <- "STARTING — background process launched; waiting for target progress."
        }
      } else {
        state$report <- "A targets run is already active."
      }
    })
    shiny::observeEvent(input$validate, {
      tryCatch({
        override <- .write_ui_targets_override(input, input$results_dir)
        .validate_ui_configuration(input$base_config, override,
                                   input$results_dir)
        state$report <- "Configuration is valid. Click Run to validate again and launch."
      }, error = function(e) {
        state$report <- paste("VALIDATION FAILED:", conditionMessage(e))
      })
    })
    shiny::observeEvent(input$stop, {
      if (!is.null(state$process) && state$process$is_alive()) {
        state$process$kill()
        .log_run_event(state$log_path, "STOP REQUESTED")
        state$report <- "Stop requested for the active targets process."
      } else {
        state$report <- "No active targets process to stop."
      }
    })
    shiny::observe({
      shiny::invalidateLater(1500, session)
      if (!is.null(state$process)) {
        if (state$process$is_alive()) {
          if (is.null(state$last_heartbeat) ||
              as.numeric(difftime(Sys.time(), state$last_heartbeat,
                                  units = "secs")) >= 15) {
            .log_run_event(state$log_path, "HEARTBEAT",
                           paste("elapsed:",
                                 .format_run_elapsed(state$started_at)))
            state$last_heartbeat <- Sys.time()
          }
          state$report <- paste("RUNNING —",
                                .format_run_elapsed(state$started_at),
                                "elapsed; background process is alive.")
        } else {
          state$process$poll_io()
          code <- state$process$get_exit_status()
          state$report <- if (isTRUE(code == 0L)) {
            paste("SUCCEEDED — completed in",
                  .format_run_elapsed(state$started_at))
          } else {
            paste("FAILED or stopped — process exit code", code,
                  "after", .format_run_elapsed(state$started_at))
          }
          state$process <- NULL
        }
      }
    })
    shiny::observeEvent(input$catalog_rows_selected, {
      selected_row <- input$catalog_rows_selected
      data <- catalog()
      if (length(selected_row) == 1L && nrow(data) >= selected_row) {
        shiny::updateSelectInput(session, "selected_file",
                                 selected = data$path[selected_row])
      }
    })
    output$status <- shiny::renderText(state$report)
    output$run_state_box <- shiny::renderUI({
      color <- if (grepl("SUCCEEDED", state$report)) "#198754" else
        if (grepl("FAILED|VALIDATION FAILED", state$report)) "#dc3545" else
          if (grepl("RUNNING|STARTING", state$report)) "#0d6efd" else "#6c757d"
      shiny::tags$span(state$report, style = paste0("color:", color, ";font-weight:600"))
    })
    output$elapsed <- shiny::renderText({
      shiny::invalidateLater(1000, session)
      .format_run_elapsed(state$started_at)
    })
    output$target_counts <- shiny::renderText({
      shiny::invalidateLater(2000, session)
      store <- file.path(input$results_dir, "_targets", "store")
      .targets_run_progress(store)
    })
    output$validation <- shiny::renderText({
      path <- state$log_path %||% file.path(input$results_dir, "shiny_run.log")
      process_log <- file.path(input$results_dir, "shiny_process.log")
      lines <- c(if (file.exists(path)) readLines(path, warn = FALSE),
                 if (file.exists(process_log)) readLines(process_log, warn = FALSE))
      statuses <- tryCatch({
        store <- file.path(input$results_dir, "_targets", "store")
        if (requireNamespace("targets", quietly = TRUE) &&
            dir.exists(store)) {
          targets::tar_meta(fields = c(name, error), store = store)
        } else NULL
      }, error = function(e) NULL)
      status_text <- if (!is.null(statuses) && nrow(statuses)) {
        paste(capture.output(print(statuses)), collapse = "\n")
      } else "Target metadata will appear after the workflow initializes."
      paste(c(tail(lines, 40), status_text), collapse = "\n")
    })
    output$run_log <- shiny::renderText({
      shiny::invalidateLater(1000, session)
      paths <- c(state$log_path, file.path(input$results_dir, "shiny_process.log"))
      paths <- paths[!is.na(paths) & file.exists(paths)]
      lines <- unlist(lapply(paths, readLines, warn = FALSE), use.names = FALSE)
      paste(tail(lines, 60), collapse = "\n")
    })
    catalog <- shiny::eventReactive(input$refresh, {
      if (!nzchar(input$browse_dir) || !dir.exists(input$browse_dir)) {
        data.frame()
      } else .result_catalog(input$browse_dir)
    }, ignoreNULL = FALSE)
    output$download_choice <- shiny::renderUI({
      data <- catalog()
      shiny::selectInput("selected_file", "File", choices = data$path,
                         selected = if (nrow(data)) data$path[1] else NULL)
    })
    selected <- shiny::reactive({
      path <- input$selected_file
      if (is.null(path) || !file.exists(path)) NULL else path
    })
    output$plot_preview <- shiny::renderUI({
      path <- selected()
      if (is.null(path)) return(NULL)
      ext <- tolower(tools::file_ext(path))
      if (ext %in% c("png", "jpg", "jpeg", "gif", "webp", "svg")) {
        mime <- if (ext == "jpg") "jpeg" else ext
        shiny::tags$img(src = paste0("data:image/", mime, ";base64,",
          base64enc::base64encode(readBin(path, "raw", n = file.info(path)$size))),
          style = "max-width:100%")
      } else NULL
    })
    output$table_preview <- DT::renderDT({
      path <- selected()
      if (is.null(path) || !tolower(tools::file_ext(path)) %in%
          c("csv", "tsv", "txt")) return(NULL)
      utils::read.delim(path,
        sep = if (tolower(tools::file_ext(path)) == "csv") "," else "\t",
        nrows = 1000, check.names = FALSE)
    }, options = list(pageLength = 15, scrollX = TRUE))
    output$download <- shiny::downloadHandler(
      filename = function() {
        path <- selected()
        if (is.null(path)) return("result")
        ext <- tolower(tools::file_ext(path))
        if (ext %in% c("csv", "tsv", "txt") &&
            input$download_format %in% c("CSV", "TSV")) {
          paste0(tools::file_path_sans_ext(basename(path)),
                 if (input$download_format == "CSV") ".csv" else ".tsv")
        } else basename(path)
      },
      content = function(file) {
        path <- selected()
        .write_result_download(path, file, input$download_format)
      }
    )
    output$catalog <- DT::renderDT({
      data <- catalog()
      if (nrow(data)) {
        data$path <- NULL
        data
      } else data
    }, selection = "single", options = list(pageLength = 15, scrollX = TRUE))

    overview_catalog <- shiny::eventReactive(input$overview_refresh, {
      if (!nzchar(input$overview_dir) || !dir.exists(input$overview_dir)) {
        data.frame(name = character(), path = character(), type = character(),
                   size = numeric(), stringsAsFactors = FALSE)
      } else .result_catalog(input$overview_dir)
    }, ignoreNULL = FALSE)
    overview_objects <- shiny::eventReactive(input$overview_refresh, {
      .discover_merged_objects(input$overview_dir)
    }, ignoreNULL = FALSE)
    selected_object <- shiny::reactive({
      objects <- overview_objects()
      if (!nrow(objects)) return(NULL)
      path <- objects$path[match(input$overview_group, objects$group)]
      if (is.na(path) || !file.exists(path)) return(NULL)
      tryCatch(readRDS(path), error = function(e) NULL)
    })
    output$overview_group_choice <- shiny::renderUI({
      objects <- overview_objects()
      shiny::selectInput("overview_group", "Comparison group",
                         choices = stats::setNames(objects$group, objects$group),
                         selected = if (nrow(objects)) objects$group[1L] else NULL)
    })
    output$overview_feature_choice <- shiny::renderUI({
      object <- selected_object()
      if (is.null(object)) return(shiny::helpText("No merged Seurat objects found."))
      features <- .explorer_features(object)
      shiny::selectizeInput("overview_feature", "Color by metadata or gene",
                            choices = features,
                            selected = if (length(features)) features[1L] else NULL,
                            options = list(maxOptions = 1000L))
    })
    output$overview_split_choice <- shiny::renderUI({
      object <- selected_object()
      if (is.null(object)) return(NULL)
      metadata <- colnames(object[[]])
      shiny::selectInput("overview_split", "Split/group by metadata",
                         choices = c("None", metadata), selected = "None")
    })
    overview_metrics <- shiny::reactive({
      .result_overview_metrics(overview_catalog())
    })
    output$overview_files <- shiny::renderText(overview_metrics()$files)
    output$overview_groups <- shiny::renderText({
      objects <- overview_objects()
      if (nrow(objects)) length(unique(objects$group)) else overview_metrics()$groups
    })
    output$overview_size <- shiny::renderText({
      sprintf("%.1f MB", overview_metrics()$size_mb)
    })
    output$overview_summary <- DT::renderDT({
      DT::datatable(.summarise_results(overview_catalog()),
                    options = list(pageLength = 15, scrollX = TRUE),
                    rownames = FALSE)
    })
    output$overview_catalog <- DT::renderDT({
      data <- overview_catalog()
      data$path <- NULL
      DT::datatable(data, options = list(pageLength = 15, scrollX = TRUE),
                    rownames = FALSE)
    })
    plot_data <- shiny::reactive({
      object <- selected_object()
      shiny::req(object, input$overview_feature)
      group_by <- if (identical(input$overview_split, "None")) NULL else
        input$overview_split
      .explorer_plot_data(object, input$overview_feature, group_by)
    })
    output$overview_umap_ui <- shiny::renderUI({
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput("overview_umap_plot", height = "650px")
      } else {
        shiny::tagList(shiny::helpText(
          "Install 'plotly' for an interactive UMAP; a static plot is shown instead."),
          shiny::plotOutput("overview_umap_static", height = "650px"))
      }
    })
    output$overview_gene_violin_ui <- shiny::renderUI({
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput("overview_gene_violin_plot", height = "450px")
      } else {
        shiny::plotOutput("overview_gene_violin_static", height = "450px")
      }
    })
    output$overview_composition_ui <- shiny::renderUI({
      if (requireNamespace("plotly", quietly = TRUE)) {
        plotly::plotlyOutput("overview_composition_plot", height = "500px")
      } else {
        shiny::plotOutput("overview_composition_static", height = "500px")
      }
    })
    output$overview_expression_ui <- shiny::renderUI({
      object <- selected_object()
      if (is.null(object)) return(shiny::helpText("Choose a comparison group first."))
      genes <- rownames(object)
      shiny::selectizeInput("overview_gene", "Gene",
                            choices = genes,
                            options = list(maxOptions = 1000L,
                                           placeholder = "Search gene"))
    })
    make_umap_plot <- function(data) {
      ggplot2::ggplot(data, ggplot2::aes(x, y, color = value, key = cell,
                                        frame = group,
                                        text = paste("Cell:", cell,
                                                     "<br>Value:", value,
                                                     "<br>Group:", group))) +
        ggplot2::geom_point(size = 0.35, alpha = 0.75) +
        ggplot2::labs(x = "UMAP 1", y = "UMAP 2", color = input$overview_feature)
    }
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$overview_umap_plot <- plotly::renderPlotly({
        data <- plot_data()
        plot <- plotly::ggplotly(make_umap_plot(data), tooltip = "text",
                                 source = "umap") |>
          plotly::layout(dragmode = "lasso")
        plotly::event_register(plot, "plotly_selected")
      })
    } else {
      output$overview_umap_static <- shiny::renderPlot({
        make_umap_plot(plot_data())
      })
    }
    output$overview_expression_values <- DT::renderDT({
      object <- selected_object()
      shiny::req(object, input$overview_gene)
      group_by <- if (identical(input$overview_split, "None")) NULL else
        input$overview_split
      DT::datatable(utils::head(.gene_expression_frame(
        object, input$overview_gene, group_by
      ), 5000L),
                    options = list(pageLength = 15, scrollX = TRUE),
                    rownames = FALSE)
    })
    make_gene_violin_plot <- function() {
      object <- selected_object()
      shiny::req(object, input$overview_gene)
      group_by <- if (identical(input$overview_split, "None")) NULL else
        input$overview_split
      data <- .gene_expression_frame(object, input$overview_gene, group_by)
      ggplot2::ggplot(data, ggplot2::aes(x = group, y = expression,
                                        fill = group)) +
        ggplot2::geom_violin(scale = "width", trim = TRUE) +
        ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA) +
        ggplot2::labs(x = input$overview_split, y = input$overview_gene) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30,
                                                           hjust = 1))
    }
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$overview_gene_violin_plot <- plotly::renderPlotly({
        plotly::ggplotly(make_gene_violin_plot()) |>
          plotly::config(displayModeBar = TRUE, scrollZoom = TRUE)
      })
    } else {
      output$overview_gene_violin_static <- shiny::renderPlot({
        make_gene_violin_plot()
      })
    }
    output$overview_composition_table <- DT::renderDT({
      DT::datatable(composition_data(),
                    options = list(pageLength = 15, scrollX = TRUE),
                    rownames = FALSE)
    })
    composition_data <- shiny::reactive({
      object <- selected_object()
      shiny::req(object, input$overview_split, input$overview_feature)
      .composition_counts(object, input$overview_feature, input$overview_split)
    })
    make_composition_plot <- function() {
      counts <- composition_data()
      if (!nrow(counts)) {
        return(ggplot2::ggplot() + ggplot2::theme_void() +
                 ggplot2::labs(title = "Select metadata to split the view"))
      }
      ggplot2::ggplot(counts, ggplot2::aes(x = level, y = cells, fill = group)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::labs(x = input$overview_feature, y = "Cells",
                      fill = input$overview_split)
    }
    if (requireNamespace("plotly", quietly = TRUE)) {
      output$overview_composition_plot <- plotly::renderPlotly({
        plotly::ggplotly(make_composition_plot()) |>
          plotly::config(displayModeBar = TRUE, scrollZoom = TRUE)
      })
    } else {
      output$overview_composition_static <- shiny::renderPlot({
        make_composition_plot()
      })
    }
    selected_umap_cells <- shiny::reactive({
      object <- selected_object()
      shiny::req(object, input$overview_feature)
      if (!requireNamespace("plotly", quietly = TRUE) ||
          .explorer_feature_kind(input$overview_feature,
                                 colnames(object[[]])) == "gene") {
        return(character())
      }
      event <- plotly::event_data("plotly_selected", source = "umap")
      if (is.null(event) || !"key" %in% names(event)) return(character())
      as.character(event$key)
    })
    output$overview_selected_cells <- DT::renderDT({
      selected <- selected_umap_cells()
      if (!length(selected)) {
        return(DT::datatable(data.frame(message =
          "Lasso-select points in the dimensional reduction plot."),
          options = list(dom = "t"), rownames = FALSE))
      }
      cell_data <- plot_data()
      cells <- utils::head(cell_data[cell_data$cell %in% selected, ,
                                     drop = FALSE], 5000L)
      DT::datatable(cells, options = list(pageLength = 15, scrollX = TRUE),
                    rownames = FALSE)
    })
    output$overview_download_ui <- shiny::renderUI({
      if (length(selected_umap_cells())) {
        shiny::downloadButton("overview_download_cells", "Download selected cells")
      }
    })
    output$overview_download_cells <- shiny::downloadHandler(
      filename = function() paste0("selected_cells_", Sys.Date(), ".tsv"),
      content = function(file) {
        if (!requireNamespace("plotly", quietly = TRUE)) {
          stop("Cell selection downloads require package 'plotly'.",
               call. = FALSE)
        }
        object <- selected_object()
        shiny::req(object, input$overview_feature)
        if (.explorer_feature_kind(input$overview_feature,
                                  colnames(object[[]])) == "gene") {
          stop("Lasso selection is available for metadata coloring, not gene expression.",
               call. = FALSE)
        }
        chosen <- plot_data()[plot_data()$cell %in% selected_umap_cells(), ,
                              drop = FALSE]
        if (!nrow(chosen)) {
          stop("Lasso-select cells in the UMAP first.", call. = FALSE)
        }
        utils::write.table(chosen, file, sep = "\t", quote = FALSE,
                           row.names = FALSE, na = "")
      }
    )
  }
  shiny::shinyApp(ui, server, ...)
}
