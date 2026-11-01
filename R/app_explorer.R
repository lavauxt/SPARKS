#' SPARKS explorer: a Loupe-style Shiny app with on-the-fly re-analysis
#'
#' Open a finished run, look at the embedding, **lasso-select cells**, name them
#' as categories, and ask questions about them -- differential expression
#' (selection vs rest, or vs another category), re-clustering, signature
#' scoring, over-representation -- without re-running the pipeline. Every
#' analysis is a cached target in \code{<results_dir>/_targets_interactive}, so
#' asking twice is free and the history survives restarts. One click rebuilds
#' the unified HTML report ([sparks_report()]), which then includes those
#' analyses.
#'
#' @param results_dir Pipeline results directory (the one containing
#'   \code{<comparison>/RData/Merged_*.rds}).
#' @param object Optional path to a Seurat \code{.rds} to open first.
#' @param max_cells Cells drawn in the embedding (random down-sampling above).
#' @param background Run analyses in a separate R process. Defaults to
#'   \code{TRUE} when SPARKS is installed (not just loaded) and \code{callr}
#'   is available.
#' @param launch.browser,port,... Passed to [shiny::runApp()].
#' @return Invisibly, the value of [shiny::runApp()].
#' @export
#' @examples
#' \dontrun{
#' sparks_explorer("results")
#'
#' demo <- tempfile("demo_"); sparks_demo_results(demo)
#' sparks_explorer(demo)
#' }
sparks_explorer <- function(results_dir, object = NULL, max_cells = 60000L,
                            background = NULL, launch.browser = interactive(),
                            port = getOption("shiny.port"), ...) {
  .xp_require()
  if (!dir.exists(results_dir)) stop("Results directory does not exist: ", results_dir, call. = FALSE)
  app <- shiny::shinyApp(
    ui = bslib::page_fluid(
      theme = bslib::bs_theme(version = 5, preset = "shiny"),
      title = "SPARKS Explorer",
      shiny::tags$style(.xp_css),
      bslib::layout_columns(
        col_widths = 12,
        shiny::div(class = "d-flex align-items-center gap-3 mb-2",
                   shiny::h4("SPARKS Explorer", class = "mb-0"),
                   shiny::span(class = "text-muted small", basename(normalizePath(results_dir))))
      ),
      sparks_explorer_ui("explorer")
    ),
    server = function(input, output, session) {
      sparks_explorer_server("explorer", results_dir = results_dir, object = object,
                             max_cells = max_cells, background = background)
    }
  )
  shiny::runApp(app, launch.browser = launch.browser, port = port, ...)
}

.xp_require <- function() {
  need <- c("shiny", "bslib", "plotly", "targets", "digest")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop("The explorer requires package(s): ", paste(miss, collapse = ", "),
         ". Install with install.packages(c(", paste0("'", miss, "'", collapse = ", "), ")).",
         call. = FALSE)
  }
}

.xp_css <- "
.xp-sel { font-size: 1.05rem; } .xp-hint { color: var(--bs-secondary-color); font-size: .85rem; }
.xp-cat { display:inline-flex; align-items:center; gap:.35rem; border:1px solid var(--bs-border-color);
          border-radius:999px; padding:.1rem .6rem; margin:0 .3rem .3rem 0; font-size:.85rem; }
.xp-status-built { color:#006300; } .xp-status-error { color:#b3261e; } .xp-status-pending { color:#8a5a00; }
.xp-plot .plotly { min-height: 560px; }
"

# ---- UI -------------------------------------------------------------------------------

#' @rdname sparks_explorer
#' @param id Module id.
#' @export
sparks_explorer_ui <- function(id = "explorer") {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(8, 4),
    bslib::navset_card_tab(
      id = ns("views"),
      bslib::nav_panel(
        "Embedding",
        shiny::div(class = "d-flex flex-wrap gap-3 align-items-end mb-2",
          shiny::selectInput(ns("dataset"), "Dataset", choices = NULL, width = "220px"),
          shiny::selectInput(ns("reduction"), "Reduction", choices = NULL, width = "140px"),
          shiny::selectizeInput(ns("color_by"), "Colour by (metadata or any gene)", choices = NULL, width = "300px"),
          shiny::sliderInput(ns("point_size"), "Point size", 1, 10, 4, step = 0.5, width = "150px"),
          shiny::sliderInput(ns("opacity"), "Opacity", 0.1, 1, 0.8, step = 0.05, width = "150px")),
        shiny::div(class = "xp-plot", plotly::plotlyOutput(ns("embedding"), height = "620px")),
        shiny::div(class = "xp-hint mt-1",
                   "Use the lasso or box tool in the plot toolbar to select cells; double-click to clear.")
      ),
      bslib::nav_panel(
        "Genes",
        shiny::div(class = "d-flex flex-wrap gap-3 align-items-end mb-2",
          shiny::selectizeInput(ns("gene"), "Gene", choices = NULL, width = "220px"),
          shiny::selectInput(ns("group_by"), "Group by", choices = NULL, width = "220px"),
          shiny::textInput(ns("dot_genes"), "Dot plot genes (comma-separated)", width = "360px")),
        bslib::layout_columns(col_widths = c(6, 6),
          plotly::plotlyOutput(ns("violin"), height = "420px"),
          plotly::plotlyOutput(ns("dotplot"), height = "420px"))
      ),
      bslib::nav_panel(
        "Composition",
        shiny::div(class = "d-flex flex-wrap gap-3 align-items-end mb-2",
          shiny::selectInput(ns("comp_stack"), "Stack by", choices = NULL, width = "220px"),
          shiny::selectInput(ns("comp_across"), "Compare across", choices = NULL, width = "220px")),
        plotly::plotlyOutput(ns("composition"), height = "520px")
      ),
      bslib::nav_panel(
        "Analyses",
        shiny::div(class = "d-flex flex-wrap gap-2 align-items-center mb-2",
          shiny::actionButton(ns("refresh"), "Refresh"),
          shiny::actionButton(ns("delete"), "Delete selected"),
          shiny::actionButton(ns("build_report"), "Build unified report", class = "btn-primary"),
          shiny::uiOutput(ns("report_link"), inline = TRUE)),
        shiny::uiOutput(ns("job_status")),
        shiny::tableOutput(ns("history")),
        shiny::selectInput(ns("view_id"), "Show result", choices = NULL, width = "100%"),
        shiny::uiOutput(ns("result_ui"))
      )
    ),
    shiny::tagList(
      bslib::card(
        bslib::card_header("Selection"),
        shiny::div(class = "xp-sel", shiny::textOutput(ns("sel_summary"))),
        shiny::div(class = "d-flex gap-2 align-items-end mt-2",
          shiny::textInput(ns("cat_name"), "Save as category", placeholder = "e.g. Region A", width = "100%"),
          shiny::actionButton(ns("save_cat"), "Save")),
        shiny::div(class = "d-flex gap-2 mt-1",
          shiny::actionButton(ns("clear_sel"), "Clear selection"),
          shiny::actionButton(ns("clear_cats"), "Clear categories")),
        shiny::uiOutput(ns("cats_ui")),
        bslib::accordion(
          open = FALSE,
          bslib::accordion_panel(
            "Select by barcodes",
            shiny::fileInput(ns("import"), "Import barcodes (CSV/TXT, e.g. from the HTML report)", width = "100%"),
            shiny::textAreaInput(ns("paste"), NULL, rows = 3, placeholder = "Paste barcodes, one per line", width = "100%"),
            shiny::actionButton(ns("apply_paste"), "Select pasted cells"))
        )
      ),
      bslib::card(
        bslib::card_header("Analyse"),
        bslib::accordion(
          open = "de",
          bslib::accordion_panel(
            "Differential expression", value = "de",
            shiny::selectInput(ns("de_vs"), "Compare selection against", choices = c("All other cells" = "rest")),
            shiny::div(class = "d-flex gap-2",
              shiny::numericInput(ns("de_min_pct"), "min.pct", 0.1, min = 0, max = 1, step = 0.05, width = "50%"),
              shiny::numericInput(ns("de_lfc"), "logfc.threshold", 0.25, min = 0, step = 0.05, width = "50%")),
            shiny::selectInput(ns("de_test"), "Test", c("wilcox", "t", "bimod"), width = "100%"),
            shiny::actionButton(ns("run_de"), "Find differentially expressed genes", class = "btn-primary")),
          bslib::accordion_panel(
            "Re-cluster selection", value = "recluster",
            shiny::sliderInput(ns("rc_res"), "Resolution", 0.1, 2, 0.5, step = 0.1, width = "100%"),
            shiny::actionButton(ns("run_recluster"), "Re-cluster the selected cells", class = "btn-primary")),
          bslib::accordion_panel(
            "Score a gene signature", value = "signature",
            shiny::selectInput(ns("sig_set"), "From configuration", choices = c("(custom)" = "")),
            shiny::textAreaInput(ns("sig_genes"), "Genes (comma / space separated)", rows = 2, width = "100%"),
            shiny::actionButton(ns("run_sig"), "Score all cells", class = "btn-primary")),
          bslib::accordion_panel(
            "Over-representation", value = "ora",
            shiny::selectInput(ns("ora_from"), "Gene list", choices = c("Select an analysis in History" = "")),
            shiny::fileInput(ns("gmt"), "Gene sets (GMT, optional)", width = "100%"),
            shiny::actionButton(ns("run_ora"), "Test enrichment", class = "btn-primary"))
        )
      )
    )
  )
}

# ---- Server -----------------------------------------------------------------------------

#' @rdname sparks_explorer
#' @param results_dir Results directory (module argument).
#' @export
sparks_explorer_server <- function(id = "explorer", results_dir, object = NULL,
                                   max_cells = 60000L, background = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    results_dir <- normalizePath(results_dir, winslash = "/", mustWork = TRUE)
    bundle <- sparks_read_results(results_dir)
    objs <- .discover_result_objects(bundle)
    if (!is.null(object)) {
      objs <- rbind(data.frame(label = tools::file_path_sans_ext(basename(object)), kind = "custom",
                               comparison = NA_character_,
                               path = normalizePath(object, winslash = "/", mustWork = TRUE),
                               stringsAsFactors = FALSE), objs)
    }
    use_bg <- if (is.null(background)) .xp_can_background() else isTRUE(background) && .xp_can_background()
    gene_sets_cfg <- sparks_config_gene_sets(bundle$config)

    rv <- shiny::reactiveValues(
      sel = character(), categories = list(), extra_meta = list(),
      job = NULL, tick = 0L, report = NULL, flash = NULL, gmt = NULL
    )
    src <- ns("embedding")

    # ---- Data -------------------------------------------------------------------------
    shiny::observe({
      ch <- stats::setNames(objs$path, paste0(objs$label, " (", objs$kind, ")"))
      shiny::updateSelectInput(session, "dataset", choices = ch)
    })
    obj_path <- shiny::reactive({ shiny::req(input$dataset); input$dataset })
    obj <- shiny::reactive({
      path <- obj_path()
      shiny::withProgress(message = "Loading object", value = 0.4, .load_object_cached(path))
    })
    shiny::observeEvent(obj_path(), {
      rv$sel <- character(); rv$categories <- list(); rv$extra_meta <- list()
      o <- obj()
      red <- names(o@reductions)
      shiny::updateSelectInput(session, "reduction", choices = red,
                               selected = if ("umap" %in% red) "umap" else red[[1L]])
      meta <- .xp_meta_columns(o)
      shiny::updateSelectizeInput(session, "color_by", server = TRUE,
        choices = c(meta$cat, meta$num, rownames(o)), selected = .xp_default_group(meta$cat))
      shiny::updateSelectizeInput(session, "gene", server = TRUE, choices = rownames(o),
                                  selected = utils::head(SeuratObject::VariableFeatures(o), 1L))
      shiny::updateSelectInput(session, "group_by", choices = meta$cat, selected = .xp_default_group(meta$cat))
      shiny::updateSelectInput(session, "comp_stack", choices = meta$cat, selected = .xp_default_group(meta$cat))
      shiny::updateSelectInput(session, "comp_across", choices = meta$cat,
                               selected = if ("orig.ident" %in% meta$cat) "orig.ident" else meta$cat[[1L]])
      shiny::updateTextInput(session, "dot_genes",
                             value = paste(utils::head(SeuratObject::VariableFeatures(o), 8L), collapse = ", "))
      shiny::updateSelectInput(session, "sig_set", choices = c("(custom)" = "", names(gene_sets_cfg)))
    })

    panel <- shiny::reactive({
      shiny::req(input$reduction)
      o <- obj()
      meta <- .xp_meta_columns(o)
      p <- sparks_panel(o, reduction = input$reduction, max_cells = max_cells,
                        meta_cols = c(meta$cat, meta$num), max_genes = 0L)
      p
    })
    # Panel plus everything created in this session (categories, re-clusters, scores).
    panel_view <- shiny::reactive({
      p <- panel()
      if (length(rv$categories)) {
        p <- sparks_panel_add_meta(p, "Selection categories",
                                   sparks_categories_to_factor(rv$categories, p$cells))
      }
      for (nm in names(rv$extra_meta)) p <- sparks_panel_add_meta(p, nm, rv$extra_meta[[nm]])
      p
    })

    gene_values <- function(p, gene) {
      o <- obj()
      if (!gene %in% rownames(o)) return(NULL)
      v <- .xp_gene_vector(o, gene, p$cells)
      v
    }

    # Colour-by choices follow session-created columns.
    shiny::observe({
      p <- panel_view()
      extra <- setdiff(names(p$meta), .xp_meta_columns(obj())$all)
      if (length(extra)) {
        o <- obj(); meta <- .xp_meta_columns(o)
        shiny::updateSelectizeInput(session, "color_by", server = TRUE,
          choices = c(extra, meta$cat, meta$num, rownames(o)), selected = shiny::isolate(input$color_by))
      }
    })

    # ---- Embedding -------------------------------------------------------------------------
    output$embedding <- plotly::renderPlotly({
      p <- panel_view()
      key <- input$color_by
      shiny::req(key)
      values <- NULL
      if (!key %in% names(p$meta)) values <- gene_values(p, key)
      fig <- sparks_fig_embedding(p, color_by = key, selected = rv$sel, point_size = input$point_size,
                                  opacity = input$opacity, source = src, values = values,
                                  value_label = key)
      shiny::isolate(rv$emb_ready <- TRUE)
      fig
    })
    # plotly warns when event_data() is first read before the plot that registers the event
    # has been drawn (the embedding is only drawn once the colour-by choices exist). So the
    # selection observers are attached after the first draw.
    shiny::observeEvent(rv$emb_ready, {
      if (isTRUE(rv$events_attached)) return()
      rv$events_attached <- TRUE
      shiny::observeEvent(plotly::event_data("plotly_selected", source = src), {
        ev <- plotly::event_data("plotly_selected", source = src)
        if (is.data.frame(ev) && nrow(ev) && "customdata" %in% names(ev)) {
          rv$sel <- unique(as.character(ev$customdata))
        }
      })
      shiny::observeEvent(plotly::event_data("plotly_deselect", source = src), rv$sel <- character())
    })

    # ---- Selection & categories -----------------------------------------------------------------
    output$sel_summary <- shiny::renderText({
      n <- length(rv$sel)
      if (!n) "No cells selected" else sprintf("%s cells selected", format(n, big.mark = ","))
    })
    shiny::observeEvent(input$clear_sel, rv$sel <- character())
    shiny::observeEvent(input$clear_cats, { rv$categories <- list() })
    shiny::observeEvent(input$save_cat, {
      nm <- trimws(input$cat_name %||% "")
      if (!length(rv$sel)) return(.xp_flash(session, "Select some cells first.", "warning"))
      if (!nzchar(nm)) nm <- paste("Category", length(rv$categories) + 1L)
      rv$categories[[nm]] <- rv$sel
      shiny::updateTextInput(session, "cat_name", value = "")
    })
    output$cats_ui <- shiny::renderUI({
      cats <- rv$categories
      if (!length(cats)) return(NULL)
      shiny::div(class = "mt-2", lapply(names(cats), function(nm) {
        shiny::span(class = "xp-cat", nm, shiny::span(class = "text-muted", length(cats[[nm]])))
      }))
    })
    shiny::observe({
      ch <- c("All other cells" = "rest", stats::setNames(names(rv$categories), names(rv$categories)))
      shiny::updateSelectInput(session, "de_vs", choices = ch, selected = shiny::isolate(input$de_vs))
    })
    shiny::observeEvent(input$import, {
      f <- input$import
      cells <- tryCatch(sparks_read_barcodes(f$datapath), error = function(e) character())
      .xp_set_selection(rv, session, cells, colnames(obj()))
    })
    shiny::observeEvent(input$apply_paste, {
      cells <- unique(trimws(unlist(strsplit(input$paste %||% "", "[\r\n,;\t ]+"))))
      .xp_set_selection(rv, session, cells[nzchar(cells)], colnames(obj()))
    })

    # ---- Genes / composition -------------------------------------------------------------------------
    output$violin <- plotly::renderPlotly({
      shiny::req(input$gene, input$group_by)
      p <- panel_view()
      v <- gene_values(p, input$gene)
      shiny::req(!is.null(v), input$group_by %in% names(p$meta))
      sparks_fig_violin(p, input$gene, input$group_by, values = v)
    })
    output$dotplot <- plotly::renderPlotly({
      shiny::req(input$group_by)
      p <- panel_view()
      genes <- unique(trimws(unlist(strsplit(input$dot_genes %||% "", "[,; ]+"))))
      genes <- intersect(genes[nzchar(genes)], rownames(obj()))
      shiny::req(length(genes) > 0L, input$group_by %in% names(p$meta))
      p <- sparks_panel_add_genes(p, obj(), genes)
      sparks_fig_dotplot(p, genes, input$group_by)
    })
    output$composition <- plotly::renderPlotly({
      shiny::req(input$comp_stack, input$comp_across)
      p <- panel_view()
      shiny::req(input$comp_stack %in% names(p$meta), input$comp_across %in% names(p$meta))
      sparks_fig_composition(.xp_composition_frame(p, input$comp_stack, input$comp_across),
                             title = paste(input$comp_stack, "by", input$comp_across))
    })

    # ---- Analyses -----------------------------------------------------------------------------------
    run_analysis <- function(spec) {
      id <- sparks_interactive_submit(results_dir, spec, obj_path())
      rv$flash <- NULL
      .xp_start_job(rv, results_dir, obj_path(), use_bg, session)
      rv$last_id <- id
      shiny::updateSelectInput(session, "view_id", selected = id)
      bslib::nav_select("views", "Analyses", session = session)
      id
    }
    shiny::observeEvent(input$run_de, {
      if (!length(rv$sel)) return(.xp_flash(session, "Select cells first (lasso on the embedding, or import barcodes).", "warning"))
      b <- if (identical(input$de_vs, "rest") || !nzchar(input$de_vs %||% "")) NULL else rv$categories[[input$de_vs]]
      lab <- sprintf("DE: selection (%d) vs %s", length(rv$sel),
                     if (is.null(b)) "rest" else sprintf("%s (%d)", input$de_vs, length(b)))
      run_analysis(sparks_spec("de", cells_a = rv$sel, cells_b = b, test = input$de_test,
                               min_pct = input$de_min_pct, logfc_threshold = input$de_lfc, label = lab))
    })
    shiny::observeEvent(input$run_recluster, {
      cells <- if (length(rv$sel)) rv$sel else colnames(obj())
      run_analysis(sparks_spec("recluster", cells = cells, resolution = input$rc_res))
    })
    shiny::observeEvent(input$sig_set, {
      if (nzchar(input$sig_set %||% "")) {
        shiny::updateTextAreaInput(session, "sig_genes", value = paste(gene_sets_cfg[[input$sig_set]], collapse = ", "))
      }
    })
    shiny::observeEvent(input$run_sig, {
      genes <- unique(trimws(unlist(strsplit(input$sig_genes %||% "", "[,; \r\n\t]+"))))
      genes <- genes[nzchar(genes)]
      if (length(genes) < 2L) return(.xp_flash(session, "Enter at least two genes.", "warning"))
      nm <- if (nzchar(input$sig_set %||% "")) input$sig_set else "Custom signature"
      run_analysis(sparks_spec("signature", genes = genes, name = nm))
    })
    shiny::observeEvent(input$gmt, {
      rv$gmt <- tryCatch(sparks_read_gmt(input$gmt$datapath), error = function(e) NULL)
    })
    shiny::observeEvent(input$run_ora, {
      st <- shiny::isolate(history())
      id <- input$ora_from
      if (!nzchar(id %||% "")) return(.xp_flash(session, "Pick a finished DE analysis as the gene list.", "warning"))
      res <- tryCatch(sparks_interactive_result(results_dir, id), error = function(e) NULL)
      if (is.null(res) || is.null(res$table)) return(.xp_flash(session, "That analysis has no gene table.", "warning"))
      tab <- res$table
      up <- tab$gene[!is.na(tab$p_val_adj) & tab$p_val_adj < 0.05 & tab$avg_log2FC > 0]
      sets <- c(gene_sets_cfg, rv$gmt %||% list())
      if (!length(sets)) return(.xp_flash(session, "No gene sets: add gene_signatures to the config or upload a GMT file.", "warning"))
      run_analysis(sparks_spec("ora", genes = up, universe = rownames(obj()), gene_sets = sets,
                               label = sprintf("ORA of %d up-regulated genes from %s", length(up), id)))
    })

    history <- shiny::reactive({
      rv$tick; input$refresh
      sparks_interactive_status(results_dir)
    })
    shiny::observe({
      st <- history()
      st <- st[st$object == obj_path() | !nzchar(st$object), , drop = FALSE]
      ch <- stats::setNames(st$id, sprintf("%s - %s [%s]", st$created, st$label, st$status))
      cur <- shiny::isolate(input$view_id)
      pick <- if (!is.null(rv$last_id) && rv$last_id %in% st$id) rv$last_id
              else if (!is.null(cur) && cur %in% st$id) cur else if (nrow(st)) st$id[[1L]] else character()
      shiny::updateSelectInput(session, "view_id", choices = ch, selected = pick)
      de <- st[st$type == "de" & st$status == "built", , drop = FALSE]
      shiny::updateSelectInput(session, "ora_from",
        choices = c("Select an analysis in History" = "", stats::setNames(de$id, de$label)))
    })
    output$history <- shiny::renderTable({
      st <- history()
      if (!nrow(st)) return(data.frame(Message = "No analyses yet. Select cells and use the Analyse panel."))
      data.frame(When = st$created, Analysis = st$label, Type = st$type, Status = st$status,
                 Seconds = ifelse(is.na(st$seconds), "", format(round(st$seconds, 1))),
                 stringsAsFactors = FALSE)
    }, striped = TRUE, width = "100%")
    output$job_status <- shiny::renderUI({
      rv$tick
      if (!is.null(rv$job)) return(shiny::div(class = "alert alert-info py-1", "Running analysis..."))
      if (!is.null(rv$flash)) return(shiny::div(class = "alert alert-warning py-1", rv$flash))
      NULL
    })
    shiny::observeEvent(input$delete, {
      if (nzchar(input$view_id %||% "")) {
        sparks_interactive_delete(results_dir, input$view_id)
        rv$tick <- rv$tick + 1L
      }
    })

    # Poll a background job.
    shiny::observe({
      if (is.null(rv$job)) return()
      shiny::invalidateLater(600, session)
      if (!rv$job$is_alive()) {
        rv$job <- NULL; rv$tick <- rv$tick + 1L
        if (any(shiny::isolate(history())$status == "pending")) {
          .xp_start_job(rv, results_dir, shiny::isolate(obj_path()), use_bg, session)
        }
      }
    })

    result <- shiny::reactive({
      rv$tick
      id <- input$view_id
      shiny::req(nzchar(id %||% ""))
      st <- history()
      row <- st[st$id == id, , drop = FALSE]
      shiny::req(nrow(row) == 1L)
      if (row$status == "error") return(list(error = row$error))
      if (row$status != "built") return(list(pending = TRUE))
      sparks_interactive_result(results_dir, id)
    })
    output$result_ui <- shiny::renderUI({
      r <- result()
      if (!is.null(r$error)) return(shiny::div(class = "alert alert-danger", r$error))
      if (isTRUE(r$pending)) return(shiny::div(class = "text-muted", "Waiting for the analysis to finish..."))
      switch(
        r$type,
        de = shiny::tagList(shiny::p(sprintf("%s cells vs %s cells (assay %s)", r$n_a, r$n_b, r$assay)),
                            plotly::plotlyOutput(ns("res_volcano"), height = "460px"),
                            .xp_table_output(ns("res_table")),
                            shiny::downloadButton(ns("dl_table"), "Download CSV")),
        recluster = shiny::tagList(shiny::p(sprintf("%d clusters at resolution %s", r$n_clusters, r$resolution)),
                                   plotly::plotlyOutput(ns("res_recluster"), height = "460px"),
                                   shiny::actionButton(ns("apply_recluster"), "Add clusters as a colour-by column"),
                                   if (!is.null(r$markers)) .xp_table_output(ns("res_table"))),
        signature = shiny::tagList(shiny::p(sprintf("%s: %d genes scored (%s)%s", r$name, length(r$genes), r$method,
                                                    if (length(r$missing)) paste("; missing:", paste(r$missing, collapse = ", ")) else "")),
                                   shiny::actionButton(ns("apply_signature"), "Add score as a colour-by column")),
        ora = shiny::tagList(shiny::p(sprintf("%d genes tested against %d universe genes", r$n_genes, r$universe)),
                             .xp_table_output(ns("res_table")),
                             shiny::downloadButton(ns("dl_table"), "Download CSV")),
        NULL
      )
    })
    result_table <- shiny::reactive({
      r <- result()
      shiny::req(is.null(r$error), is.null(r$pending))
      if (!is.null(r$table)) r$table else r$markers
    })
    output$res_volcano <- plotly::renderPlotly({
      sparks_fig_volcano(result_table(), title = NULL)
    })
    output$res_recluster <- plotly::renderPlotly({
      r <- result(); shiny::req(r$type == "recluster")
      sparks_fig_recluster(r)
    })
    .xp_render_table(output, "res_table", shiny::reactive({
      t <- result_table(); utils::head(t, 500L)
    }))
    output$dl_table <- shiny::downloadHandler(
      filename = function() paste0("sparks_", input$view_id, ".csv"),
      content = function(file) utils::write.csv(result_table(), file, row.names = FALSE))
    shiny::observeEvent(input$apply_recluster, {
      r <- result(); shiny::req(r$type == "recluster")
      rv$extra_meta[[paste0("recluster_", substr(input$view_id, 1L, 6L))]] <- r$clusters
      .xp_flash(session, "Added. Pick it in 'Colour by' on the Embedding tab.", "success")
    })
    shiny::observeEvent(input$apply_signature, {
      r <- result(); shiny::req(r$type == "signature")
      rv$extra_meta[[paste0("score_", gsub("[^A-Za-z0-9]+", "_", r$name))]] <- r$scores
      .xp_flash(session, "Added. Pick it in 'Colour by' on the Embedding tab.", "success")
    })

    # ---- Report ----------------------------------------------------------------------------------------
    shiny::observeEvent(input$build_report, {
      out <- file.path(results_dir, "SPARKS_all_results.html")
      shiny::withProgress(message = "Building the unified report", value = 0.3, {
        res <- tryCatch(sparks_report(results_dir, output_file = out, quiet = TRUE),
                        error = function(e) e)
      })
      if (inherits(res, "error")) {
        rv$flash <- paste("Report failed:", conditionMessage(res))
      } else {
        rv$report <- out
        rv$flash <- NULL
      }
      rv$tick <- rv$tick + 1L
    })
    output$report_link <- shiny::renderUI({
      if (is.null(rv$report)) return(NULL)
      shiny::downloadLink(ns("dl_report"), "Download SPARKS_all_results.html")
    })
    output$dl_report <- shiny::downloadHandler(
      filename = function() "SPARKS_all_results.html",
      content = function(file) file.copy(rv$report, file, overwrite = TRUE))

    # Exposed for tests.
    list(
      selection = shiny::reactive(rv$sel),
      categories = shiny::reactive(rv$categories),
      history = history,
      panel = panel_view,
      extra_meta = shiny::reactive(rv$extra_meta),
      result = result,
      job_running = shiny::reactive(!is.null(rv$job))
    )
  })
}

# ---- Helpers -----------------------------------------------------------------------------------

.xp_can_background <- function() {
  requireNamespace("callr", quietly = TRUE) &&
    nzchar(system.file(package = "SPARKS")) &&
    !(requireNamespace("pkgload", quietly = TRUE) && isTRUE(tryCatch(pkgload::is_dev_package("SPARKS"), error = function(e) FALSE)))
}

.xp_flash <- function(session, text, type = "warning") {
  shiny::showNotification(text, type = switch(type, success = "message", type), duration = 5, session = session)
}

# Which metadata columns are worth offering in the UI.
.xp_meta_columns <- function(object) {
  m <- object[[]]
  is_cat <- vapply(m, function(v) {
    n <- length(unique(v))
    (is.character(v) || is.factor(v) || is.logical(v) || (is.numeric(v) && n <= 25L)) && n > 1L && n <= 200L
  }, logical(1))
  is_num <- vapply(m, is.numeric, logical(1)) & !is_cat
  list(cat = colnames(m)[is_cat], num = colnames(m)[is_num], all = colnames(m))
}

.xp_default_group <- function(cats) {
  pref <- c("cell_type", "celltype", "SingleR_main", "seurat_clusters", "orig.ident")
  hit <- pref[pref %in% cats]
  if (length(hit)) hit[[1L]] else cats[[1L]]
}

.xp_gene_vector <- function(object, gene, cells) {
  assay <- .analysis_assay(object)
  mat <- tryCatch(SeuratObject::LayerData(object, assay = assay, layer = "data"),
                  error = function(e) NULL)
  if (is.null(mat) || !gene %in% rownames(mat)) {
    mat <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  }
  as.numeric(mat[gene, cells])
}

.xp_set_selection <- function(rv, session, cells, all_cells) {
  found <- intersect(cells, all_cells)
  if (!length(found)) {
    .xp_flash(session, sprintf("None of the %d barcodes are in this object.", length(cells)), "warning")
    return(invisible())
  }
  rv$sel <- found
  if (length(found) < length(cells)) {
    .xp_flash(session, sprintf("%d of %d barcodes were not found in this object.", length(cells) - length(found), length(cells)), "warning")
  }
}

.xp_composition_frame <- function(panel, stack, across) {
  a <- .panel_meta_values(panel, stack)
  b <- .panel_meta_values(panel, across)
  tab <- table(a, b)
  df <- as.data.frame(prop.table(tab, margin = 2L), stringsAsFactors = FALSE)
  names(df) <- c("Var1", "Var2", "Freq")
  df$Count <- as.integer(tab)
  df
}

.xp_start_job <- function(rv, results_dir, object_path, use_bg, session) {
  if (!is.null(rv$job)) return(invisible())   # pending specs are picked up when it ends
  if (use_bg) {
    rv$job <- callr::r_bg(function(rd, op) SPARKS::sparks_interactive_run(rd, op, background = FALSE),
                          args = list(rd = results_dir, op = object_path), supervise = TRUE)
  } else {
    shiny::withProgress(message = "Running analysis", value = 0.5, {
      tryCatch(sparks_interactive_run(results_dir, object_path),
               error = function(e) rv$flash <- conditionMessage(e))
    })
    rv$tick <- rv$tick + 1L
  }
  invisible()
}

.xp_table_output <- function(id) {
  if (requireNamespace("DT", quietly = TRUE)) DT::DTOutput(id) else shiny::tableOutput(id)
}

.xp_render_table <- function(output, id, data) {
  if (requireNamespace("DT", quietly = TRUE)) {
    output[[id]] <- DT::renderDT(DT::datatable(data(), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) |>
                                    .xp_round_dt())
  } else {
    output[[id]] <- shiny::renderTable(data(), striped = TRUE)
  }
}

.xp_round_dt <- function(dt) {
  d <- dt$x$data
  num <- names(d)[vapply(d, is.numeric, logical(1))]
  if (length(num)) dt <- DT::formatSignif(dt, num, digits = 3)
  dt
}

#' Embedding of a re-clustering result
#' @param analysis A \code{sparks_analysis} of type \code{"recluster"}.
#' @param mode Colour mode.
#' @return plotly object.
#' @export
sparks_fig_recluster <- function(analysis, mode = "light") {
  .require_plotly()
  emb <- analysis$embedding
  cl <- analysis$clusters[rownames(emb)]
  lv <- sort(unique(cl))
  lv <- lv[order(suppressWarnings(as.numeric(lv)), lv)]
  cols <- .level_colors(lv, mode)
  p <- plotly::plot_ly()
  for (k in lv) {
    i <- which(cl == k)
    p <- plotly::add_trace(p, type = if (length(cl) > 3000L) "scattergl" else "scatter", mode = "markers",
                           x = emb[i, 1L], y = emb[i, 2L], name = paste("Cluster", k),
                           marker = list(size = 5, color = cols[[k]], opacity = 0.85),
                           text = rownames(emb)[i], hoverinfo = "text+name")
  }
  .fig_config(.fig_base(p, NULL, "UMAP 1", "UMAP 2", TRUE, mode), "recluster")
}
