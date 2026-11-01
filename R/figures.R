# Shared, data-in figure builders.
#
# Every function takes plain data frames / a `sparks_panel` and returns a
# plotly object (or NULL when there is nothing to draw). The same builders feed
# the Shiny explorer and the figures embedded in the unified HTML report, so a
# plot looks the same in both places (the alevinQC "data in, figure out" pattern).
#
# Colours follow a validated categorical / sequential / diverging scheme. Figures
# are built with the *light* hex values; the report and the app swap them for the
# dark-mode steps through `.sparks_theme_map()` when the theme changes.

# ---- Palette -------------------------------------------------------------------

.SPARKS_CAT_LIGHT <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                       "#e87ba4", "#008300", "#6250d6", "#e34948")
.SPARKS_CAT_DARK  <- c("#3987e5", "#d95926", "#199e70", "#c98500",
                       "#d55181", "#008300", "#9085e9", "#e66767")
.SPARKS_SEQ  <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
                  "#184f95", "#0d366b")
.SPARKS_DIV  <- c("#e34948", "#f0a3a2", "#f0efec", "#9ec5f4", "#2a78d6")
.SPARKS_CHROME <- list(
  light = list(surface = "#fcfcfb", text = "#0b0b0b", secondary = "#52514e",
               muted = "#898781", grid = "#e1e0d9", axis = "#c3c2b7",
               mid = "#f0efec"),
  dark  = list(surface = "#1a1a19", text = "#f0efec", secondary = "#c3c2b7",
               muted = "#898781", grid = "#2c2c2a", axis = "#383835",
               mid = "#383835")
)

#' Categorical colours for \code{n} levels
#'
#' The first eight levels use the validated categorical palette in fixed order
#' (never cycled by rank). Beyond eight, colours are recycled with progressively
#' lighter tints -- a deliberate, documented deviation, because clusterings
#' routinely have more than eight groups; figures with that many groups always
#' carry direct labels so identity never rests on colour alone.
#'
#' @param n Number of levels.
#' @param mode \code{"light"} or \code{"dark"}.
#' @return Character vector of \code{n} hex colours.
#' @export
sparks_colors <- function(n, mode = c("light", "dark")) {
  mode <- match.arg(mode)
  base <- if (mode == "light") .SPARKS_CAT_LIGHT else .SPARKS_CAT_DARK
  if (n <= 0L) return(character())
  idx <- (seq_len(n) - 1L) %% 8L + 1L
  round_no <- (seq_len(n) - 1L) %/% 8L
  vapply(seq_len(n), function(i) .tint(base[[idx[[i]]]], c(0, 0.38, 0.62)[min(round_no[[i]], 2L) + 1L],
                                       toward = if (mode == "light") "#ffffff" else "#000000"),
         character(1))
}

# Mix `col` with `toward` by `amount` (0 = col, 1 = toward).
.tint <- function(col, amount, toward = "#ffffff") {
  if (amount <= 0) return(col)
  a <- grDevices::col2rgb(col)[, 1L]
  b <- grDevices::col2rgb(toward)[, 1L]
  m <- round(a + (b - a) * amount)
  grDevices::rgb(m[[1L]], m[[2L]], m[[3L]], maxColorValue = 255)
}

#' Light -> dark colour lookup used when the theme changes
#' @return Named character vector (names = light hex, values = dark hex).
#' @keywords internal
.sparks_theme_map <- function() {
  light <- c(.SPARKS_CAT_LIGHT, unlist(.SPARKS_CHROME$light, use.names = FALSE))
  dark  <- c(.SPARKS_CAT_DARK,  unlist(.SPARKS_CHROME$dark,  use.names = FALSE))
  keep <- !duplicated(light)
  stats::setNames(dark[keep], light[keep])
}

# Named level -> colour assignment, stable across plots (colour follows the entity).
.level_colors <- function(levels, mode = "light") {
  levels <- unique(as.character(levels))
  stats::setNames(sparks_colors(length(levels), mode), levels)
}

.fig_font <- "system-ui, -apple-system, 'Segoe UI', sans-serif"

.require_plotly <- function() {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Interactive figures require package 'plotly'. ",
         "Install it with install.packages('plotly').", call. = FALSE)
  }
}

.fig_base <- function(p, title = NULL, xlab = NULL, ylab = NULL,
                      legend = TRUE, mode = "light", margin = NULL) {
  ch <- .SPARKS_CHROME[[mode]]
  axis <- function(label) list(
    title = list(text = label %||% "", font = list(color = ch$secondary)),
    gridcolor = ch$grid, zerolinecolor = ch$axis, linecolor = ch$axis,
    tickfont = list(color = ch$muted), automargin = TRUE
  )
  plotly::layout(
    p,
    title = if (!is.null(title)) list(text = title, x = 0, xanchor = "left",
                                      font = list(size = 14, color = ch$text)),
    paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
    font = list(family = .fig_font, size = 12, color = ch$secondary),
    xaxis = axis(xlab), yaxis = axis(ylab),
    showlegend = legend,
    legend = list(orientation = "h", y = -0.22, font = list(color = ch$secondary)),
    margin = margin %||% list(l = 56, r = 16, t = if (is.null(title)) 16 else 44, b = 48),
    hoverlabel = list(font = list(family = .fig_font))
  )
}

.fig_config <- function(p, filename = "sparks_figure") {
  plotly::config(p, displaylogo = FALSE, responsive = TRUE,
                 modeBarButtonsToRemove = c("select2d", "autoScale2d"),
                 toImageButtonOptions = list(format = "png", filename = filename,
                                             scale = 2))
}

# ---- QC ---------------------------------------------------------------------------

#' Cells per sample before / after filtering
#' @param stages Data frame from \code{bundle$qc$stages} (columns: stage, label,
#'   sample, cells).
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_qc_stages <- function(stages, mode = "light") {
  .require_plotly()
  if (!is.data.frame(stages) || !nrow(stages)) return(NULL)
  d <- stages[stages$stage %in% c("before_filtering", "after_filtering"), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  d$sample <- ifelse(is.na(d$sample), d$label, d$sample)
  samples <- unique(d$sample)
  cols <- c(before_filtering = sparks_colors(2L, mode)[[1L]],
            after_filtering = sparks_colors(2L, mode)[[2L]])
  names_ <- c(before_filtering = "Before filtering", after_filtering = "After filtering")
  p <- plotly::plot_ly()
  for (st in c("before_filtering", "after_filtering")) {
    s <- d[d$stage == st, , drop = FALSE]
    if (!nrow(s)) next
    s <- s[match(samples, s$sample), , drop = FALSE]
    p <- plotly::add_trace(p, type = "bar", x = samples, y = s$cells,
                           name = names_[[st]], marker = list(color = cols[[st]]),
                           hovertemplate = "%{x}<br>%{y} cells<extra></extra>")
  }
  p <- plotly::layout(p, barmode = "group", bargap = 0.35)
  .fig_config(.fig_base(p, "Cells per sample", NULL, "Cells", TRUE, mode), "qc_cells")
}

#' Cells per sample
#' @param counts Data frame with \code{sample} and \code{cells}.
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_cells_per_sample <- function(counts, mode = "light") {
  .require_plotly()
  if (!is.data.frame(counts) || !all(c("sample", "cells") %in% names(counts)) ||
      !nrow(counts)) return(NULL)
  p <- plotly::plot_ly(counts, type = "bar", x = ~as.character(sample), y = ~cells,
                       marker = list(color = sparks_colors(1L, mode)),
                       text = ~cells, textposition = "outside", cliponaxis = FALSE,
                       hovertemplate = "%{x}<br>%{y} cells<extra></extra>")
  p <- plotly::layout(p, xaxis = list(type = "category"), bargap = 0.4)
  .fig_config(.fig_base(p, NULL, NULL, "Cells", FALSE, mode,
                        margin = list(l = 56, r = 16, t = 24, b = 48)), "cells_per_sample")
}

#' Singlets and doublets per sample
#' @param doublets Data frame from \code{bundle$qc$doublets}.
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_doublets <- function(doublets, mode = "light") {
  .require_plotly()
  if (!is.data.frame(doublets) || !nrow(doublets)) return(NULL)
  cls <- unique(doublets$class)
  cls <- c(intersect("singlet", cls), setdiff(cls, "singlet"))
  cols <- stats::setNames(sparks_colors(length(cls), mode), cls)
  samples <- unique(doublets$sample)
  ch <- .SPARKS_CHROME[[mode]]
  p <- plotly::plot_ly()
  for (k in cls) {
    d <- doublets[doublets$class == k, , drop = FALSE]
    d <- d[match(samples, d$sample), , drop = FALSE]
    p <- plotly::add_trace(
      p, type = "bar", x = samples, y = ifelse(is.na(d$cells), 0L, d$cells), name = k,
      marker = list(color = cols[[k]], line = list(color = ch$surface, width = 2)),
      hovertemplate = paste0("%{x}<br>", k, ": %{y}<extra></extra>"))
  }
  p <- plotly::layout(p, barmode = "stack", bargap = 0.35)
  .fig_config(.fig_base(p, "scDblFinder calls per sample", NULL, "Cells", TRUE, mode),
              "qc_doublets")
}

# ---- DEG ------------------------------------------------------------------------

#' Significant genes per group
#' @param counts A \code{DEG_Counts_*} table (\code{Group}, \code{DEG_Count}).
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_deg_counts <- function(counts, mode = "light") {
  .require_plotly()
  if (!is.data.frame(counts) || !all(c("Group", "DEG_Count") %in% names(counts)) ||
      !nrow(counts)) return(NULL)
  counts <- counts[order(counts$DEG_Count), , drop = FALSE]
  col <- sparks_colors(1L, mode)
  p <- plotly::plot_ly(counts, type = "bar", orientation = "h",
                       x = ~DEG_Count, y = ~as.character(Group),
                       marker = list(color = col),
                       text = ~DEG_Count, textposition = "outside", cliponaxis = FALSE,
                       hovertemplate = "%{y}: %{x} DEGs<extra></extra>")
  p <- plotly::layout(p, yaxis = list(type = "category", categoryorder = "array",
                                      categoryarray = as.character(counts$Group)),
                      bargap = 0.4)
  .fig_config(.fig_base(p, "Significant genes per group", "Genes", NULL, FALSE, mode,
                        margin = list(l = 72, r = 40, t = 44, b = 40)), "deg_counts")
}

#' Volcano plot
#'
#' Up / down use the two poles of the diverging pair; non-significant genes are
#' neutral grey. The \code{label_top} most significant genes per direction are
#' labelled directly.
#'
#' @param deg Data frame with fold-change, adjusted p-value and gene columns.
#' @param lfc_col,p_col,gene_col Column names. If \code{lfc_col} is missing,
#'   \code{logFC} is tried; if \code{p_col} is missing, \code{padj} / \code{p_val} are.
#' @param lfc_cut,p_cut Significance thresholds.
#' @param label_top Genes to label per direction.
#' @param title Plot title.
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_volcano <- function(deg, lfc_col = "avg_log2FC", p_col = "p_val_adj",
                               gene_col = "gene", lfc_cut = 0.25, p_cut = 0.05,
                               label_top = 8L, title = NULL, mode = "light") {
  .require_plotly()
  if (!is.data.frame(deg) || !nrow(deg)) return(NULL)
  pick <- function(col, alts) if (col %in% names(deg)) col else
    intersect(alts, names(deg))[1L]
  lfc_col <- pick(lfc_col, c("logFC", "avg_logFC", "log2FoldChange"))
  p_col <- pick(p_col, c("padj", "FDR", "adj.P.Val", "p_val", "PValue"))
  gene_col <- pick(gene_col, c("gene", "Gene", "symbol", "feature"))
  if (is.na(lfc_col) || is.na(p_col)) return(NULL)
  genes <- if (is.na(gene_col)) as.character(seq_len(nrow(deg))) else as.character(deg[[gene_col]])
  lfc <- suppressWarnings(as.numeric(deg[[lfc_col]]))
  pv <- suppressWarnings(as.numeric(deg[[p_col]]))
  ok <- is.finite(lfc) & is.finite(pv)
  d <- data.frame(gene = genes[ok], lfc = lfc[ok], p = pv[ok], stringsAsFactors = FALSE)
  if (!nrow(d)) return(NULL)
  d$nlp <- -log10(pmax(d$p, 1e-300))
  d$class <- ifelse(d$p < p_cut & d$lfc >= lfc_cut, "Up",
                    ifelse(d$p < p_cut & d$lfc <= -lfc_cut, "Down", "Not significant"))
  cols <- c(Up = .SPARKS_DIV[[1L]], Down = .SPARKS_DIV[[5L]],
            "Not significant" = .SPARKS_CHROME[[mode]]$muted)
  trace_type <- if (nrow(d) > 4000L) "scattergl" else "scatter"
  p <- plotly::plot_ly()
  for (k in c("Not significant", "Down", "Up")) {
    s <- d[d$class == k, , drop = FALSE]
    if (!nrow(s)) next
    p <- plotly::add_trace(
      p, type = trace_type, mode = "markers", x = s$lfc, y = s$nlp, name = paste0(k, " (", nrow(s), ")"),
      text = s$gene, marker = list(color = cols[[k]], size = 6, opacity = if (k == "Not significant") 0.45 else 0.85),
      hovertemplate = "%{text}<br>log2FC %{x:.2f}<br>-log10 p %{y:.1f}<extra></extra>")
  }
  lab <- do.call(rbind, lapply(c("Up", "Down"), function(k) {
    s <- d[d$class == k, , drop = FALSE]
    utils::head(s[order(-s$nlp), , drop = FALSE], label_top)
  }))
  annotations <- NULL
  if (!is.null(lab) && nrow(lab)) {
    annotations <- lapply(seq_len(nrow(lab)), function(i) list(
      x = lab$lfc[[i]], y = lab$nlp[[i]], text = lab$gene[[i]], showarrow = FALSE,
      xanchor = if (lab$lfc[[i]] > 0) "left" else "right", yanchor = "bottom",
      font = list(size = 10, color = .SPARKS_CHROME[[mode]]$secondary)))
  }
  p <- plotly::layout(
    p, annotations = annotations,
    shapes = list(
      list(type = "line", x0 = -lfc_cut, x1 = -lfc_cut, y0 = 0, y1 = 1, yref = "paper",
           line = list(color = .SPARKS_CHROME[[mode]]$axis, width = 1, dash = "dot")),
      list(type = "line", x0 = lfc_cut, x1 = lfc_cut, y0 = 0, y1 = 1, yref = "paper",
           line = list(color = .SPARKS_CHROME[[mode]]$axis, width = 1, dash = "dot")),
      list(type = "line", x0 = 0, x1 = 1, xref = "paper", y0 = -log10(p_cut), y1 = -log10(p_cut),
           line = list(color = .SPARKS_CHROME[[mode]]$axis, width = 1, dash = "dot"))))
  .fig_config(.fig_base(p, title, "log2 fold change", "-log10 adjusted p", TRUE, mode), "volcano")
}

# ---- Composition -------------------------------------------------------------------

#' Cell-type composition per sample / condition
#'
#' @param prop A \code{Proportions_*} table (\code{Var1} = group, \code{Var2} =
#'   sample/condition, \code{Freq} = proportion, \code{Count}).
#' @param mode Colour mode.
#' @param title Plot title.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_composition <- function(prop, mode = "light", title = NULL) {
  .require_plotly()
  need <- c("Var1", "Var2", "Freq")
  if (!is.data.frame(prop) || !all(need %in% names(prop)) || !nrow(prop)) return(NULL)
  prop$Var1 <- as.character(prop$Var1)
  prop$Var2 <- as.character(prop$Var2)
  groups <- unique(prop$Var1)
  cols <- .level_colors(groups, mode)
  ch <- .SPARKS_CHROME[[mode]]
  x_levels <- unique(prop$Var2)
  p <- plotly::plot_ly()
  for (g in groups) {
    s <- prop[prop$Var1 == g, , drop = FALSE]
    s <- s[match(x_levels, s$Var2), , drop = FALSE]
    cnt <- if ("Count" %in% names(s)) s$Count else rep(NA, nrow(s))
    p <- plotly::add_trace(
      p, type = "bar", x = x_levels, y = ifelse(is.na(s$Freq), 0, s$Freq * 100),
      name = g, customdata = cnt,
      marker = list(color = cols[[g]], line = list(color = ch$surface, width = 2)),
      hovertemplate = paste0("%{x}<br>", g, ": %{y:.1f}% (%{customdata} cells)<extra></extra>"))
  }
  p <- plotly::layout(p, barmode = "stack", bargap = 0.3,
                      yaxis = list(range = c(0, 100), ticksuffix = "%"))
  .fig_config(.fig_base(p, title, NULL, "Share of cells", TRUE, mode), "composition")
}

# ---- Enrichment ---------------------------------------------------------------------

#' GSEA dot plot: pathways x groups
#'
#' Colour is the normalised enrichment score on the diverging scale (zero at the
#' neutral midpoint); size is -log10 adjusted p.
#'
#' @param gsea Stacked GSEA tables (\code{pathway}, \code{group}, \code{NES},
#'   \code{padj}).
#' @param top_n Pathways kept (ranked by their strongest |NES|).
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_gsea <- function(gsea, top_n = 20L, mode = "light") {
  .require_plotly()
  need <- c("pathway", "group", "NES", "padj")
  if (!is.data.frame(gsea) || !all(need %in% names(gsea)) || !nrow(gsea)) return(NULL)
  gsea$NES <- suppressWarnings(as.numeric(gsea$NES))
  gsea$padj <- suppressWarnings(as.numeric(gsea$padj))
  gsea <- gsea[is.finite(gsea$NES) & is.finite(gsea$padj), , drop = FALSE]
  if (!nrow(gsea)) return(NULL)
  strength <- tapply(abs(gsea$NES), gsea$pathway, max)
  keep <- names(sort(strength, decreasing = TRUE))[seq_len(min(top_n, length(strength)))]
  d <- gsea[gsea$pathway %in% keep, , drop = FALSE]
  d$pathway <- factor(d$pathway, levels = rev(keep))
  d$nlp <- -log10(pmax(d$padj, 1e-300))
  lim <- max(abs(d$NES), 1)
  p <- plotly::plot_ly(
    d, type = "scatter", mode = "markers", x = ~as.character(group), y = ~as.character(pathway),
    marker = list(
      size = pmin(8 + 6 * d$nlp, 30), sizemode = "diameter",
      color = d$NES, cmin = -lim, cmax = lim,
      colorscale = list(c(0, .SPARKS_DIV[[5L]]), c(0.5, .SPARKS_CHROME[[mode]]$mid), c(1, .SPARKS_DIV[[1L]])),
      colorbar = list(title = list(text = "NES"), thickness = 10, len = 0.6),
      line = list(color = .SPARKS_CHROME[[mode]]$surface, width = 1)),
    text = ~paste0(pathway, "<br>", group, "<br>NES ", round(NES, 2), "<br>padj ", signif(padj, 2)),
    hoverinfo = "text")
  p <- plotly::layout(p, xaxis = list(type = "category"),
                      yaxis = list(type = "category", categoryorder = "array",
                                   categoryarray = levels(d$pathway)))
  .fig_config(.fig_base(p, "Gene-set enrichment (NES)", NULL, NULL, FALSE, mode,
                        margin = list(l = 220, r = 16, t = 44, b = 48)), "gsea")
}

# ---- Heatmap -------------------------------------------------------------------------

#' Heatmap of a numeric matrix
#' @param mat Numeric matrix (rows x columns).
#' @param title Plot title.
#' @param diverging Use the diverging scale centred on zero (else sequential).
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_heatmap <- function(mat, title = NULL, diverging = FALSE, mode = "light") {
  .require_plotly()
  if (!is.matrix(mat) || !length(mat)) return(NULL)
  ch <- .SPARKS_CHROME[[mode]]
  scale <- if (diverging) {
    list(c(0, .SPARKS_DIV[[5L]]), c(0.25, .SPARKS_DIV[[4L]]), c(0.5, ch$mid),
         c(0.75, .SPARKS_DIV[[2L]]), c(1, .SPARKS_DIV[[1L]]))
  } else {
    n <- length(.SPARKS_SEQ)
    lapply(seq_len(n), function(i) list((i - 1) / (n - 1), .SPARKS_SEQ[[i]]))
  }
  lim <- if (diverging) max(abs(mat), na.rm = TRUE) else NULL
  p <- plotly::plot_ly(
    z = mat, x = colnames(mat), y = rownames(mat), type = "heatmap",
    colorscale = scale, zmin = if (diverging) -lim, zmax = if (diverging) lim,
    colorbar = list(thickness = 10, len = 0.6),
    xgap = 1, ygap = 1)
  p <- plotly::layout(p, yaxis = list(autorange = "reversed"))
  .fig_config(.fig_base(p, title, NULL, NULL, FALSE, mode,
                        margin = list(l = 90, r = 16, t = 44, b = 80)), "heatmap")
}

# ---- Panel figures (embedding, violin, dot plot) ---------------------------------------

#' Embedding scatter coloured by metadata or gene expression
#'
#' @param panel A [sparks_panel()].
#' @param color_by A metadata column, or a gene present in the panel.
#' @param selected Optional character vector of cell barcodes to highlight.
#' @param point_size,opacity Marker styling.
#' @param mode Colour mode.
#' @param source Plotly event source id (for lasso selection in Shiny).
#' @param values Optional numeric vector (one per panel cell) to colour by instead
#'   of \code{color_by} -- used for genes that are not stored in the panel.
#' @param value_label Legend title for \code{values}.
#' @return plotly object.
#' @export
sparks_fig_embedding <- function(panel, color_by = NULL, selected = NULL,
                                 point_size = 4, opacity = 0.8, mode = "light",
                                 source = "sparks_embedding", values = NULL,
                                 value_label = color_by) {
  .require_plotly()
  stopifnot(inherits(panel, "sparks_panel"))
  n <- length(panel$cells)
  ch <- .SPARKS_CHROME[[mode]]
  type <- if (n > 3000L) "scattergl" else "scatter"
  axis_names <- paste0(toupper(panel$reduction), c(" 1", " 2"))
  base <- plotly::plot_ly(source = source)
  is_gene <- is.null(values) && !is.null(color_by) && color_by %in% names(panel$genes) &&
    !color_by %in% names(panel$meta)
  is_num <- is.null(values) && !is.null(color_by) && color_by %in% names(panel$meta) &&
    identical(panel$meta[[color_by]]$type, "num")
  if (!is.null(values) || is_gene || is_num) {
    v <- if (!is.null(values)) as.numeric(values) else if (is_gene)
      .panel_gene_values(panel, color_by) else panel$meta[[color_by]]$values
    color_by <- value_label %||% color_by %||% "value"
    ord <- order(v, na.last = FALSE)
    n_seq <- length(.SPARKS_SEQ)
    p <- plotly::add_trace(
      base, type = type, mode = "markers", x = panel$x[ord], y = panel$y[ord],
      customdata = panel$cells[ord], text = signif(v[ord], 3L),
      marker = list(
        size = point_size, opacity = opacity, color = v[ord],
        colorscale = lapply(seq_len(n_seq), function(i) list((i - 1) / (n_seq - 1), .SPARKS_SEQ[[i]])),
        colorbar = list(title = list(text = color_by), thickness = 10, len = 0.6)),
      hovertemplate = paste0("%{customdata}<br>", color_by, ": %{text}<extra></extra>"),
      showlegend = FALSE)
  } else {
    groups <- if (!is.null(color_by) && color_by %in% names(panel$meta))
      .panel_meta_values(panel, color_by) else factor(rep("All cells", n))
    lv <- levels(groups)
    cols <- .level_colors(lv, mode)
    p <- base
    for (g in lv) {
      idx <- which(groups == g)
      p <- plotly::add_trace(
        p, type = type, mode = "markers", x = panel$x[idx], y = panel$y[idx], name = g,
        customdata = panel$cells[idx],
        marker = list(size = point_size, opacity = opacity, color = cols[[g]]),
        hovertemplate = paste0("%{customdata}<br>", g, "<extra></extra>"))
    }
    # Direct labels at group centroids (identity never rests on colour alone).
    if (length(lv) > 1L && length(lv) <= 40L) {
      ann <- lapply(lv, function(g) {
        idx <- which(groups == g)
        list(x = stats::median(panel$x[idx]), y = stats::median(panel$y[idx]),
             text = g, showarrow = FALSE,
             font = list(size = 11, color = ch$text),
             bgcolor = if (mode == "light") "rgba(252,252,251,0.7)" else "rgba(26,26,25,0.7)")
      })
      p <- plotly::layout(p, annotations = ann)
    }
  }
  if (length(selected)) {
    idx <- which(panel$cells %in% selected)
    if (length(idx)) {
      p <- plotly::add_trace(
        p, type = "scatter", mode = "markers", x = panel$x[idx], y = panel$y[idx],
        name = paste0("Selection (", length(idx), ")"), hoverinfo = "skip",
        marker = list(size = point_size + 3, color = "rgba(0,0,0,0)",
                      line = list(color = ch$text, width = 1.5)))
    }
  }
  p <- plotly::layout(p, dragmode = "lasso",
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, title = list(text = axis_names[[1L]])),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, title = list(text = axis_names[[2L]]),
                                   scaleanchor = "x"))
  p <- .fig_base(p, NULL, axis_names[[1L]], axis_names[[2L]], TRUE, mode,
                 margin = list(l = 48, r = 16, t = 8, b = 40))
  p <- plotly::layout(p, yaxis = list(showgrid = FALSE, zeroline = FALSE, scaleanchor = "x"),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE))
  p <- plotly::event_register(p, "plotly_selected")
  p <- plotly::event_register(p, "plotly_deselect")
  .fig_config(p, "embedding")
}

#' Expression of one gene across groups (violin)
#' @param panel A [sparks_panel()].
#' @param gene Gene in the panel.
#' @param group_by Categorical metadata column.
#' @param mode Colour mode.
#' @param values Optional expression vector (one per panel cell) for genes not stored in the panel.
#' @return plotly object.
#' @export
sparks_fig_violin <- function(panel, gene, group_by, mode = "light", values = NULL) {
  .require_plotly()
  stopifnot(inherits(panel, "sparks_panel"))
  v <- if (!is.null(values)) as.numeric(values) else .panel_gene_values(panel, gene)
  g <- .panel_meta_values(panel, group_by)
  lv <- levels(g)
  cols <- .level_colors(lv, mode)
  p <- plotly::plot_ly()
  for (k in lv) {
    vals <- v[g == k]
    if (!length(vals)) next
    p <- plotly::add_trace(
      p, type = "violin", y = vals, x = rep(k, length(vals)), name = k,
      points = FALSE, spanmode = "hard", scalemode = "width", width = 0.8,
      fillcolor = cols[[k]], line = list(color = cols[[k]], width = 1),
      meanline = list(visible = TRUE), opacity = 0.9,
      hoverinfo = "y+name")
  }
  p <- plotly::layout(p, xaxis = list(type = "category", categoryarray = lv))
  .fig_config(.fig_base(p, paste0(gene, " by ", group_by), group_by,
                        "Expression (log-normalised)", FALSE, mode), paste0("violin_", gene))
}

#' Per-group mean expression and detection rate of a gene set
#' @param panel A [sparks_panel()].
#' @param genes Genes in the panel.
#' @param group_by Categorical metadata column.
#' @return Data frame with \code{gene}, \code{group}, \code{mean_expr},
#'   \code{pct_expr} (percent of cells with expression > 0) and \code{scaled}
#'   (per-gene z-score of the group means).
#' @export
sparks_dot_data <- function(panel, genes, group_by) {
  stopifnot(inherits(panel, "sparks_panel"))
  genes <- intersect(genes, names(panel$genes))
  g <- .panel_meta_values(panel, group_by)
  if (!length(genes)) {
    return(data.frame(gene = character(), group = character(), mean_expr = numeric(),
                      pct_expr = numeric(), scaled = numeric()))
  }
  rows <- lapply(genes, function(gene) {
    v <- .panel_gene_values(panel, gene)
    mean_expr <- tapply(v, g, mean)
    pct <- tapply(v > 0, g, function(z) 100 * mean(z))
    sc <- if (stats::sd(mean_expr, na.rm = TRUE) > 0)
      as.numeric(scale(as.numeric(mean_expr))) else rep(0, length(mean_expr))
    data.frame(gene = gene, group = names(mean_expr), mean_expr = as.numeric(mean_expr),
               pct_expr = as.numeric(pct), scaled = sc, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Dot plot of genes x groups
#' @inheritParams sparks_dot_data
#' @param mode Colour mode.
#' @return plotly object or \code{NULL}.
#' @export
sparks_fig_dotplot <- function(panel, genes, group_by, mode = "light") {
  .require_plotly()
  d <- sparks_dot_data(panel, genes, group_by)
  if (!nrow(d)) return(NULL)
  n_seq <- length(.SPARKS_SEQ)
  lim <- max(abs(d$scaled), 1)
  p <- plotly::plot_ly(
    d, type = "scatter", mode = "markers", x = ~gene, y = ~group,
    marker = list(
      size = 4 + 0.22 * d$pct_expr, sizemode = "diameter", color = d$scaled,
      cmin = -lim, cmax = lim,
      colorscale = lapply(seq_len(n_seq), function(i) list((i - 1) / (n_seq - 1), .SPARKS_SEQ[[i]])),
      colorbar = list(title = list(text = "Scaled mean"), thickness = 10, len = 0.6),
      line = list(color = .SPARKS_CHROME[[mode]]$surface, width = 1)),
    text = ~paste0(gene, " / ", group, "<br>mean ", signif(mean_expr, 3),
                   "<br>", round(pct_expr), "% of cells"),
    hoverinfo = "text")
  p <- plotly::layout(p, xaxis = list(type = "category", categoryarray = unique(d$gene)),
                      yaxis = list(type = "category", autorange = "reversed"))
  .fig_config(.fig_base(p, NULL, NULL, group_by, FALSE, mode,
                        margin = list(l = 100, r = 16, t = 16, b = 80)), "dotplot")
}

# ---- Serialisation ----------------------------------------------------------------------

#' Serialise a plotly figure to a JSON spec (\code{data}, \code{layout}, \code{config})
#' @param p plotly object.
#' @return JSON string, or \code{NULL} for a \code{NULL} figure.
#' @keywords internal
.figure_spec_json <- function(p) {
  if (is.null(p)) return(NULL)
  built <- plotly::plotly_build(p)$x
  spec <- list(data = built$data, layout = built$layout, config = built$config)
  as.character(jsonlite::toJSON(spec, auto_unbox = TRUE, null = "null", na = "null",
                                digits = 6, force = TRUE))
}
