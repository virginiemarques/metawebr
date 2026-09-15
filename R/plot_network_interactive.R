#' Plot an interactive food web
#'
#' Draws a metaweb, or a local food web, as an interactive network with the
#' layout of [plot_tree_network()]: the height of each node is its trophic level
#' (computed on the full metaweb), nodes are coloured by rounded trophic level
#' and sized by number of prey. Horizontal positions are random: call
#' [set.seed()] first for a reproducible layout.
#'
#' Click a taxon to highlight its links and fade the rest of the web. Hold Ctrl
#' (Cmd on macOS) while clicking to select several taxa, or pick a taxon by name
#' in the menu above the plot. Click an empty area to clear the selection. Hover
#' a taxon or a link for details; scroll to zoom and drag to pan.
#'
#' Every non-zero value of `MW` is drawn as a link, so binarise a probability
#' metaweb first (see [correct_metaweb_proba()]). Requires the visNetwork
#' package.
#'
#' @inheritParams plot_tree_network
#' @param links Links highlighted for the selected taxa: `"both"` (prey and
#'   predators), `"prey"` or `"predators"`.
#' @param height,width Size of the plot, in CSS units (e.g. `"700px"` or
#'   `"100%"`).
#'
#' @return A visNetwork htmlwidget. Print it to show it in the viewer or in an
#'   R Markdown document, or save it as a stand-alone HTML file with
#'   `htmlwidgets::saveWidget()`.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' web <- plot_network_interactive(mw_bin, sp_to_keep = marseille)
#' web
#' htmlwidgets::saveWidget(web, "foodweb.html")
#' }
#'
#' @export
plot_network_interactive <- function(MW, sp_to_keep = NULL, links = c("both", "prey", "predators"),
                                     height = "700px", width = "100%") {
  rlang::check_installed("visNetwork", reason = "to draw interactive food webs.")
  links <- match.arg(links)
  MW <- check_metaweb(MW)
  MW_sub <- subset_metaweb(MW, sp_to_keep)

  nodes <- trophic_graph_nodes(MW_sub, trophic_levels(MW))
  labels <- gsub("_", " ", nodes$name)
  spread <- max(1200, 25 * nrow(nodes))  # canvas width, in pixels

  vis_nodes <- data.frame(
    id = nodes$name,
    label = labels,
    group = nodes$Trophic_Level,
    value = nodes$Degree_in,
    x = nodes$x * spread,
    y = -nodes$TL * 250,  # the vis.js y axis points down
    title = sprintf("<b>%s</b><br>Trophic level: %.2f<br>Prey: %d<br>Predators: %d",
                    labels, nodes$TL, colSums(MW_sub > 0), rowSums(MW_sub > 0)),
    stringsAsFactors = FALSE
  )

  # Links go from prey (rows) to predators (columns)
  idx <- which(MW_sub > 0, arr.ind = TRUE)
  vis_edges <- data.frame(
    from = rownames(MW_sub)[idx[, 1]],
    to = colnames(MW_sub)[idx[, 2]],
    title = sprintf("%s &rarr; %s", labels[idx[, 1]], labels[idx[, 2]]),
    stringsAsFactors = FALSE
  )

  # "from" steps up to the prey of the selected taxa, "to" to their predators
  degree <- switch(links,
                   both = list(from = 1, to = 1),
                   prey = list(from = 1, to = 0),
                   predators = list(from = 0, to = 1))

  web <- visNetwork::visNetwork(vis_nodes, vis_edges, height = height, width = width)
  for (level in names(cols_troph)) {
    col <- cols_troph[[level]]
    web <- visNetwork::visGroups(web, groupname = level, shape = "dot",
                                 color = list(background = col, border = "#3a3a3a",
                                              highlight = list(background = col, border = "#000000"),
                                              hover = list(background = col, border = "#000000")))
  }

  levels_present <- sort(unique(nodes$Trophic_Level))
  legend_nodes <- data.frame(label = levels_present, shape = "dot",
                             color = unname(cols_troph[levels_present]), stringsAsFactors = FALSE)

  web <- visNetwork::visNodes(web, borderWidth = 1, borderWidthSelected = 3,
                              scaling = list(min = 8, max = 30), font = list(size = 16))
  web <- visNetwork::visEdges(web, smooth = FALSE, selectionWidth = 1.5,
                              arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5)),
                              color = list(color = "rgba(120,120,120,0.35)",
                                           highlight = "#1c5cab", hover = "#1c5cab"))
  web <- visNetwork::visOptions(web,
                                highlightNearest = list(enabled = TRUE, algorithm = "hierarchical",
                                                        degree = degree,
                                                        hideColor = "rgba(200,200,200,0.25)"),
                                nodesIdSelection = list(enabled = TRUE, useLabels = TRUE,
                                                        main = "Select a taxon"))
  web <- visNetwork::visInteraction(web, multiselect = TRUE, hover = TRUE, tooltipDelay = 100)
  web <- visNetwork::visPhysics(web, enabled = FALSE)
  visNetwork::visLegend(web, useGroups = FALSE, addNodes = legend_nodes, position = "right",
                        main = "Trophic level", width = 0.12)
}
