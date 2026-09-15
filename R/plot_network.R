#' Plot a food web by trophic level
#'
#' Draws a metaweb, or a local food web, with the height of each node given by
#' its trophic level. Trophic levels are computed on the full metaweb, so
#' heights are comparable between sub-webs. Nodes are coloured by rounded
#' trophic level and sized by number of prey. Horizontal positions are random:
#' call [set.seed()] first for a reproducible layout.
#'
#' @param MW Metaweb, usually binary (prey in rows, predators in columns).
#' @param sp_to_keep Optional taxa to keep, to plot a local web. The producer
#'   nodes are always kept.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_tree_network <- function(MW, sp_to_keep = NULL) {
  MW <- check_metaweb(MW)
  MW_sub <- subset_metaweb(MW, sp_to_keep)
  plot_trophic_graph(MW_sub, trophic_levels(MW))
}

#' Plot the prey of focal taxa
#'
#' Plots one or more focal taxa with their prey, up to `order_interaction`
#' steps down the food chain (`2` adds the prey of the prey), optionally within
#' a local web. The layout is the same as in [plot_tree_network()].
#'
#' @inheritParams plot_tree_network
#' @param sp_focus Names of the focal taxa.
#' @param order_interaction Number of trophic steps to include below the focal
#'   taxa.
#' @param sp_to_keep Optional taxa to keep before looking for prey (species and
#'   interactions not in this vector are removed).
#'
#' @return A ggplot object.
#'
#' @export
plot_network_focus_sp <- function(MW, sp_focus, order_interaction = 1, sp_to_keep = NULL) {
  MW <- check_metaweb(MW)
  check_count(order_interaction, "order_interaction")
  MW_sub <- subset_metaweb(MW, sp_to_keep)

  found <- intersect(sp_focus, rownames(MW_sub))
  where <- if (is.null(sp_to_keep)) "the metaweb" else "the metaweb filtered with `sp_to_keep`"
  if (length(found) == 0) {
    stop(sprintf("None of `sp_focus` are in %s: %s.", where, format_taxa(sp_focus)), call. = FALSE)
  }
  if (length(found) < length(sp_focus)) {
    warning(sprintf("Focal taxa not in %s, ignored: %s.", where,
                    format_taxa(setdiff(sp_focus, found))), call. = FALSE)
  }

  # Focal taxa and their prey ("in" neighbours) up to order_interaction steps
  graph <- igraph::graph_from_adjacency_matrix(MW_sub, weighted = TRUE)
  neighbours <- igraph::neighborhood(graph, order = order_interaction, nodes = found, mode = "in")
  keep <- rownames(MW_sub)[sort(unique(unlist(lapply(neighbours, as.integer))))]

  plot_trophic_graph(MW_sub[keep, keep, drop = FALSE], trophic_levels(MW))
}

# Metaweb restricted to sp_to_keep (plus producers), with checks
subset_metaweb <- function(MW, sp_to_keep) {
  if (is.null(sp_to_keep)) return(MW)
  absent <- setdiff(sp_to_keep, rownames(MW))
  if (length(absent) == length(unique(sp_to_keep))) {
    stop("None of `sp_to_keep` are in the metaweb.", call. = FALSE)
  }
  if (length(absent) > 0) {
    warning(sprintf("%d taxa of `sp_to_keep` are not in the metaweb: %s.",
                    length(absent), format_taxa(absent)), call. = FALSE)
  }
  filter_adjacency(MW, sp_to_keep)
}

# Colours of the rounded trophic levels
cols_troph <- c(`1` = "#F3E79A", `2` = "#F9C483", `3` = "#ED7C97", `4` = "#A653A8", `5` = "#704D9E")

# One row per node of `web`, with its trophic level from `TL` (output of
# trophic_levels() on the full metaweb), rounded level, number of prey and a
# random horizontal position. Shared by the static and interactive plots.
trophic_graph_nodes <- function(web, TL) {
  taxa <- rownames(web)
  tl <- TL[taxa, "TL"]
  data.frame(
    name = taxa,
    TL = tl,
    Trophic_Level = as.character(pmin(pmax(round(tl), 1), 5)),
    Degree_in = unname(colSums(web)),
    x = runif(length(taxa)),
    stringsAsFactors = FALSE
  )
}

# Draw `web` with node heights from `TL` (output of trophic_levels() on the
# full metaweb)
plot_trophic_graph <- function(web, TL) {
  nodes <- trophic_graph_nodes(web, TL)
  graph <- igraph::graph_from_adjacency_matrix(web, weighted = TRUE)
  igraph::V(graph)$Trophic_Level <- nodes$Trophic_Level
  igraph::V(graph)$Degree_in <- nodes$Degree_in

  layout_matrix <- cbind(nodes$x, nodes$TL)  # random x, trophic level as y

  ggraph::ggraph(graph, layout = layout_matrix) +
    ggraph::geom_edge_link(edge_colour = "grey50", edge_alpha = 0.35,
                           arrow = arrow(ends = "last", angle = 20, length = unit(0.12, "inches"),
                                         type = "closed"),
                           end_cap = ggraph::circle(2.5, "mm")) +
    ggraph::geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
    ggraph::geom_node_text(aes(label = name), family = "serif", repel = TRUE) +
    ggraph::theme_graph(base_family = "sans") +
    scale_fill_manual(values = cols_troph, name = "Trophic level") +
    scale_size_continuous(name = "Number of prey")
}

# Keep the rows and columns of `keep_vec` (and, by default, the producers)
filter_adjacency <- function(adj_matrix, keep_vec, keep_producers = TRUE) {
  if (keep_producers) {
    keep_vec <- c(keep_vec, producer_names)
  }
  keep <- intersect(keep_vec, rownames(adj_matrix))
  adj_matrix[keep, keep, drop = FALSE]
}
