#' Title: plot_tree_network
#'
#' Plot network
#'
#' @description
#'
#'
#' @param x Description of the first parameter.
#' @param y Description of the second parameter (if applicable).
#' @param ... Other optional parameters passed to methods.
#'
#' @details
#'
#' @return
#' Description of the object that the function returns.
#' If the function doesn't return anything meaningful, you can say `NULL`.
#'
#' @examples
#'
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text theme_graph graph_from_adjacency_matrix
#' @importFrom NetIndices TrophInd
#' @export
#'

plot_tree_network <- function(MW){

  # Define the colors
  cols_troph <- c(
    "#F3E79A",  # Level 1
    "#F9C483",  # Level 2
    "#ED7C97",  # Level 3
    "#A653A8",  # Level 4
    "#704D9E"   # Level 5
  )
  names(cols_troph) <- c("1", "2", "3", "4", "5")

  # Create layout matrix
  Troph <- TrophInd(Flow = MW,
                    Tij = t(MW))

  graph <- graph_from_adjacency_matrix(data.matrix(MW),weighted=TRUE)

  layout.matrix<-matrix(nrow=length(V(graph)),ncol=2)  # Rows equal to the number of vertices
  layout.matrix[,1]<-runif(length(V(graph))) # randomly assign along x-axis
  layout.matrix[,2] <- Troph$TL # y-axis value based on trophic level

  Degree_in <- colSums(MW)
  Degree <- colSums(MW) + rowSums(MW)
  Trophic_Level <- as.character(round(Troph$TL,0))

  # Plot
  ggraph(graph, layout = layout.matrix)+
    geom_edge_link(aes(edge_alpha = 0.1), edge_colour = "grey66", arrow=arrow(ends="last", angle=20, length=unit(0.15, "inches"), type="closed"), show.legend=F)+
    geom_node_point(aes(fill = Trophic_Level, size = Degree_in), shape = 21) +
    geom_node_text(aes(label = name), family = "serif", repel="true")+
    theme_graph()+
    scale_fill_manual(values = cols_troph)


}
