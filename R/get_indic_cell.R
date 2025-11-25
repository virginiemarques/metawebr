#' Title: metaweb_mod_parameters
#'
#' Find the parameter for the metaweb
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
#' @importFrom parallel mclapply
#' @importFrom NetIndices TrophInd
#' @importFrom igraph graph.adjacency cluster_walktrap modularity diameter degree closeness articulation.points transitivity graph.coreness average.path.length shortest.paths distances
#'
#'
#' @export


get_indic_cells <- function(P_A_data, Lniche, mc.cores=1) {

  reseau_cell <-  parallel::mclapply(1:nrow(P_A_data),mc.cores=mc.cores,function(i){
    if((i%%100)==0) cat("i=",i,"\n")
    Names <- colnames(P_A_data)[which(P_A_data[i, ] > 0)]
    Calc_indic_proba(x=Lniche[Names,Names])
  })
  reseau_cell
} # get_indic_cells


#' Title: Calc_indic_proba
#'
#'
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
#'
#' @import NetIndices
#' @import igraph
#' @import gtools
#'
#'
# Helper function
Calc_indic_proba  <- function(x=reseau_cell[[1]], bin_threshold=1) {

  #cat("### Indicator proba calculation ###", "\n")
  Species <- nrow(x)
  Connectance_p <- sum(x)/Species^2

  #cat("### Binary calculation indicators ###","\n")
  bin_net <- round(x,3); bin_net[bin_net<bin_threshold] <- 0 ############## transformation from weighted web to binary web
  bin_net[bin_net>=bin_threshold] <- 1
  binary_indic <- get_binary_indic(web=bin_net, S=Species)

  web=data.matrix(bin_net)
  #cat("### Path calculation indicators ###","\n")
  Path_stat <- get_path_stats(web=web)

  #cat("### Igraph indicator calculation ###", "\n")
  igraph_res <- Calc_indic_igraph(web=web)

  R=0

  c(Species=Species,Connectance_p=Connectance_p,binary_indic,igraph_res,Path_stat, redundancy=R)

}  # end of function Calc_ind

#' Title: get_binary_indic
#'
#' @import NetIndices
#' @import igraph
#' @import gtools
#'
#'
# Helper function

get_binary_indic <- function(web=mat, S=Species){

  TL <- NetIndices::TrophInd(Flow =web,Tij = t(web))
  length_chain <- ceiling(round(max(TL[,1]),1))
  TL_moy  <- mean(TL[,1])
  Omn_moy <- mean(TL[,2])
  if ( nrow(TL)>50){
    Plankt <- length(which(TL[,1] > 2.4 & TL[,1] < 3.7 )) / nrow(TL) #proportion of plantivores
    HTI <- length(which(TL[,1] > 4 )) / nrow(TL) #High Trophic level Indicator (Bourdeau et al. 2016)
    MTI <- mean(TL[which(TL[,1] > 3.25 ),1])
  }else{ Plankt <- HTI <- MTI <- NA}

  Link <- sum(web); Link_max <- nrow(web)^2
  Connectance <- Link/Link_max
  b <- log2(Link) / (log2(S)-1) #parameter of the power law between L and S #Carpentier et al. 2021

  web <- web[-grep("Producteur",rownames(web), fixed=T),-grep("Producteur", colnames(web), fixed=T)]
  if( length(web)>1){
    Nprey <- colSums(web)  ### Nombre de proie qu'a chaque espèce
    Npred <- rowSums(web)  ### Nombre de prédateur qu'a chaque espece
    S <- nrow(web)
  }else{ Npred <- Nprey <- NA ; S <- 1}
  Ntop  <- ((sum(Npred == 0 & Nprey>=1)) / S)  #proportion of top predator
  Nbas  <- ((sum(Npred >= 0 & Nprey<=1)) / S) #proportion of herbivorous+planktivorous
  Nint  <- 1 - (Ntop + Nbas)

  Nprey <- Nprey[Nprey!=0] #espèces avec des proies
  Npred <- Npred[Npred!=0] #espèces avec des prédateurs

  Vul <-  mean(Npred) ; Vulsd<- sd(Npred)
  Gen <- mean(Nprey); Gensd <- sd(Nprey)

  c(Link=Link,Link_max=Link_max,Connectance=Connectance, b_power_law=b, Ntop=Ntop,Nbas=Nbas,Nint=Nint,Vul=Vul,Vulsd=Vulsd,Gen=Gen,Gensd=Gensd,TL_moy=TL_moy,
    length_chain=length_chain,Omn_moy=Omn_moy, Planktivores=Plankt, HTI=HTI, MTI=MTI)

} # end of get_binary_indic

#' Title: Calc_indic_igraph
#'
#' @import NetIndices
#' @import igraph
#' @import gtools
#'
#'
# Helper function


# Calc_indic_igraph<- function(web=mat){
#
#   web <-  igraph::graph.adjacency(web,weighted=NULL)
#
#   Mod <- modularity(walktrap.community(web))
#   Diam <- diameter(web)
#
#   Din <- igraph::degree(web, mode=c("in"))
#   Din_stat <- c(Din_mean=mean(Din),Din_min = min(Din),Din_max = max(Din))
#
#   Dout <- igraph::degree(web, mode=c("out"))
#   Dout_stat <- c(Dout_mean=mean(Dout),Dout_min = min(Dout),Dout_max = max(Dout))
#
#   #Bet <- betweenness(web)
#   #Bet_stat <- c(Betweenness_mean=mean(Bet),Betweenness_min = min(Bet),Betweenness_max = max(Bet))
#
#   Clos <- closeness(web,normalized = F)
#   Clos_stat <- c(Closeness_mean=mean(Clos, na.rm=T),Closeness_min = min(Clos, na.rm=T),Closeness_max = max(Clos, na.rm=T))
#
#   Nb_artpt <- length(articulation.points(web))
#   Trans <-  transitivity(web) # clustering coef
#
#   Cor <- graph.coreness(web)
#   Cor_stat <- c(Coreness_mean=mean(Cor, na.rm=T),Coreness_min=min(Cor, na.rm=T),Coreness_max=max(Cor, na.rm=T))
#
#   #Centrality <-  evcent(web)$vector
#   #c(Modularity=Mod, Diameter= Diam, Din_stat,Dout_stat,Bet_stat,Clos_stat,Nb_articulate_point = Nb_artpt,Transitivity=Trans,Cor_stat)
#   c(Modularity=Mod,Diameter= Diam, Din_stat,Dout_stat,Clos_stat,Nb_articulate_point = Nb_artpt,Transitivity=Trans,Cor_stat)
#
# } # end of function Calc_indic_igraph


# Chatgpt's suggestions
Calc_indic_igraph <- function(web = mat) {
  web <- igraph::graph.adjacency(web, weighted = NULL)

  # Community detection + modularity
  wc <- igraph::cluster_walktrap(web)
  Mod <- igraph::modularity(wc)

  # Diameter
  Diam <- igraph::diameter(web)

  # In-degree stats
  Din <- igraph::degree(web, mode = "in")
  Din_stat <- c(Din_mean = mean(Din), Din_min = min(Din), Din_max = max(Din))

  # Out-degree stats
  Dout <- igraph::degree(web, mode = "out")
  Dout_stat <- c(Dout_mean = mean(Dout), Dout_min = min(Dout), Dout_max = max(Dout))

  # Closeness
  Clos <- igraph::closeness(web, normalized = FALSE)
  Clos_stat <- c(
    Closeness_mean = mean(Clos, na.rm = TRUE),
    Closeness_min = min(Clos, na.rm = TRUE),
    Closeness_max = max(Clos, na.rm = TRUE)
  )

  # Articulation points
  Nb_artpt <- length(igraph::articulation.points(web))

  # Transitivity
  Trans <- igraph::transitivity(web)

  # Coreness
  Cor <- igraph::graph.coreness(web)
  Cor_stat <- c(
    Coreness_mean = mean(Cor, na.rm = TRUE),
    Coreness_min = min(Cor, na.rm = TRUE),
    Coreness_max = max(Cor, na.rm = TRUE)
  )

  # Combine all results
  c(
    Modularity = Mod,
    Diameter = Diam,
    Din_stat,
    Dout_stat,
    Clos_stat,
    Nb_articulate_point = Nb_artpt,
    Transitivity = Trans,
    Cor_stat
  )
}

#' Title: get_path_stats
#'
#' @import NetIndices
#' @import igraph
#' @import gtools
#'
#'

get_path_stats <- function(web=Path_net){
  web <-  igraph::graph.adjacency(web,mode="directed", weighted=NULL)
  Averpth_length <- igraph::average.path.length(web) #mean shortest path to each vertices

  b <- table(igraph::shortest.paths(web,mode="out"))
  b <- b[-which(names(b)==0 | names(b)==Inf)]
  c <- rep(NA,5) ; names(c) <- seq(1,5,1)
  d <- b/ sum(b)
  for (i in 1:length(d)){ c[i] <- d[i]}
  names(c) <- paste("Shortest_path_",names(c),sep="")

  c(Mean_path_length=Averpth_length,c)

} # end of get_path_stats


## Chatgpt's suggestion to improve below
#. get_path_stats <- function(web = Path_net) {
#.   # Convert adjacency matrix to igraph object
#.   web_graph <- igraph::graph.adjacency(web, mode = "directed", weighted = NULL)
#.
#.   # Compute average path length
#.   mean_path_length <- igraph::average.path.length(web_graph)
#.
#.   # Compute all shortest path distances (mode = "out")
#.   dist_mat <- igraph::distances(web_graph, mode = "out")
#.
#.   # Flatten the distance matrix and remove 0 and Inf distances
#.   dist_vec <- as.vector(dist_mat)
#.   dist_vec <- dist_vec[!dist_vec %in% c(0, Inf)]
#.
#.   # Compute frequency table
#.   freq_table <- table(dist_vec)
#.
#.   # Normalize frequencies (relative proportions)
#.   rel_freq <- freq_table / sum(freq_table)
#.
#.   # Prepare output vector for first 5 shortest path lengths
#.   out_vec <- rep(NA, 5)
#.   names(out_vec) <- paste0("Shortest_path_", 1:5)
#.
#.   # Fill in the observed proportions
#.   matched_lengths <- intersect(names(rel_freq), names(out_vec))
#.   out_vec[matched_lengths] <- rel_freq[matched_lengths]
#.
#.   # Combine with mean path length
#.   c(Mean_path_length = mean_path_length, out_vec)
#. }




