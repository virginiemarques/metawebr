#' Compute food-web indices for each site
#'
#' Extracts the local food web of each site (the sub-network of the metaweb
#' formed by the taxa present) and computes network indices for it.
#' Optionally compares the indices with null models (see [get_indic_null()]).
#'
#' @param P_A_data Site x taxon matrix with site names as row names and taxon
#'   names as column names, as returned by [make_species_station_matrix()].
#'   Values > 0 are presences.
#' @param Lniche Binary metaweb (prey in rows, predators in columns), e.g. from
#'   [correct_metaweb_proba()].
#' @param mc.cores Number of cores. Values > 1 use forked processes, which are
#'   not available on Windows.
#' @param null_model `"none"` (default) for observed indices only, or a null
#'   model passed to [get_indic_null()]: `"equiprobable"` or `"frequency"`.
#' @param n_null Number of null communities per site, when `null_model` is not
#'   `"none"`.
#'
#' @details Taxa of `P_A_data` absent from `Lniche` are dropped with a warning.
#'
#' Indices:
#' * `Species`: number of nodes, producers included. `Connectance_p`: sum of
#'   link values / `Species`^2.
#' * `Link`, `Link_max`, `Connectance`, `b_power_law`: number of links, nodes^2,
#'   their ratio, and \eqn{\log_2 L / (\log_2 S - 1)}.
#' * `Ntop`, `Nbas`, `Nint`: fractions of top (no predator, at least one prey),
#'   basal (at most one prey) and intermediate taxa, producers excluded.
#' * `Vul`, `Vulsd`, `Gen`, `Gensd`: mean and sd of the number of predators
#'   (vulnerability) and prey (generality), producers excluded.
#' * `TL_moy`, `length_chain`, `Omn_moy`: mean trophic level, maximum trophic
#'   level rounded up, mean omnivory index.
#' * `Planktivores`, `HTI`, `MTI`: fraction of taxa with a trophic level in
#'   (2.4, 3.7), fraction above 4, mean trophic level of taxa above 3.25.
#'   Only for webs with more than 50 nodes, otherwise `NA`.
#' * `Modularity` (walktrap communities), `Diameter`, `Din_*` / `Dout_*`
#'   (in- and out-degree mean, min, max), `Closeness_*`,
#'   `Nb_articulate_point`, `Transitivity`, `Coreness_*`.
#' * `Mean_path_length` and `Shortest_path_1` to `Shortest_path_5`: mean
#'   shortest path and fraction of shortest paths of each length (`NA` when
#'   there is none).
#' * `redundancy`: not implemented yet, always 0.
#'
#' @return If `null_model = "none"`, a data frame with one row per site: a
#'   `site` column followed by the indices. Otherwise a list with `indices`
#'   (that data frame) and `null_model` (the output of [get_indic_null()]).
#'
#' @export
get_indic_cells <- function(P_A_data, Lniche, mc.cores = 1,
                            null_model = c("none", "equiprobable", "frequency"),
                            n_null = 100) {
  null_model <- match.arg(null_model)
  Lniche <- check_metaweb(Lniche, "Lniche")
  pa <- prepare_sites(P_A_data, Lniche)
  indices <- compute_site_indices(pa, Lniche, mc.cores)
  if (null_model == "none") return(indices)

  list(indices = indices,
       null_model = null_model_indices(pa, Lniche, method = null_model, n_null = n_null,
                                       indices = NULL, mc.cores = mc.cores, observed = indices))
}

# Logical site x taxon matrix restricted to the taxa of the metaweb
prepare_sites <- function(P_A_data, Lniche) {
  if (is.data.frame(P_A_data)) P_A_data <- data.matrix(P_A_data)
  if (!is.matrix(P_A_data) || !is.numeric(P_A_data)) {
    stop("`P_A_data` must be a numeric site x taxon matrix.", call. = FALSE)
  }
  if (is.null(colnames(P_A_data))) {
    stop("`P_A_data` must have taxon names as column names.", call. = FALSE)
  }
  if (anyNA(P_A_data)) {
    stop("`P_A_data` contains missing values.", call. = FALSE)
  }
  if (is.null(rownames(P_A_data))) {
    rownames(P_A_data) <- paste0("site_", seq_len(nrow(P_A_data)))
  }

  pa <- P_A_data > 0
  absent <- !colnames(pa) %in% rownames(Lniche)
  dropped <- colnames(pa)[absent & colSums(pa) > 0]
  if (length(dropped) > 0) {
    n_sites <- sum(rowSums(pa[, absent, drop = FALSE]) > 0)
    warning(sprintf("%d taxa are not in the metaweb and were dropped (they occur in %d of %d sites): %s.",
                    length(dropped), n_sites, nrow(pa), format_taxa(dropped)), call. = FALSE)
  }
  pa[, !absent, drop = FALSE]
}

# Indices of the food web formed by `taxa`, or NULL for an empty web
web_indices <- function(taxa, Lniche) {
  if (length(taxa) == 0) return(NULL)
  Calc_indic_proba(Lniche[taxa, taxa, drop = FALSE])
}

# Row-bind index vectors, filling empty webs (NULL) with NA
bind_indices <- function(res) {
  template <- Find(Negate(is.null), res)
  if (is.null(template)) stop("All food webs are empty.", call. = FALSE)
  empty <- setNames(rep(NA_real_, length(template)), names(template))
  do.call(rbind, lapply(res, function(x) if (is.null(x)) empty else x))
}

compute_site_indices <- function(pa, Lniche, mc.cores = 1) {
  res <- run_parallel(seq_len(nrow(pa)), function(i) {
    web_indices(colnames(pa)[pa[i, ]], Lniche)
  }, mc.cores)
  data.frame(site = rownames(pa), bind_indices(res), row.names = NULL, check.names = FALSE)
}

# Indices of one food web
Calc_indic_proba <- function(x, bin_threshold = 1) {
  Species <- nrow(x)
  Connectance_p <- sum(x) / Species^2

  bin_net <- (round(x, 3) >= bin_threshold) * 1
  binary_indic <- get_binary_indic(web = bin_net, S = Species)

  g <- igraph::graph_from_adjacency_matrix(bin_net, mode = "directed")
  igraph_res <- Calc_indic_igraph(g)
  Path_stat <- get_path_stats(g)

  c(Species = Species, Connectance_p = Connectance_p, binary_indic, igraph_res, Path_stat,
    redundancy = 0)
}

# Trophic-level and link-based indices of a binary web (prey in rows)
get_binary_indic <- function(web, S = nrow(web)) {
  TL <- trophic_levels(web)
  length_chain <- ceiling(round(max(TL$TL), 1))
  TL_moy <- mean(TL$TL)
  Omn_moy <- mean(TL$OI)
  if (nrow(TL) > 50) {
    # Binary webs often give trophic levels exactly on a cut-off (e.g. 3.25);
    # rounding avoids classifying them by floating-point noise
    tl <- round(TL$TL, 10)
    Plankt <- sum(tl > 2.4 & tl < 3.7) / nrow(TL)  # proportion of planktivores
    HTI <- sum(tl > 4) / nrow(TL)  # High Trophic level Indicator
    MTI <- mean(TL$TL[tl > 3.25])
  } else {
    Plankt <- HTI <- MTI <- NA
  }

  Link <- sum(web)
  Link_max <- nrow(web)^2
  Connectance <- Link / Link_max
  b <- log2(Link) / (log2(S) - 1)  # exponent of the link-species scaling

  # Top / basal / intermediate taxa and degrees, without the producer nodes
  is_taxon <- !rownames(web) %in% producer_names
  web <- web[is_taxon, is_taxon, drop = FALSE]
  if (nrow(web) > 1) {
    Nprey <- colSums(web)  # number of prey of each taxon
    Npred <- rowSums(web)  # number of predators of each taxon
    S <- nrow(web)
  } else {
    Npred <- Nprey <- NA
    S <- 1
  }
  Ntop <- sum(Npred == 0 & Nprey >= 1) / S
  Nbas <- sum(Npred >= 0 & Nprey <= 1) / S
  Nint <- 1 - (Ntop + Nbas)

  Nprey <- Nprey[Nprey != 0]
  Npred <- Npred[Npred != 0]

  out <- c(Link = Link, Link_max = Link_max, Connectance = Connectance, b_power_law = b,
           Ntop = Ntop, Nbas = Nbas, Nint = Nint,
           Vul = mean(Npred), Vulsd = sd(Npred), Gen = mean(Nprey), Gensd = sd(Nprey),
           TL_moy = TL_moy, length_chain = length_chain, Omn_moy = Omn_moy,
           Planktivores = Plankt, HTI = HTI, MTI = MTI)
  out[is.nan(out)] <- NA
  out
}

# mean, min, max ignoring NA (NA when nothing is left)
summary_stats <- function(x, prefix) {
  x <- x[!is.na(x)]
  out <- if (length(x) > 0) c(mean(x), min(x), max(x)) else rep(NA_real_, 3)
  setNames(out, paste0(prefix, c("_mean", "_min", "_max")))
}

# igraph-based indices of a directed graph
Calc_indic_igraph <- function(g) {
  wc <- igraph::cluster_walktrap(g)
  c(Modularity = igraph::modularity(wc),
    Diameter = igraph::diameter(g),
    summary_stats(igraph::degree(g, mode = "in"), "Din"),
    summary_stats(igraph::degree(g, mode = "out"), "Dout"),
    summary_stats(igraph::closeness(g, normalized = FALSE), "Closeness"),
    Nb_articulate_point = length(igraph::articulation_points(g)),
    Transitivity = igraph::transitivity(g),
    summary_stats(igraph::coreness(g), "Coreness"))
}

# Mean shortest path and fraction of shortest paths of length 1 to 5
get_path_stats <- function(g) {
  d <- igraph::distances(g, mode = "out")
  d <- d[is.finite(d) & d > 0]
  prop <- tabulate(d, nbins = 5) / length(d)
  prop[!is.finite(prop) | prop == 0] <- NA
  c(Mean_path_length = if (length(d) > 0) mean(d) else NaN,
    setNames(prop, paste0("Shortest_path_", 1:5)))
}
