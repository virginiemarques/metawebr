#' Compare local food-web indices with null models
#'
#' For each site, builds `n_null` random communities with the same number of
#' taxa as the site, drawn from the regional pool (the taxa present in at least
#' one site of `P_A_data` and in the metaweb). Their food webs are extracted
#' from the metaweb and the same indices as in [get_indic_cells()] are
#' computed. The observed value of each index is then compared with its null
#' distribution.
#'
#' Two null models are available:
#' * `"equiprobable"`: every taxon of the pool has the same chance to be drawn.
#' * `"frequency"`: taxa are drawn with a probability proportional to the
#'   number of sites where they occur, so common taxa are drawn more often.
#'
#' Producer nodes present at a site are kept in all its null communities.
#' All random draws are made before any parallel computation, so call
#' [set.seed()] first for reproducible results, whatever `mc.cores`.
#'
#' @inheritParams get_indic_cells
#' @param method Null model: `"equiprobable"` or `"frequency"`.
#' @param n_null Number of null communities per site.
#' @param indices Indices to report. By default all indices except those fixed
#'   by site richness (`Species`, `Link_max`, `redundancy`).
#' @param observed Optional output of [get_indic_cells()] for the same
#'   `P_A_data`, to avoid recomputing the observed indices.
#'
#' @return A data frame with one row per site and index:
#'   * `site`, `index`;
#'   * `observed`: observed value;
#'   * `null_mean`, `null_sd`: mean and standard deviation across null
#'     communities;
#'   * `ses`: standardised effect size, (observed - null_mean) / null_sd.
#'     Negative values mean lower than expected by chance. `NA` when the null
#'     distribution has no variance;
#'   * `p_value`: two-sided rank-based p-value,
#'     2 * min(P(null <= observed), P(null >= observed)), with
#'     P = (count + 1) / (n_null + 1).
#'
#' @seealso [plot_indic_null()] to plot the results.
#'
#' @examples
#' web <- matrix(0, 5, 5, dimnames = rep(list(c("A", "B", "C", "D", "E")), 2))
#' web["A", c("B", "C")] <- 1
#' web["B", "D"] <- 1
#' web["C", c("D", "E")] <- 1
#' sites <- rbind(s1 = c(A = 1, B = 1, C = 0, D = 1, E = 0),
#'                s2 = c(A = 1, B = 0, C = 1, D = 1, E = 1),
#'                s3 = c(A = 1, B = 1, C = 1, D = 0, E = 1))
#' set.seed(1)
#' res <- get_indic_null(sites, web, n_null = 20, indices = c("Link", "Connectance"))
#' res
#'
#' @export
get_indic_null <- function(P_A_data, Lniche,
                           method = c("equiprobable", "frequency"),
                           n_null = 100, indices = NULL, mc.cores = 1, observed = NULL) {
  method <- match.arg(method)
  Lniche <- check_metaweb(Lniche, "Lniche")
  pa <- prepare_sites(P_A_data, Lniche)
  if (is.null(observed)) {
    observed <- compute_site_indices(pa, Lniche, mc.cores)
  }
  null_model_indices(pa, Lniche, method, n_null, indices, mc.cores, observed)
}

null_model_indices <- function(pa, Lniche, method, n_null, indices, mc.cores, observed) {
  check_count(n_null, "n_null", min = 2)
  if (!is.data.frame(observed) || !identical(as.character(observed$site), rownames(pa))) {
    stop("`observed` must be the output of get_indic_cells() for the same sites.", call. = FALSE)
  }
  obs <- as.matrix(observed[, setdiff(names(observed), "site"), drop = FALSE])
  if (is.null(indices)) {
    indices <- setdiff(colnames(obs), c("Species", "Link_max", "redundancy"))
  }
  unknown <- setdiff(indices, colnames(obs))
  if (length(unknown) > 0) {
    stop(sprintf("Unknown index(es): %s.", paste(unknown, collapse = ", ")), call. = FALSE)
  }

  # Regional pool and draw weights
  is_producer <- colnames(pa) %in% producer_names
  pool_pa <- pa[, !is_producer, drop = FALSE]
  occurrences <- colSums(pool_pa)
  pool <- colnames(pool_pa)[occurrences > 0]
  if (length(pool) < 2) {
    stop("The regional pool must contain at least two taxa.", call. = FALSE)
  }
  weights <- if (method == "frequency") occurrences[pool] else NULL

  # All random communities are drawn here, in the main process
  draws <- lapply(seq_len(nrow(pa)), function(i) {
    richness <- sum(pool_pa[i, ])
    producers <- colnames(pa)[is_producer & pa[i, ]]
    lapply(seq_len(n_null), function(k) {
      c(pool[sample.int(length(pool), richness, prob = weights)], producers)
    })
  })

  null_values <- run_parallel(seq_along(draws), function(i) {
    bind_indices(lapply(draws[[i]], web_indices, Lniche = Lniche))
  }, mc.cores)

  out <- lapply(seq_len(nrow(pa)), function(i) {
    summarise_null(rownames(pa)[i], obs[i, indices], null_values[[i]][, indices, drop = FALSE])
  })
  do.call(rbind, out)
}

# Compare observed indices (named vector) with null values (n_null x index matrix)
summarise_null <- function(site, observed, null) {
  n_valid <- colSums(!is.na(null))
  null_mean <- colMeans(null, na.rm = TRUE)
  null_mean[is.nan(null_mean)] <- NA
  null_sd <- apply(null, 2, sd, na.rm = TRUE)

  ses <- (observed - null_mean) / null_sd
  ses[!is.finite(ses)] <- NA

  below <- colSums(sweep(null, 2, observed, "<="), na.rm = TRUE)
  above <- colSums(sweep(null, 2, observed, ">="), na.rm = TRUE)
  p_value <- pmin(1, 2 * (pmin(below, above) + 1) / (n_valid + 1))
  p_value[is.na(observed) | n_valid == 0] <- NA

  data.frame(site = site, index = names(observed), observed = unname(observed),
             null_mean = unname(null_mean), null_sd = unname(null_sd),
             ses = unname(ses), p_value = unname(p_value), row.names = NULL)
}

#' Plot null-model effect sizes
#'
#' Plots the standardised effect size (SES) of each index for each site.
#' Dashed lines mark SES = +/-1.96; coloured points are significant at `alpha`.
#'
#' @param null_res Output of [get_indic_null()], or the `null_model` element
#'   returned by [get_indic_cells()].
#' @param indices Indices to plot, in facet order. Default: all of `null_res`.
#' @param alpha Significance level used to colour the points.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_vline geom_point facet_wrap
#'   scale_colour_manual labs theme_bw
#' @export
plot_indic_null <- function(null_res, indices = NULL, alpha = 0.05) {
  required <- c("site", "index", "ses", "p_value")
  if (!is.data.frame(null_res) || !all(required %in% names(null_res))) {
    stop("`null_res` must be the output of get_indic_null().", call. = FALSE)
  }
  if (!is.null(indices)) {
    null_res <- null_res[null_res$index %in% indices, , drop = FALSE]
    null_res$index <- factor(null_res$index, levels = indices)
  }
  if (nrow(null_res) == 0) stop("None of `indices` are in `null_res`.", call. = FALSE)

  null_res$site <- factor(null_res$site, levels = rev(unique(null_res$site)))
  null_res$significant <- factor(!is.na(null_res$p_value) & null_res$p_value < alpha,
                                 levels = c(FALSE, TRUE))

  ggplot(null_res, aes(x = ses, y = site, colour = significant)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_vline(xintercept = c(-1.96, 1.96), colour = "grey70", linetype = "dashed") +
    geom_point(size = 2.5, na.rm = TRUE) +
    facet_wrap(~index, scales = "free_x") +
    scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#C0392B"),
                        labels = c(`FALSE` = "not significant", `TRUE` = sprintf("p < %s", alpha)),
                        name = NULL, drop = FALSE) +
    labs(x = "Standardised effect size (SES)", y = NULL) +
    theme_bw()
}
