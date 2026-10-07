#' Binarise a metaweb
#'
#' Turns interaction probabilities into 0/1 links. By default the threshold is
#' chosen with [get_proba()].
#'
#' @param MW Metaweb of interaction probabilities (prey in rows, predators in
#'   columns), e.g. from [correct_metaweb_fish()].
#' @param df_traits Traits with taxon names as row names and a `TrophicLevel`
#'   column, for `method = "trophic_level"`. Not needed when `threshold` is
#'   given.
#' @param threshold Optional probability threshold. `NULL` (default) chooses it
#'   with [get_proba()].
#' @param plot_diag Draw the threshold diagnostic plot?
#' @param method,links How [get_proba()] chooses the threshold.
#'
#' @return A 0/1 matrix with the same dimensions as `MW`; links with a
#'   probability >= threshold are set to 1.
#'
#' @export
correct_metaweb_proba <- function(MW, df_traits = NULL, threshold = NULL, plot_diag = TRUE,
                                  method = c("trophic_level", "tss"), links = NULL) {
  MW <- check_metaweb(MW)
  if (is.null(threshold)) {
    threshold <- get_proba(MW, df_traits, plot_diag = plot_diag, method = method, links = links)
    message(sprintf("Threshold value chosen is %g", threshold))
  } else if (!is.numeric(threshold) || length(threshold) != 1 || is.na(threshold) ||
             threshold < 0 || threshold > 1) {
    stop("`threshold` must be a single number between 0 and 1.", call. = FALSE)
  }
  out <- (MW >= threshold) * 1
  attr(out, "corrections") <- NULL  # links recorded by correct_metaweb_fish() no longer apply
  out
}

#' Choose the probability threshold of a metaweb
#'
#' Tries thresholds from 0.1 to 0.995 and returns the best one, by one of two
#' criteria:
#'
#' * `method = "trophic_level"` (default): for each threshold, the metaweb is
#'   binarised, the trophic level of every taxon is computed from the binary
#'   web, and the absolute differences with the observed trophic levels
#'   (`TrophicLevel`, e.g. from FishBase) are summed. The threshold with the
#'   smallest error is returned. This only constrains the average position of
#'   taxa in the web, not which links exist.
#' * `method = "tss"`: the threshold maximising the true skill statistic
#'   (sensitivity + specificity - 1) against observed links (`links`), as in
#'   [validate_metaweb()]. Use links that were not used to build the metaweb.
#'
#' @inheritParams correct_metaweb_proba
#' @param df_traits Traits with taxon names as row names and a `TrophicLevel`
#'   column, for `method = "trophic_level"`. Taxa are matched to the metaweb by
#'   name.
#' @param method `"trophic_level"` or `"tss"`.
#' @param links For `method = "tss"`: observed links, a data frame with the
#'   columns `predator` and `prey` (see [validate_metaweb()]), e.g. from
#'   [get_diet_links()].
#'
#' @return The selected threshold (a single number).
#'
#' @export
get_proba <- function(MW, df_traits = NULL, plot_diag = TRUE, method = c("trophic_level", "tss"),
                      links = NULL) {
  method <- match.arg(method)
  MW <- check_metaweb(MW)
  thres <- threshold_grid

  if (method == "tss") {
    if (is.null(links)) stop("`method = \"tss\"` needs observed `links`.", call. = FALSE)
    pairs <- validation_pairs(MW, links)
    score <- tss_curve(pairs$predicted, pairs$observed)
    best <- thres[which.max(score)]
    ylab <- "TSS (observed links)"
  } else {
    if (is.null(df_traits)) stop("`method = \"trophic_level\"` needs `df_traits`.", call. = FALSE)
    check_traits(df_traits, "TrophicLevel", arg = "df_traits")
    common <- intersect(rownames(MW), rownames(df_traits)[!is.na(df_traits$TrophicLevel)])
    if (length(common) == 0) {
      stop("No taxon of `MW` has a trophic level in `df_traits` (row names must match).", call. = FALSE)
    }
    observed_tl <- df_traits[common, "TrophicLevel"]
    rounded <- round(MW, 4)
    score <- vapply(thres, function(th) {
      inferred <- trophic_levels((rounded >= th) * 1)
      sum(abs(inferred[common, "TL"] - observed_tl))
    }, numeric(1))
    best <- thres[which.min(score)]
    ylab <- "Sum of |observed TL - inferred TL|"
  }

  if (isTRUE(plot_diag)) {
    plot(thres, score, pch = 1, xlab = "Binary threshold", ylab = ylab)
    best_score <- score[thres == best]
    points(best, best_score, col = "red", pch = 20)
    abline(h = best_score, col = "red", lty = 2, lwd = 0.5)
  }
  best
}
