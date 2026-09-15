#' Binarise a metaweb
#'
#' Turns interaction probabilities into 0/1 links. By default the threshold is
#' chosen with [get_proba()].
#'
#' @param MW Metaweb of interaction probabilities (prey in rows, predators in
#'   columns), e.g. from [correct_metaweb_fish()].
#' @param df_traits Traits with taxon names as row names and a `TrophicLevel`
#'   column. Not needed when `threshold` is given.
#' @param threshold Optional probability threshold. `NULL` (default) chooses it
#'   with [get_proba()].
#' @param plot_diag Draw the threshold diagnostic plot?
#'
#' @return A 0/1 matrix with the same dimensions as `MW`; links with a
#'   probability >= threshold are set to 1.
#'
#' @export
correct_metaweb_proba <- function(MW, df_traits, threshold = NULL, plot_diag = TRUE) {
  MW <- check_metaweb(MW)
  if (is.null(threshold)) {
    threshold <- get_proba(MW, df_traits, plot_diag = plot_diag)
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
#' Tries thresholds from 0.1 to 0.995. For each, the metaweb is binarised, the
#' trophic level of every taxon is computed from the binary web, and the
#' absolute differences with the observed trophic levels (`TrophicLevel`, e.g.
#' from FishBase) are summed. The threshold with the smallest error is returned.
#'
#' @inheritParams correct_metaweb_proba
#' @param df_traits Traits with taxon names as row names and a `TrophicLevel`
#'   column. Taxa are matched to the metaweb by name.
#'
#' @return The selected threshold (a single number).
#'
#' @export
get_proba <- function(MW, df_traits, plot_diag = TRUE) {
  MW <- check_metaweb(MW)
  check_traits(df_traits, "TrophicLevel", arg = "df_traits")
  common <- intersect(rownames(MW), rownames(df_traits)[!is.na(df_traits$TrophicLevel)])
  if (length(common) == 0) {
    stop("No taxon of `MW` has a trophic level in `df_traits` (row names must match).", call. = FALSE)
  }
  observed_tl <- df_traits[common, "TrophicLevel"]

  thres <- c(seq(0.1, 0.975, 0.025), 0.99, 0.995)
  rounded <- round(MW, 4)
  TLm <- vapply(thres, function(th) {
    inferred <- trophic_levels((rounded >= th) * 1)
    sum(abs(inferred[common, "TL"] - observed_tl))
  }, numeric(1))

  best <- thres[which.min(TLm)]
  if (isTRUE(plot_diag)) {
    plot(thres, TLm, pch = 1, xlab = "Binary threshold",
         ylab = "Sum of |observed TL - inferred TL|")
    points(best, min(TLm), col = "red", pch = 20)
    abline(h = min(TLm), col = "red", lty = 2, lwd = 0.5)
  }
  best
}
