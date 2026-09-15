#' Predict the probability of every trophic interaction
#'
#' Applies the calibrated allometric niche model (Gravel et al. 2013) to every
#' pair of taxa in `data_traits`.
#'
#' @param data_traits Data frame with one row per taxon and taxon names as row
#'   names, as returned by [get_traits()]. Must contain `length_column`.
#' @param path_pars Model parameters (`a0`, `a1`, `b0`, `b1`): the path to a
#'   file written by [metaweb_mod_parameters()], the list it returns, or a
#'   numeric vector of length 4.
#' @param length_column Name of the column holding body length (cm).
#'
#' @return A square matrix of interaction probabilities with **prey in rows and
#'   predators in columns**: `mw[i, j]` is the probability that taxon `j` eats
#'   taxon `i`.
#'
#' @references Gravel, D., Poisot, T., Albouy, C., Velez, L. & Mouillot, D.
#'   (2013). Inferring food web structure from predator-prey body size
#'   relationships. *Methods in Ecology and Evolution*, 4, 1083-1090.
#'
#' @examples
#' traits <- data.frame(
#'   CommonLengthEstim = c(100, 30, 5),
#'   row.names = c("Gadus_morhua", "Clupea_harengus", "Gasterosteus_aculeatus")
#' )
#' pars <- c(a0 = -0.56, a1 = 0.96, b0 = 0.12, b1 = 0.14)
#' round(apply_model_metaweb(traits, pars), 2)
#'
#' @export
apply_model_metaweb <- function(data_traits,
                                path_pars,
                                length_column = "CommonLengthEstim") {
  check_traits(data_traits, length_column)
  pars <- as_pars(path_pars)

  len <- data_traits[[length_column]]
  if (!is.numeric(len)) {
    stop(sprintf("Column `%s` must be numeric.", length_column), call. = FALSE)
  }
  bad <- is.na(len) | len <= 0
  if (any(bad)) {
    stop(sprintf("%d taxa have a missing or non-positive `%s`: %s. Fill or remove them first (see clean_traits()).",
                 sum(bad), length_column, format_taxa(rownames(data_traits)[bad])),
         call. = FALSE)
  }

  size <- setNames(log10(len), rownames(data_traits))
  outer(size, size, pLMFitted, Pars = pars)
}

# Conditional probability that a predator of log10 length MPred eats a prey of
# log10 length MPrey (vectorised).
pLMFitted <- function(MPrey, MPred, Pars) {
  Pars <- as_pars(Pars)
  o <- Pars[["a0"]] + Pars[["a1"]] * MPred   # optimal prey size
  r <- Pars[["b0"]] + Pars[["b1"]] * MPred   # niche range
  exp(-(o - MPrey)^2 / 2 / r^2)
}
