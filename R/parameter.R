#' Calibrate the allometric niche model
#'
#' Estimates the four parameters of the niche model of Gravel et al. (2013)
#' from observed predator-prey body lengths, by maximum likelihood with
#' generalized simulated annealing ([GenSA::GenSA()]).
#'
#' Before fitting, duplicated records (same predator, prey and lengths) are
#' removed and each predator-prey pair is randomly down-sampled to at most
#' `observation_resample` records, so that heavily sampled pairs don't dominate
#' the fit. Call [set.seed()] first for reproducible results.
#'
#' The model: a predator of log10 length \eqn{M} eats prey of log10 length
#' around \eqn{o = a_0 + a_1 M}, within a Gaussian niche of width
#' \eqn{r = b_0 + b_1 M}.
#'
#' @param data_path `NULL` (default) to use [df_interaction_fish], the path to
#'   a CSV file, or a data frame. Must contain the columns `predator`, `prey`,
#'   `standardised_predator_length` and `si_prey_length` (lengths in cm).
#' @param path_par_output Optional path where the fitted parameters are written
#'   with [utils::write.table()].
#' @param observation_resample Maximum number of records kept per
#'   predator-prey pair.
#' @param max_time Maximum running time of the optimisation, in seconds.
#' @param verbose Print the optimisation progress?
#'
#' @return A list with
#'   * `calibration_data`: starting values, from a linear regression of prey
#'     size on predator size;
#'   * `calibrated_data`: fitted parameters, a named vector `a0`, `a1`, `b0`, `b1`;
#'   * `neg_log_likelihood`: value of the objective at the optimum;
#'   * `n_obs`: number of records used.
#'
#'   The list can be passed directly to [apply_model_metaweb()].
#'
#' @references Gravel, D., Poisot, T., Albouy, C., Velez, L. & Mouillot, D.
#'   (2013). Inferring food web structure from predator-prey body size
#'   relationships. *Methods in Ecology and Evolution*, 4, 1083-1090.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' params <- metaweb_mod_parameters(path_par_output = "params_metaweb.txt")
#' params$calibrated_data
#' }
#'
#' @export
metaweb_mod_parameters <- function(data_path = NULL,
                                   path_par_output = NULL,
                                   observation_resample = 50,
                                   max_time = 1000,
                                   verbose = TRUE) {
  data_end <- read_interaction_data(data_path)
  check_count(observation_resample, "observation_resample")
  data <- prepare_calibration_data(data_end, observation_resample)

  # Starting values from a linear regression
  fit <- lm(MPrey ~ MPred, data = data)
  pars <- c(a0 = unname(fit$coefficients[1]), a1 = unname(fit$coefficients[2]),
            b0 = sd(fit$residuals), b1 = 0)
  par_lo <- c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
  par_hi <- c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)

  # Terms of the likelihood that don't depend on the parameters, computed once
  mean_prey <- mean(unique(data$MPrey))
  sd_prey <- sd(unique(data$MPrey))
  pM <- dnorm(data$MPrey, mean = mean_prey, sd = sd_prey)

  estim <- GenSA::GenSA(par = pars, fn = model, lower = par_lo, upper = par_hi,
                        control = list(verbose = verbose, max.time = max_time, smooth = FALSE),
                        data = data, mean_prey = mean_prey, sd_prey = sd_prey, pM = pM)

  fitted <- setNames(estim$par, c("a0", "a1", "b0", "b1"))
  if (!is.null(path_par_output)) {
    write.table(fitted, file = path_par_output)
  }

  list(calibration_data = data.frame(pars),
       calibrated_data = fitted,
       neg_log_likelihood = estim$value,
       n_obs = nrow(data))
}

# Deduplicate and down-sample the interaction records; returns log10 lengths.
prepare_calibration_data <- function(data_end, observation_resample = 50) {
  cols <- c("predator", "prey", "standardised_predator_length", "si_prey_length")
  missing_cols <- setdiff(cols, names(data_end))
  if (length(missing_cols) > 0) {
    stop(sprintf("Interaction data is missing column(s): %s.", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  pred_len <- data_end$standardised_predator_length
  prey_len <- data_end$si_prey_length
  ok <- is.finite(pred_len) & pred_len > 0 & is.finite(prey_len) & prey_len > 0
  if (any(!ok)) {
    warning(sprintf("%d records with a missing or non-positive length were removed.", sum(!ok)),
            call. = FALSE)
    data_end <- data_end[ok, , drop = FALSE]
  }

  # Unique observations (keeps the first occurrence)
  key <- paste(data_end$predator, data_end$prey,
               data_end$standardised_predator_length, data_end$si_prey_length)
  data_u <- data_end[!duplicated(key), , drop = FALSE]

  # Down-sample over-represented predator-prey pairs. Pairs are visited
  # predator by predator (in order of appearance), then prey by prey.
  pair <- paste(data_u$predator, data_u$prey, sep = "\r")
  visit_order <- order(match(data_u$predator, unique(data_u$predator)), seq_along(pair))
  rows_by_pair <- split(seq_along(pair), factor(pair, levels = unique(pair[visit_order])))

  drop_rows <- unlist(lapply(rows_by_pair, function(idx) {
    n <- length(idx)
    if (n > observation_resample) idx[sample.int(n, n - observation_resample)] else NULL
  }), use.names = FALSE)
  if (length(drop_rows) > 0) data_u <- data_u[-drop_rows, , drop = FALSE]

  data.frame(MPrey = log10(data_u$si_prey_length),
             MPred = log10(data_u$standardised_predator_length))
}

# Negative log-likelihood of the niche model (objective minimised by GenSA).
# mean_prey, sd_prey and pM don't depend on `pars`; metaweb_mod_parameters()
# passes them precomputed.
model <- function(pars, data,
                  mean_prey = mean(unique(data$MPrey)),
                  sd_prey = sd(unique(data$MPrey)),
                  pM = dnorm(data$MPrey, mean = mean_prey, sd = sd_prey)) {
  MPred <- data$MPred
  MPrey <- data$MPrey

  # Optimum and range
  o <- pars[1] + pars[2] * MPred
  r <- pars[3] + pars[4] * MPred

  # Conditional probability of interaction given prey size
  pLM <- exp(-(o - MPrey)^2 / 2 / r^2)

  # Integrated denominator
  pL <- r / (r^2 + sd_prey^2)^0.5 * exp(-(o - mean_prey)^2 / 2 / (r^2 + sd_prey^2))

  # Posterior probability
  pML <- pLM * pM / pL
  pML[pML <= 0] <- .Machine$double.xmin  # avoid log(0)

  -sum(log(pML))
}
