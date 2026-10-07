#' Cross-validate the body-size niche model
#'
#' Measures how well the niche model predicts the prey of predators in studies
#' it was **not** calibrated on. The interaction records are split by study
#' (`group`, by default the `reference` column of [df_interaction_fish]). For
#' each study in turn, the model is calibrated with [metaweb_mod_parameters()]
#' on all other studies, then used to predict the links of the left-out study.
#'
#' Predictions are made as for a metaweb: one size per species, the mean
#' log10 length of each predator in the left-out study and of each prey. Every
#' predator of the left-out study is paired with every prey species of the
#' pool (`prey_pool`): pairs recorded in the data are links, the others are
#' taken as non-links. These are pseudo-absences: a prey not found in a
#' predator's stomachs may still be eaten, so specificity, TSS and AUC are
#' conservative.
#'
#' * `prey_pool = "all"` (default): all prey species of `data`, with their
#'   mean size over all studies. Every study can be scored, including those
#'   on a single predator. Some of these prey may not live where the predator
#'   was sampled, so this tests whether body size separates a predator's prey
#'   from other prey, not its local diet.
#' * `prey_pool = "study"`: only the prey recorded in the left-out study,
#'   with their sizes in that study. Closer to a local food web, but studies
#'   on a single predator have no non-links and can't be scored.
#'
#' Leave-one-study-out is stricter than random folds: records of the same study
#' share methods, region and species, so random folds would put near-copies of
#' the test records in the calibration data. Studies contributing fewer than
#' `min_links` links are not scored (their records are still used for
#' calibration in the other folds).
#'
#' With `in_sample = TRUE`, the same pairs are also predicted with the model
#' calibrated on all studies. The difference between in-sample and
#' out-of-sample scores shows how much the calibration overfits.
#'
#' Calibration is run once per study (plus once for `in_sample`), each for up
#' to `max_time` seconds. Call [set.seed()] first for reproducible results.
#'
#' @inheritParams metaweb_mod_parameters
#' @param data `NULL` (default) to use [df_interaction_fish], the path to a
#'   CSV file, or a data frame with the columns `predator`, `prey`,
#'   `standardised_predator_length`, `si_prey_length` and `group`.
#' @param group Column defining the folds: one study (or region, or predator)
#'   per value.
#' @param min_links Minimum number of observed links (distinct predator-prey
#'   species pairs) for a left-out study to be scored.
#' @param in_sample Also score the pairs with the model calibrated on all
#'   studies?
#' @param prey_pool `"all"` or `"study"`: prey paired with each left-out
#'   predator (see Details).
#'
#' @return A list with
#'   * `folds`: one row per study: `group`, `n_train` and `n_test` (records),
#'     `n_predators`, `n_prey`, `n_links`, `n_pairs`, `AUC`, `TSS` (maximum
#'     over thresholds), `threshold` (giving that TSS), `AUC_in_sample`, and
#'     the fitted parameters `a0`, `a1`, `b0`, `b1`. Unscored studies have `NA`
#'     scores;
#'   * `overall`: AUC, maximum TSS and its threshold over the pairs of all
#'     scored studies pooled, out of sample and (if `in_sample`) in sample;
#'   * `pairs`: every scored pair: `group`, `predator`, `prey`, `observed`
#'     (0/1), `proba` (out of sample) and `proba_in_sample`.
#'
#' @seealso [validate_metaweb()] to validate a corrected metaweb against
#'   independent links.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' cv <- cv_size_model(max_time = 30)
#' cv$folds
#' cv$overall
#' }
#'
#' @export
cv_size_model <- function(data = NULL, group = "reference", observation_resample = 50,
                          max_time = 60, min_links = 5, in_sample = TRUE,
                          prey_pool = c("all", "study")) {
  prey_pool <- match.arg(prey_pool)
  data <- read_interaction_data(data)
  if (!group %in% names(data)) stop(sprintf("`data` has no column `%s`.", group), call. = FALSE)
  check_count(min_links, "min_links")
  pred_len <- data$standardised_predator_length
  prey_len <- data$si_prey_length
  ok <- is.finite(pred_len) & pred_len > 0 & is.finite(prey_len) & prey_len > 0 & !is.na(data[[group]])
  data <- data[ok, , drop = FALSE]
  groups <- unique(as.character(data[[group]]))
  if (length(groups) < 2) stop("`data` must contain at least two groups.", call. = FALSE)

  fit <- function(d) {
    metaweb_mod_parameters(d, observation_resample = observation_resample,
                           max_time = max_time, verbose = FALSE)$calibrated_data
  }
  pars_all <- if (in_sample) fit(data)

  folds <- vector("list", length(groups))
  pairs <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    in_test <- data[[group]] == groups[k]
    test <- data[in_test, , drop = FALSE]
    p <- if (prey_pool == "all") species_pairs(test, prey_data = data, link_data = data) else species_pairs(test)
    n_links <- sum(p$observed)
    row <- data.frame(group = groups[k], n_train = sum(!in_test), n_test = nrow(test),
                      n_predators = length(unique(p$predator)), n_prey = length(unique(p$prey)),
                      n_links = n_links, n_pairs = nrow(p), AUC = NA_real_, TSS = NA_real_,
                      threshold = NA_real_, AUC_in_sample = NA_real_,
                      a0 = NA_real_, a1 = NA_real_, b0 = NA_real_, b1 = NA_real_)
    if (n_links >= min_links && n_links < nrow(p)) {
      pars <- fit(data[!in_test, , drop = FALSE])
      p$proba <- pLMFitted(p$prey_size, p$predator_size, pars)
      p$proba_in_sample <- if (in_sample) pLMFitted(p$prey_size, p$predator_size, pars_all) else NA_real_
      best <- best_tss(p$proba, p$observed)
      row[, c("AUC", "TSS", "threshold")] <- list(auc(p$proba, p$observed), best[["TSS"]], best[["threshold"]])
      if (in_sample) row$AUC_in_sample <- auc(p$proba_in_sample, p$observed)
      row[, c("a0", "a1", "b0", "b1")] <- as.list(pars)
      pairs[[k]] <- data.frame(group = groups[k], p[, c("predator", "prey", "observed", "proba", "proba_in_sample")])
    }
    folds[[k]] <- row
  }
  folds <- do.call(rbind, folds)
  pairs <- do.call(rbind, pairs)
  if (is.null(pairs)) stop("No group has enough links to be scored; lower `min_links`.", call. = FALSE)

  best <- best_tss(pairs$proba, pairs$observed)
  overall <- data.frame(model = "out_of_sample", AUC = auc(pairs$proba, pairs$observed),
                        TSS = best[["TSS"]], threshold = best[["threshold"]])
  if (in_sample) {
    best <- best_tss(pairs$proba_in_sample, pairs$observed)
    overall <- rbind(overall, data.frame(model = "in_sample", AUC = auc(pairs$proba_in_sample, pairs$observed),
                                         TSS = best[["TSS"]], threshold = best[["threshold"]]))
  }
  rownames(pairs) <- NULL
  list(folds = folds, overall = overall, pairs = pairs)
}

#' Validate a metaweb against observed links
#'
#' Compares a metaweb (probabilities or 0/1 links) with independent
#' observations of who eats whom, e.g. FishBase diet records
#' ([get_diet_links()]), GloBI, or your own stomach-content data.
#'
#' Without an `interaction` column, `links` lists observed links only. Each
#' predator with at least one observed prey in the metaweb is then paired with
#' every taxon of the metaweb: pairs not in `links` are taken as non-links
#' (pseudo-absences; diet studies rarely list every prey, so specificity is
#' conservative). With an `interaction` column (1 = link, 0 = no link), only
#' the listed pairs are used.
#'
#' Producer nodes are ignored. Taxa are matched by name, with spaces or
#' underscores.
#'
#' Validate against data that were not used to build the metaweb: if the
#' herbivores were set from FishBase diet data ([set_herbivores()] with
#' `source = "fishbase"`), validating against FishBase diet links favours the
#' corrected metaweb, and a warning says so.
#'
#' @param MW Metaweb (prey in rows, predators in columns), with probabilities
#'   or 0/1 links.
#' @param links Data frame with the columns `predator` and `prey`, and
#'   optionally `interaction` (0/1).
#' @param threshold Probability threshold for the presence/absence scores.
#'   `NULL` (default): the threshold maximising TSS for a probability metaweb,
#'   0.5 for a 0/1 metaweb.
#' @param exclude_self Ignore cannibalism (a taxon eating itself)?
#'
#' @return A list with
#'   * `summary`: one row with `n_predators`, `n_links`, `n_pairs`, `AUC`
#'     (`NA` for a 0/1 metaweb), `threshold`, `sensitivity` (share of observed
#'     links predicted), `specificity` (share of non-links not predicted),
#'     `TSS` (sensitivity + specificity - 1) and `precision` (share of
#'     predicted links that are observed);
#'   * `per_predator`: one row per predator: `predator`, `n_observed`,
#'     `n_predicted` (predicted prey among the tested pairs), `n_hits`,
#'     `sensitivity` and `missed` (observed prey not predicted);
#'   * `pairs`: every tested pair with `observed` and `predicted` values.
#'
#' @seealso [cv_size_model()] to cross-validate the size model, [get_proba()]
#'   with `method = "tss"` to choose the threshold from observed links.
#'
#' @examples
#' web <- matrix(c(0, 0.9, 0.2,
#'                 0, 0, 0.8,
#'                 0, 0, 0), 3, byrow = TRUE,
#'               dimnames = rep(list(c("A", "B", "C")), 2))
#' obs <- data.frame(predator = c("B", "C"), prey = c("A", "B"))
#' validate_metaweb(web, obs)$summary
#'
#' @export
validate_metaweb <- function(MW, links, threshold = NULL, exclude_self = TRUE) {
  MW <- check_metaweb(MW)
  pairs <- validation_pairs(MW, links, exclude_self)
  if (identical(attr(links, "source"), "fishbase_fooditems") &&
      "diet_items" %in% attr(MW, "data_sources")) {
    warning(paste("The herbivores of `MW` were set from FishBase diet items, and `links` come from",
                  "FishBase diet items too: the validation is not independent."), call. = FALSE)
  }

  binary <- all(MW %in% c(0, 1))
  if (is.null(threshold)) {
    threshold <- if (binary) 0.5 else best_tss(pairs$predicted, pairs$observed)[["threshold"]]
  }
  scores <- confusion_scores(pairs$predicted >= threshold, pairs$observed == 1)
  summary <- data.frame(n_predators = length(unique(pairs$predator[pairs$observed == 1])),
                        n_links = sum(pairs$observed), n_pairs = nrow(pairs),
                        AUC = if (binary) NA_real_ else auc(pairs$predicted, pairs$observed),
                        threshold = threshold, t(scores))

  hit <- pairs$predicted >= threshold
  per_predator <- do.call(rbind, lapply(split(seq_len(nrow(pairs)), pairs$predator), function(i) {
    obs <- pairs$observed[i] == 1
    data.frame(predator = pairs$predator[i[1]], n_observed = sum(obs), n_predicted = sum(hit[i]),
               n_hits = sum(obs & hit[i]),
               sensitivity = if (any(obs)) sum(obs & hit[i]) / sum(obs) else NA_real_,
               missed = paste(pairs$prey[i][obs & !hit[i]], collapse = ", "))
  }))
  per_predator <- per_predator[order(per_predator$sensitivity, -per_predator$n_observed), , drop = FALSE]
  rownames(per_predator) <- NULL
  list(summary = summary, per_predator = per_predator, pairs = pairs)
}

#' Fish-fish links from FishBase diet records
#'
#' Lists the predator-prey links among `taxa` recorded in the FishBase
#' food-item table, where the prey is identified to species. Genus- and
#' family-level taxa get the links of all their FishBase species. Use the
#' result with [validate_metaweb()] or [get_proba()] (`method = "tss"`).
#'
#' FishBase lists the prey species of relatively few fish (about 4,500
#' records), so many predators will have no link: they are simply not
#' tested by [validate_metaweb()].
#'
#' @param taxa Taxon names (`Genus_species`, genera or families), e.g. the row
#'   names of a metaweb.
#'
#' @return A data frame with the columns `predator`, `prey` and `n_records`,
#'   with the attribute `source = "fishbase_fooditems"`.
#'
#' @examples
#' \dontrun{
#' links <- get_diet_links(rownames(df_traits_obis))
#' }
#'
#' @export
get_diet_links <- function(taxa) {
  taxa <- unique(gsub(" ", "_", setdiff(taxa, producer_names)))
  codes <- taxon_spec_codes(taxa)
  code_map <- data.frame(taxon = rep(names(codes), lengths(codes)),
                         SpecCode = unlist(codes, use.names = FALSE))
  food <- fb_table("fooditems")
  food <- food[!is.na(food$PreySpecCode), c("SpecCode", "PreySpecCode")]

  m <- merge(food, setNames(code_map, c("predator", "SpecCode")), by = "SpecCode")
  m <- merge(m, setNames(code_map, c("prey", "PreySpecCode")), by = "PreySpecCode")
  out <- if (nrow(m) == 0) {
    data.frame(predator = character(), prey = character(), n_records = integer())
  } else {
    stats::aggregate(list(n_records = rep(1L, nrow(m))), m[, c("predator", "prey")], sum)
  }
  out <- out[order(out$predator, out$prey), , drop = FALSE]
  rownames(out) <- NULL
  message(sprintf("%d links for %d predators found in FishBase diet records.",
                  nrow(out), length(unique(out$predator))))
  attr(out, "source") <- "fishbase_fooditems"
  out
}

# Interaction data from NULL (package data), a CSV path or a data frame
read_interaction_data <- function(data_path) {
  if (is.null(data_path)) {
    metawebr::df_interaction_fish
  } else if (is.data.frame(data_path)) {
    data_path
  } else {
    read.csv(data_path)
  }
}

# All pairs of the predators of `d` with the prey of `prey_data`, with mean
# log10 sizes, and whether the pair is recorded in `link_data`
species_pairs <- function(d, prey_data = d, link_data = d) {
  pred_size <- tapply(log10(d$standardised_predator_length), d$predator, mean)
  prey_size <- tapply(log10(prey_data$si_prey_length), prey_data$prey, mean)
  p <- expand.grid(predator = names(pred_size), prey = names(prey_size), stringsAsFactors = FALSE)
  p <- p[p$predator != p$prey, , drop = FALSE]
  p$predator_size <- unname(pred_size[p$predator])
  p$prey_size <- unname(prey_size[p$prey])
  p$observed <- as.integer(paste(p$predator, p$prey) %in% paste(link_data$predator, link_data$prey))
  p
}

# Pairs of the metaweb to test, with observed (0/1) and predicted values
validation_pairs <- function(MW, links, exclude_self = TRUE) {
  if (!is.data.frame(links) || !all(c("predator", "prey") %in% names(links))) {
    stop("`links` must be a data frame with columns `predator` and `prey`.", call. = FALSE)
  }
  taxa <- setdiff(rownames(MW), producer_names)
  links$predator <- gsub(" ", "_", as.character(links$predator))
  links$prey <- gsub(" ", "_", as.character(links$prey))
  inside <- links$predator %in% taxa & links$prey %in% taxa
  if (exclude_self) inside <- inside & links$predator != links$prey
  if (any(!inside)) {
    message(sprintf("%d links with a taxon absent from the metaweb (or cannibalism) ignored.", sum(!inside)))
  }
  links <- links[inside, , drop = FALSE]

  if ("interaction" %in% names(links)) {
    pairs <- data.frame(predator = links$predator, prey = links$prey, observed = as.integer(links$interaction))
    pairs <- pairs[!duplicated(pairs[, c("predator", "prey")]), , drop = FALSE]
  } else {
    predators <- unique(links$predator)
    pairs <- expand.grid(prey = taxa, predator = predators, stringsAsFactors = FALSE)[, c("predator", "prey")]
    if (exclude_self) pairs <- pairs[pairs$predator != pairs$prey, , drop = FALSE]
    pairs$observed <- as.integer(paste(pairs$predator, pairs$prey) %in% paste(links$predator, links$prey))
  }
  if (sum(pairs$observed) == 0) stop("No observed link involves two taxa of the metaweb.", call. = FALSE)
  pairs$predicted <- unclass(MW)[cbind(match(pairs$prey, rownames(MW)), match(pairs$predator, colnames(MW)))]
  rownames(pairs) <- NULL
  pairs
}

# Area under the ROC curve (Mann-Whitney), NA without both classes
auc <- function(score, observed) {
  pos <- observed == 1
  n1 <- sum(pos)
  n0 <- sum(!pos)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

confusion_scores <- function(predicted, observed) {
  tp <- sum(predicted & observed)
  fn <- sum(!predicted & observed)
  tn <- sum(!predicted & !observed)
  fp <- sum(predicted & !observed)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  c(sensitivity = sens, specificity = spec, TSS = sens + spec - 1,
    precision = if (tp + fp > 0) tp / (tp + fp) else NA_real_)
}

# Thresholds tried by get_proba() and the TSS-based functions
threshold_grid <- c(seq(0.1, 0.975, 0.025), 0.99, 0.995)

# TSS at each threshold of `grid`, and the threshold with the highest TSS
tss_curve <- function(score, observed, grid = threshold_grid) {
  vapply(grid, function(th) confusion_scores(score >= th, observed == 1)[["TSS"]], numeric(1))
}

best_tss <- function(score, observed, grid = threshold_grid) {
  tss <- tss_curve(score, observed, grid)
  if (all(is.na(tss))) return(c(TSS = NA_real_, threshold = NA_real_))
  c(TSS = max(tss, na.rm = TRUE), threshold = grid[which.max(tss)])
}
