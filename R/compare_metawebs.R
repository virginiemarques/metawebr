#' Compare two metawebs for manual quality checks
#'
#' Lists the interactions that differ between two metawebs, e.g. before and
#' after [correct_metaweb_fish()], so you can check them by hand or export them
#' to a spreadsheet.
#'
#' When `after` was returned by [correct_metaweb_fish()] and `before` is the
#' metaweb it corrected, the correction columns give the last correction that
#' changed each link: `"herbivory"`, `"water_position"`, `"cannibalism"`,
#' `"depth"`, `"small_fish"`, `"producers"` or `"big_fish"`, and the source
#' columns the data it was based on (e.g. `"diet_items"` or `"trophic_level"`
#' for herbivory, `"DemersPelag"` for water position). They are `NA` for
#' unchanged links and for changes made by any other step, such as binarisation
#' or manual edits.
#'
#' Taxa present in only one of the metawebs, such as the producer nodes added
#' by [correct_metaweb_fish()], count as having links of 0 in the other.
#'
#' @param before,after Metawebs (prey in rows, predators in columns), with
#'   probabilities or 0/1 links.
#' @param data_traits Optional traits with taxon names as row names. Their
#'   `CommonLengthEstim`, `TrophicLevel`, `DemersPelag`, `Herbivore`,
#'   `Herbivore_source`, `DepthMin` and `DepthMax` columns, when present, are
#'   added to the table.
#' @param level `"pair"` (default) for one row per species pair, `"link"` for
#'   one row per changed link (each direction separately), or `"taxon"` for one
#'   row per taxon.
#' @param min_proba Pair and link levels: changed links below `min_proba` both
#'   before and after are not listed, to hide negligible probabilities. Use 0 to
#'   list every change.
#'
#' @return With `level = "pair"`, a data frame with one row per pair of taxa
#'   where at least one direction changed, each pair appearing once. Columns:
#'   `taxon_1`, `taxon_2`; `t1_eats_t2_before`, `t1_eats_t2_after`,
#'   `t1_eats_t2_correction`, `t1_eats_t2_source` for taxon 1 eating taxon 2;
#'   the same four
#'   `t2_eats_t1_*` columns for the other direction (`NA` when taxon 1 and 2 are
#'   the same, i.e. cannibalism); then the traits of both taxa (`t1_*`, `t2_*`).
#'   Rows are sorted by the highest probability, before correction, of their
#'   changed links.
#'
#'   With `level = "link"`, a data frame sorted by decreasing `before`, with the
#'   columns `prey`, `predator`, `before`, `after`, `change` (`"removed"`,
#'   `"added"` or `"modified"`), `correction` and `source`, then the traits of
#'   the prey
#'   (`prey_*`) and of the predator (`predator_*`).
#'
#'   With `level = "taxon"`, a data frame with one row per taxon: `taxon`,
#'   `prey_before`, `prey_after`, `predators_before` and `predators_after` (sums
#'   of link values: numbers of links for 0/1 metawebs, expected numbers for
#'   probabilities), one `prey_change_<correction>` column per correction (change
#'   in prey, negative when prey were lost; `prey_change_other` for changes
#'   without a recorded correction), then the traits. All links are used at this
#'   level, whatever `min_proba`.
#'
#' @examples
#' traits <- data.frame(
#'   CommonLengthEstim = c(100, 30, 5),
#'   TrophicLevel = c(4.1, 2.2, 3.3),
#'   DemersPelag = c("benthopelagic", "benthopelagic", "demersal"),
#'   row.names = c("Gadus_morhua", "Sarpa_salpa", "Gasterosteus_aculeatus")
#' )
#' mw <- apply_model_metaweb(traits, c(-0.56, 0.96, 0.12, 0.14))
#' mw_corr <- suppressMessages(correct_metaweb_fish(mw, traits))
#'
#' compare_metawebs(mw, mw_corr, traits)
#' compare_metawebs(mw, mw_corr, traits, level = "taxon")
#'
#' @export
compare_metawebs <- function(before, after, data_traits = NULL,
                             level = c("pair", "link", "taxon"), min_proba = 0.01) {
  level <- match.arg(level)
  log <- attr(after, "corrections")
  before <- check_metaweb(before, "before")
  after <- check_metaweb(after, "after")
  if (!is.null(data_traits)) check_traits(data_traits)
  if (!is.numeric(min_proba) || length(min_proba) != 1 || is.na(min_proba) || min_proba < 0) {
    stop("`min_proba` must be a single number >= 0.", call. = FALSE)
  }

  taxa <- union(rownames(after), rownames(before))
  B <- expand_metaweb(before, taxa)
  A <- expand_metaweb(after, taxa)

  switch(level,
         pair = pair_changes(taxa, B, A, log, data_traits, min_proba),
         link = link_changes(taxa, B, A, log, data_traits, min_proba),
         taxon = taxon_changes(taxa, B, A, link_changes(taxa, B, A, log, NULL, 0), data_traits))
}

# Square matrix over `taxa`, with 0 for taxa absent from `m`
expand_metaweb <- function(m, taxa) {
  out <- matrix(0, length(taxa), length(taxa), dimnames = list(taxa, taxa))
  out[rownames(m), colnames(m)] <- m
  out
}

# Links at row/column positions `prey` and `predator` of B (before) and A
# (after), with the recorded correction and whether they changed
links_at <- function(taxa, B, A, prey, predator, log) {
  cells <- cbind(prey, predator)
  links <- data.frame(prey = taxa[prey], predator = taxa[predator], before = B[cells], after = A[cells])
  matched <- match_corrections(links, log)
  links$correction <- matched$correction
  links$source <- matched$source
  links$changed <- abs(links$after - links$before) > 1e-12
  links
}

# One row per changed link
link_changes <- function(taxa, B, A, log, data_traits, min_proba) {
  idx <- which(abs(A - B) > 1e-12, arr.ind = TRUE)
  links <- links_at(taxa, B, A, idx[, 1], idx[, 2], log)
  links$change <- ifelse(links$before == 0, "added", ifelse(links$after == 0, "removed", "modified"))
  links <- links[links$before >= min_proba | links$after >= min_proba,
                 c("prey", "predator", "before", "after", "change", "correction", "source"), drop = FALSE]
  links <- links[order(-links$before, links$prey, links$predator), , drop = FALSE]
  rownames(links) <- NULL
  if (!is.null(data_traits)) {
    links <- cbind(links,
                   taxon_traits(data_traits, links$prey, "prey_"),
                   taxon_traits(data_traits, links$predator, "predator_"))
  }
  links
}

# One row per species pair with at least one changed direction
pair_changes <- function(taxa, B, A, log, data_traits, min_proba) {
  changed <- abs(A - B) > 1e-12
  idx <- which((changed | t(changed)) & upper.tri(changed, diag = TRUE), arr.ind = TRUE)
  i <- idx[, 1]
  j <- idx[, 2]  # i <= j: taxon_1 = taxa[i], taxon_2 = taxa[j]

  eats_12 <- links_at(taxa, B, A, prey = j, predator = i, log)  # taxon 1 eats taxon 2
  eats_21 <- links_at(taxa, B, A, prey = i, predator = j, log)  # taxon 2 eats taxon 1
  self <- i == j
  eats_21[self, c("before", "after")] <- NA
  eats_21$correction[self] <- NA
  eats_21$source[self] <- NA
  eats_21$changed[self] <- FALSE

  relevant_12 <- eats_12$changed & pmax(eats_12$before, eats_12$after) >= min_proba
  relevant_21 <- eats_21$changed & pmax(eats_21$before, eats_21$after) >= min_proba
  keep <- relevant_12 | relevant_21
  top <- pmax(ifelse(relevant_12, eats_12$before, -Inf), ifelse(relevant_21, eats_21$before, -Inf))

  out <- data.frame(taxon_1 = taxa[i], taxon_2 = taxa[j],
                    t1_eats_t2_before = eats_12$before, t1_eats_t2_after = eats_12$after,
                    t1_eats_t2_correction = eats_12$correction, t1_eats_t2_source = eats_12$source,
                    t2_eats_t1_before = eats_21$before, t2_eats_t1_after = eats_21$after,
                    t2_eats_t1_correction = eats_21$correction, t2_eats_t1_source = eats_21$source)
  out <- out[keep, , drop = FALSE]
  out <- out[order(-top[keep], out$taxon_1, out$taxon_2), , drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(data_traits)) {
    out <- cbind(out,
                 taxon_traits(data_traits, out$taxon_1, "t1_"),
                 taxon_traits(data_traits, out$taxon_2, "t2_"))
  }
  out
}

# Correction recorded by correct_metaweb_fish() for each link, and its source,
# used only when the values before and after match the recorded ones
match_corrections <- function(links, log) {
  none <- list(correction = rep(NA_character_, nrow(links)), source = rep(NA_character_, nrow(links)))
  if (!is.data.frame(log) || nrow(links) == 0) return(none)
  i <- match(paste(links$prey, links$predator, sep = "\r"), paste(log$prey, log$predator, sep = "\r"))
  same <- !is.na(i) &
    abs(links$before - log$before[i]) < 1e-12 &
    abs(links$after - log$after[i]) < 1e-12
  source <- if ("source" %in% names(log)) log$source[i] else NA_character_
  list(correction = ifelse(same, log$correction[i], NA_character_),
       source = ifelse(same, source, NA_character_))
}

# Traits of `taxa` (NA for taxa not in data_traits, e.g. producers)
taxon_traits <- function(data_traits, taxa, prefix) {
  cols <- intersect(c("CommonLengthEstim", "TrophicLevel", "DemersPelag", "Herbivore", "Herbivore_source",
                      "DepthMin", "DepthMax"), names(data_traits))
  out <- data_traits[match(taxa, rownames(data_traits)), cols, drop = FALSE]
  names(out) <- paste0(prefix, cols)
  rownames(out) <- NULL
  out
}

# One row per taxon: prey and predators before/after, change in prey per correction
taxon_changes <- function(taxa, B, A, links, data_traits) {
  out <- data.frame(taxon = taxa,
                    prey_before = unname(colSums(B)), prey_after = unname(colSums(A)),
                    predators_before = unname(rowSums(B)), predators_after = unname(rowSums(A)))

  step <- ifelse(is.na(links$correction), "other", links$correction)
  steps <- c("herbivory", "water_position", "cannibalism", "depth", "small_fish", "producers", "big_fish", "other")
  for (s in intersect(steps, step)) {
    sel <- step == s
    delta <- tapply(links$after[sel] - links$before[sel], factor(links$predator[sel], levels = taxa), sum)
    delta[is.na(delta)] <- 0
    out[[paste0("prey_change_", s)]] <- as.numeric(delta)
  }

  if (!is.null(data_traits)) out <- cbind(out, taxon_traits(data_traits, taxa, ""))
  out
}
