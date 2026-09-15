#' Set which taxa are treated as herbivores
#'
#' By default, [correct_metaweb_fish()] treats taxa with a trophic level below
#' `threshold` as herbivores: they get no fish prey and feed on primary
#' producers. Trophic levels from FishBase can misclassify some taxa, so
#' `set_herbivores()` stores the choice in a logical `Herbivore` column of the
#' trait table, which [correct_metaweb_fish()] uses when present.
#'
#' If `data_traits` has no `Herbivore` column yet, it is first created from the
#' trophic level (`TrophicLevel < threshold`).
#'
#' @param data_traits Traits with taxon names as row names.
#' @param add Taxa to treat as herbivores.
#' @param remove Taxa not to treat as herbivores.
#' @param threshold Trophic level below which taxa are herbivores, used when
#'   creating the `Herbivore` column.
#'
#' @return `data_traits` with a logical `Herbivore` column.
#'
#' @examples
#' traits <- data.frame(
#'   TrophicLevel = c(2.0, 2.2, 3.4),
#'   row.names = c("Sarpa_salpa", "Scarus_iseri", "Coris_julis")
#' )
#' set_herbivores(traits, remove = "Scarus_iseri")
#'
#' @export
set_herbivores <- function(data_traits, add = NULL, remove = NULL, threshold = 2.4) {
  check_traits(data_traits)
  both <- intersect(add, remove)
  if (length(both) > 0) {
    stop(sprintf("Taxa both in `add` and `remove`: %s.", format_taxa(both)), call. = FALSE)
  }
  unknown <- setdiff(c(add, remove), rownames(data_traits))
  if (length(unknown) > 0) {
    stop(sprintf("Taxa not in the row names of `data_traits`: %s.", format_taxa(unknown)), call. = FALSE)
  }

  if (!"Herbivore" %in% names(data_traits)) {
    check_traits(data_traits, "TrophicLevel")
    data_traits$Herbivore <- rownames(data_traits) %in%
      herbivore_taxa(data_traits, rownames(data_traits), threshold)
  }
  data_traits[add, "Herbivore"] <- TRUE
  data_traits[remove, "Herbivore"] <- FALSE
  data_traits
}

# Taxa (among `taxa`) treated as herbivores: the `Herbivore` column when
# present and not NA, otherwise a trophic level below `threshold`.
herbivore_taxa <- function(data_traits, taxa, threshold = 2.4) {
  tl <- if ("TrophicLevel" %in% names(data_traits)) data_traits[taxa, "TrophicLevel"] else rep(NA_real_, length(taxa))
  is_herb <- !is.na(tl) & tl < threshold
  if ("Herbivore" %in% names(data_traits)) {
    flag <- as.logical(data_traits[taxa, "Herbivore"])
    is_herb <- ifelse(is.na(flag), is_herb, flag)
  }
  taxa[is_herb]
}
