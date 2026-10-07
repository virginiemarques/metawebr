#' Set which taxa are treated as herbivores
#'
#' [correct_metaweb_fish()] treats herbivores differently: they get no fish
#' prey and feed on primary producers. `set_herbivores()` stores which taxa are
#' herbivores in a logical `Herbivore` column of the trait table, which
#' [correct_metaweb_fish()] uses when present, and where the decision came from
#' in a `Herbivore_source` column. The source is recorded for every link the
#' herbivory correction changes (see [compare_metawebs()]).
#'
#' The `Herbivore` column is built from:
#' * `source = "trophic_level"`: taxa with a trophic level below `threshold`
#'   (source `"trophic_level"`);
#' * `source = "fishbase"`: FishBase diet data, see
#'   [get_herbivores_fishbase()] (sources `"diet_items"` and
#'   `"feeding_type"`), and the trophic level for taxa without diet data
#'   (source `"trophic_level"`).
#'
#' With `source = NULL` (default), an existing `Herbivore` column is kept, and
#' a missing one is built from the trophic level. Then `add` and `remove` are
#' applied (source `"manual"`).
#'
#' @param data_traits Traits with taxon names as row names.
#' @param add Taxa to treat as herbivores.
#' @param remove Taxa not to treat as herbivores.
#' @param threshold Trophic level below which taxa are herbivores, used for
#'   `source = "trophic_level"` and for taxa without FishBase diet data.
#' @param source `NULL`, `"trophic_level"` or `"fishbase"`: how to build the
#'   `Herbivore` column (see Details).
#' @param ... Passed to [get_herbivores_fishbase()] when `source = "fishbase"`.
#'
#' @return `data_traits` with a logical `Herbivore` column and a
#'   `Herbivore_source` column.
#'
#' @examples
#' traits <- data.frame(
#'   TrophicLevel = c(2.0, 2.2, 3.4),
#'   row.names = c("Sarpa_salpa", "Scarus_iseri", "Coris_julis")
#' )
#' set_herbivores(traits, remove = "Scarus_iseri")
#' \dontrun{
#' set_herbivores(df_traits_obis, source = "fishbase")
#' }
#'
#' @export
set_herbivores <- function(data_traits, add = NULL, remove = NULL, threshold = 2.4,
                           source = NULL, ...) {
  check_traits(data_traits)
  if (!is.null(source)) source <- match.arg(source, c("trophic_level", "fishbase"))
  both <- intersect(add, remove)
  if (length(both) > 0) {
    stop(sprintf("Taxa both in `add` and `remove`: %s.", format_taxa(both)), call. = FALSE)
  }
  unknown <- setdiff(c(add, remove), rownames(data_traits))
  if (length(unknown) > 0) {
    stop(sprintf("Taxa not in the row names of `data_traits`: %s.", format_taxa(unknown)), call. = FALSE)
  }

  taxa <- rownames(data_traits)
  if (is.null(source) && "Herbivore" %in% names(data_traits)) {
    if (!"Herbivore_source" %in% names(data_traits)) {
      data_traits$Herbivore_source <- ifelse(is.na(data_traits$Herbivore), NA_character_, "Herbivore column")
    }
  } else if (identical(source, "fishbase")) {
    diet <- get_herbivores_fishbase(taxa, ...)
    data_traits$Herbivore <- diet$Herbivore
    data_traits$Herbivore_source <- diet$Herbivore_source
    no_diet <- is.na(data_traits$Herbivore)
    if ("TrophicLevel" %in% names(data_traits)) {
      tl <- data_traits$TrophicLevel
      fill <- no_diet & !is.na(tl)
      data_traits$Herbivore[fill] <- tl[fill] < threshold
      data_traits$Herbivore_source[fill] <- "trophic_level"
    }
    counts <- table(data_traits$Herbivore_source[data_traits$Herbivore %in% TRUE])
    message(sprintf("%d herbivores (%s).", sum(data_traits$Herbivore, na.rm = TRUE),
                    if (length(counts) > 0) paste(names(counts), counts, sep = ": ", collapse = ", ") else "none"))
  } else {
    check_traits(data_traits, "TrophicLevel")
    data_traits$Herbivore <- taxa %in% herbivore_taxa(data_traits[, "TrophicLevel", drop = FALSE], taxa, threshold)
    data_traits$Herbivore_source <- ifelse(is.na(data_traits$TrophicLevel), NA_character_, "trophic_level")
  }

  data_traits[add, "Herbivore"] <- TRUE
  data_traits[remove, "Herbivore"] <- FALSE
  data_traits[c(add, remove), "Herbivore_source"] <- "manual"
  data_traits
}

#' Identify herbivores from FishBase diet data
#'
#' Classifies taxa as herbivores (including detritivores) with two FishBase
#' tables, in this order:
#'
#' 1. **Diet items** (`fooditems` table): a taxon with at least `min_records`
#'    food-item records is a herbivore when the share of records that are
#'    `plant_items` (by default plants and detritus) is at least
#'    `min_plant_share` (source `"diet_items"`).
#' 2. **Feeding type** (`ecology` table): otherwise, a taxon is a herbivore
#'    when its main feeding type is one of `feeding_types`, and not a
#'    herbivore when it is one of `animal_feeding_types` (source
#'    `"feeding_type"`). By default, `"browsing on substrate"` counts as
#'    herbivory: it groups surgeonfishes, parrotfishes, rabbitfishes and
#'    damselfishes, and 75% of its species have a FishBase trophic level below
#'    2.4. Other feeding types (`"variable"`, `"other"`, ...) are not
#'    informative and get `NA`.
#'
#' Taxa without informative data get `NA`; [set_herbivores()] then uses their
#' trophic level. Genera and families pool the
#' records of all their FishBase species. Unlike the trophic level, these data
#' don't depend on a cut-off on a continuous estimate, so they misclassify
#' fewer omnivores and detritivores. Food-item records are counted, not
#' weighted by diet volume.
#'
#' @param taxa Taxon names: species as `Genus_species`, genera or families.
#' @param min_plant_share Minimum share of plant and detritus food items.
#' @param min_records Minimum number of food-item records to use them.
#' @param plant_items Values of the FishBase `FoodI` column counted as plant
#'   food.
#' @param feeding_types Values of the FishBase `FeedingType` column counted as
#'   herbivory.
#' @param animal_feeding_types Values of `FeedingType` counted as not
#'   herbivory.
#'
#' @return A data frame with taxon names as row names and the columns `taxon`,
#'   `n_diet_records`, `plant_share`, `feeding_type`, `Herbivore` (logical) and
#'   `Herbivore_source`.
#'
#' @seealso [set_herbivores()], which uses it with `source = "fishbase"`.
#'
#' @export
get_herbivores_fishbase <- function(taxa, min_plant_share = 0.5, min_records = 3,
                                    plant_items = c("plants", "detritus"),
                                    feeding_types = c("grazing on aquatic plants", "browsing on substrate"),
                                    animal_feeding_types = c("hunting macrofauna (predator)",
                                                             "filtering plankton",
                                                             "selective plankton feeding",
                                                             "feeding on a host (parasite)",
                                                             "picking parasites off a host (cleaner)",
                                                             "feeding on dead animals (scavenger)",
                                                             "feeding on the prey of a host (commensal)")) {
  check_count(min_records, "min_records")
  if (!is.numeric(min_plant_share) || length(min_plant_share) != 1 || is.na(min_plant_share) ||
      min_plant_share < 0 || min_plant_share > 1) {
    stop("`min_plant_share` must be a single number between 0 and 1.", call. = FALSE)
  }
  taxa <- unique(gsub(" ", "_", taxa))
  codes <- taxon_spec_codes(taxa)
  food <- fb_table("fooditems")[, c("SpecCode", "FoodI")]
  eco <- fb_table("ecology")[, c("SpecCode", "FeedingType")]

  n_records <- vapply(codes, function(cc) sum(food$SpecCode %in% cc), integer(1))
  n_plant <- vapply(codes, function(cc) sum(food$SpecCode %in% cc & food$FoodI %in% plant_items), integer(1))
  feeding <- vapply(codes, function(cc) mode_or_na(eco$FeedingType[eco$SpecCode %in% cc]), character(1))

  share <- ifelse(n_records > 0, n_plant / n_records, NA_real_)
  use_diet <- n_records >= min_records
  informative <- feeding %in% c(feeding_types, animal_feeding_types)
  herb <- ifelse(use_diet, share >= min_plant_share,
                 ifelse(informative, feeding %in% feeding_types, NA))
  src <- ifelse(use_diet, "diet_items", ifelse(informative, "feeding_type", NA_character_))

  data.frame(taxon = taxa, n_diet_records = unname(n_records), plant_share = unname(share),
             feeding_type = unname(feeding), Herbivore = unname(herb), Herbivore_source = unname(src),
             row.names = taxa)
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

# Where the herbivore status of each taxon comes from, named by taxon
herbivore_sources <- function(data_traits, taxa) {
  src <- rep("trophic_level", length(taxa))
  if ("Herbivore" %in% names(data_traits)) {
    flag <- as.logical(data_traits[taxa, "Herbivore"])
    given <- if ("Herbivore_source" %in% names(data_traits)) {
      as.character(data_traits[taxa, "Herbivore_source"])
    } else {
      rep("Herbivore column", length(taxa))
    }
    src <- ifelse(is.na(flag), src, ifelse(is.na(given), "Herbivore column", given))
  }
  setNames(src, taxa)
}
