#' Get fish traits from FishBase
#'
#' Downloads from FishBase (with \pkg{rfishbase}) the traits needed to build a
#' metaweb, for the species of a presence table:
#'
#' * `CommonLengthEstim`: common length (cm), or 0.6 x maximum length when
#'   FishBase has no common length;
#' * `TrophicLevel`: FishBase trophic level estimate;
#' * `DemersPelag`: position in the water column.
#'
#' @param data_presence Data frame with a column of species names, as
#'   `Genus_species` (spaces are converted to underscores).
#' @param column_species Name of that column.
#'
#' @return A data frame with one row per species and species names as row
#'   names: a column named after `column_species`, then the three traits.
#'   Species not found in FishBase have missing traits, with a warning.
#'
#' @seealso [get_traits_higher_edna()] and [infer_traits()] for genus- or
#'   family-level taxa, [clean_traits()] to fill missing values.
#'
#' @examples
#' \dontrun{
#' presence <- data.frame(species = c("Gadus_morhua", "Clupea_harengus"))
#' get_traits(presence)
#' }
#'
#' @export
get_traits <- function(data_presence, column_species = "species") {
  species <- presence_taxa(data_presence, column_species)
  info <- fishbase_traits_table()

  idx <- match(species, info$Species)
  warn_not_found(species[is.na(idx)])

  out <- data.frame(species, info[idx, c("CommonLengthEstim", "TrophicLevel", "DemersPelag")],
                    row.names = species)
  names(out)[1] <- column_species
  out
}

#' Get fish traits for species- and genus-level eDNA assignments
#'
#' Like [get_traits()], for eDNA data where some sequences are assigned to a
#' species and others only to a genus. Species get their FishBase traits;
#' genera get the mean length and trophic level, and the most common water
#' position, of all FishBase species of the genus.
#'
#' @param data_presence Data frame of detections.
#' @param column_species Column with species-level assignments (`Genus_species`).
#' @param column_taxon Column with the assigned taxon (species or genus). Taxa
#'   that are not already in `column_species` are treated as genera.
#'
#' @return A data frame with taxon names as row names and the columns `taxon`,
#'   `CommonLengthEstim`, `TrophicLevel` and `DemersPelag`.
#'
#' @seealso [infer_traits()], which also handles family-level taxa.
#'
#' @export
get_traits_higher_edna <- function(data_presence, column_species = "Species", column_taxon = "taxon") {
  species <- presence_taxa(data_presence, column_species)
  genera <- setdiff(presence_taxa(data_presence, column_taxon), species)
  info <- fishbase_traits_table()

  species_traits <- aggregate_traits(info, species, "Species")
  genus_traits <- aggregate_traits(info, genera, "Genus")
  warn_not_found(c(species[species_traits$n_species == 0], genera[genus_traits$n_species == 0]))

  # Genera absent from FishBase are dropped, as before
  genus_traits <- genus_traits[genus_traits$n_species > 0, , drop = FALSE]
  rbind(species_traits, genus_traits)[, c("taxon", "CommonLengthEstim", "TrophicLevel", "DemersPelag")]
}

# Unique, non-empty taxon names of a column, with underscores instead of spaces
presence_taxa <- function(data_presence, column) {
  if (!is.data.frame(data_presence)) {
    stop("`data_presence` must be a data frame.", call. = FALSE)
  }
  if (!column %in% names(data_presence)) {
    stop(sprintf("`data_presence` has no column `%s`.", column), call. = FALSE)
  }
  taxa <- gsub(" ", "_", as.character(data_presence[[column]]))
  unique(taxa[!is.na(taxa) & taxa != ""])
}

warn_not_found <- function(taxa) {
  if (length(taxa) > 0) {
    warning(sprintf("%d taxa were not found in FishBase and have missing traits: %s.",
                    length(taxa), format_taxa(taxa)), call. = FALSE)
  }
}

# One row per FishBase species with taxonomy and the traits used by the package
# (tables read through the cache, see fb_table())
fishbase_traits_table <- function() {
  taxa <- fb_table("taxa")[, c("SpecCode", "Species", "Genus", "Family")]
  species <- fb_table("species")[, c("SpecCode", "CommonLength", "Length", "DemersPelag",
                                     "DepthRangeShallow", "DepthRangeDeep")]
  estimate <- fb_table("estimate")[, c("SpecCode", "Troph")]

  info <- merge(taxa, species, by = "SpecCode", all.x = TRUE)
  info <- merge(info, estimate, by = "SpecCode", all.x = TRUE)

  data.frame(SpecCode = info$SpecCode,
             Species = gsub(" ", "_", info$Species),
             Genus = info$Genus,
             Family = info$Family,
             CommonLengthEstim = ifelse(is.na(info$CommonLength), 0.6 * info$Length, info$CommonLength),
             TrophicLevel = info$Troph,
             DemersPelag = info$DemersPelag,
             DepthMin = info$DepthRangeShallow,
             DepthMax = info$DepthRangeDeep)
}

# FishBase species codes making up each taxon (species, genus or family name),
# as a list named by taxon; empty for taxa unknown to FishBase
taxon_spec_codes <- function(taxa, info = fishbase_traits_table()) {
  taxa <- gsub(" ", "_", taxa)
  lapply(setNames(taxa, taxa), function(x) {
    codes <- info$SpecCode[info$Species == x]
    if (length(codes) == 0) codes <- info$SpecCode[info$Genus %in% x]
    if (length(codes) == 0) codes <- info$SpecCode[info$Family %in% x]
    codes
  })
}

# Traits of `taxa` at a taxonomic `level` (a column of `info`): mean length and
# trophic level, most common water position, and number of species averaged.
aggregate_traits <- function(info, taxa, level) {
  if (length(taxa) == 0) {
    return(data.frame(taxon = character(), n_species = integer(), CommonLengthEstim = numeric(),
                      TrophicLevel = numeric(), DemersPelag = character()))
  }
  rows <- info[info[[level]] %in% taxa, , drop = FALSE]
  groups <- split(rows, factor(rows[[level]], levels = taxa))
  data.frame(
    taxon = taxa,
    n_species = vapply(groups, nrow, integer(1)),
    CommonLengthEstim = vapply(groups, function(g) mean_or_na(g$CommonLengthEstim), numeric(1)),
    TrophicLevel = vapply(groups, function(g) mean_or_na(g$TrophicLevel), numeric(1)),
    DemersPelag = vapply(groups, function(g) mode_or_na(g$DemersPelag), character(1)),
    row.names = taxa
  )
}

mean_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0) mean(x) else NA_real_
}

median_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0) median(x) else NA_real_
}

# Most common value; ties go to the first in alphabetical order
mode_or_na <- function(x) {
  x <- as.character(x[!is.na(x)])
  if (length(x) > 0) names(which.max(table(x))) else NA_character_
}
