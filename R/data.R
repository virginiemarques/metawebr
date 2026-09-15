#' Predator-prey body-size records for fish
#'
#' Observed feeding events between fish predators and their prey, with body
#' lengths. Used by [metaweb_mod_parameters()] to calibrate the niche model.
#'
#' @format A data frame with 5847 rows and 61 columns. The package uses four of
#'   them: `predator` and `prey` (scientific names),
#'   `standardised_predator_length` (predator length, cm) and `si_prey_length`
#'   (prey length, cm).
"df_interaction_fish"

#' Depth range and water position of marine fish species
#'
#' @format A data frame with 11786 rows and 4 columns:
#' \describe{
#'   \item{Species}{Scientific name.}
#'   \item{Depth_min, Depth_max}{Depth range (m).}
#'   \item{Env_2}{Position in the water column (FishBase categories).}
#' }
"df_depth_env"

#' Fish occurrences at Mediterranean and Atlantic reef sites
#'
#' Example community data from OBIS: the 40 fish species (bony fish and
#' elasmobranchs) with the most occurrence records since 2000 at three
#' Mediterranean rocky-reef sites (Marseille, Medes Islands, Port-Cros) and
#' three tropical western Atlantic coral-reef sites (Florida Keys, Belize Barrier
#' Reef, Bonaire). Each site is a bounding box a few tens of kilometres wide.
#' Built by `data-raw/df_presence_obis.R`.
#'
#' @format A data frame with 240 rows and 6 columns:
#' \describe{
#'   \item{station}{Site name.}
#'   \item{region}{`"Mediterranean"` or `"Atlantic_reef"`.}
#'   \item{lon, lat}{Centre of the site bounding box (decimal degrees).}
#'   \item{species}{Species name, as `Genus_species`.}
#'   \item{count}{Number of OBIS occurrence records since 2000.}
#' }
#' @source OBIS (2026). Ocean Biodiversity Information System.
#'   Intergovernmental Oceanographic Commission of UNESCO. obis.org, accessed
#'   14 September 2026.
#' @seealso [df_traits_obis] for the traits of these species.
"df_presence_obis"

#' FishBase traits of the species in df_presence_obis
#'
#' Output of `get_traits(df_presence_obis, column_species = "species")`, stored
#' so that examples, tests and the vignette run offline. One species,
#' *Equetus punctatus*, was not found in FishBase and has missing traits.
#'
#' @format A data frame with 136 rows, species names as row names, and 4
#'   columns:
#' \describe{
#'   \item{species}{Species name.}
#'   \item{CommonLengthEstim}{Common length (cm), or 0.6 x maximum length.}
#'   \item{TrophicLevel}{FishBase trophic level estimate.}
#'   \item{DemersPelag}{Position in the water column.}
#' }
#' @source FishBase (Froese, R. & Pauly, D., eds.), retrieved with rfishbase
#'   5.0.1 on 14 September 2026.
"df_traits_obis"
