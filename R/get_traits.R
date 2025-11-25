#' Title: get_traits
#'
#' Find the parameter for the metaweb
#'
#' @description
#'
#'# End goal:
# 1. taxa
# 2. Common length
# 3. trophic_level
# 4. Env_2
#'
#' @param x Description of the first parameter.
#' @param y Description of the second parameter (if applicable).
#' @param ... Other optional parameters passed to methods.
#'
#' @details
#'
#' @return
#' Description of the object that the function returns.
#' If the function doesn't return anything meaningful, you can say `NULL`.
#'
#' @examples
#'
#'
#' @import rfishbase
#' @import dplyr
#'
#' @export
#'

get_traits <- function(data_presence, column_species = "species"){

  # Get all fishbase - first fetch all species names
  all_fishbase <- rfishbase::load_taxa()
  # Then get info on species table
  estimate_table <- rfishbase::estimate()
  species_table <- rfishbase::species() # Env_2 = DemersPelag
  info_fishbase <- species_table |>
    left_join(all_fishbase[,c("SpecCode", "Species")], by = c("SpecCode")) |>
    # Add estimate table
    left_join(estimate_table[,c("SpecCode", "Troph")]) |>
    # Add a common_length_estimated with 0.60 of length if NA
    mutate(CommonLengthEstim = ifelse(is.na(CommonLength), 0.6*Length, CommonLength)) |>
    mutate(Species = gsub(" ", "_", Species)) |>
    dplyr::select(Species, CommonLengthEstim, TrophicLevel = Troph, DemersPelag)

  # Intersect with data_presence
  data_traits <- data.frame(species = data_presence[,column_species])
  data_traits <- data_traits |>   distinct(!!sym(column_species)) |> filter(!!sym(column_species) != "")
  # Add a test here

  data_traits_completed <- data_traits |>
    left_join(info_fishbase, by = setNames("Species", column_species))

  rownames(data_traits_completed) <- data_traits_completed[,column_species]

  return(data_traits_completed)
}

#' Title: get_traits_higher_edna
#'
#' Find the parameter for the metaweb - for both species and genus-level here
#'
#' @description
#'
#'# End goal:
# 1. taxa
# 2. Common length
# 3. trophic_level
# 4. Env_2
#'
#' @param x Description of the first parameter.
#' @param y Description of the second parameter (if applicable).
#' @param ... Other optional parameters passed to methods.
#'
#' @details
#'
#' @return
#' Description of the object that the function returns.
#' If the function doesn't return anything meaningful, you can say `NULL`.
#'
#' @examples
#'
#'
#' @import rfishbase
#' @import dplyr
#' @import stringr
#'
#' @export
#'


get_traits_higher_edna <- function(data_presence, column_species = "Species", column_taxon = "taxon"){

  # Get all fishbase - first fetch all species names
  all_fishbase <- rfishbase::load_taxa()
  # Then get info on species table
  estimate_table <- rfishbase::estimate()
  species_table <- rfishbase::species() # Env_2 = DemersPelag
  info_fishbase <- species_table |>
    left_join(all_fishbase[,c("SpecCode", "Species")], by = c("SpecCode")) |>
    # Add estimate table
    left_join(estimate_table[,c("SpecCode", "Troph")]) |>
    # Add a common_length_estimated with 0.60 of length if NA
    mutate(CommonLengthEstim = ifelse(is.na(CommonLength), 0.6*Length, CommonLength)) |>
    mutate(Species = gsub(" ", "_", Species)) |>
    dplyr::select(Species, CommonLengthEstim, TrophicLevel = Troph, DemersPelag) |>
    mutate(Genus = stringr::word(Species, 1, 1, "_"))

  # Intersect with data_presence
  data_traits_species <- data.frame(species = data_presence[,column_species])
  data_traits_species <- data_traits_species |>
                            distinct(species) |>
                            filter(species != "")

  # Get the taxon but remove those that were species
  data_traits_taxon <- data.frame(taxon = data_presence[,column_taxon])
  data_traits_taxon <- data_traits_taxon |>
                            distinct(taxon) |>
                            filter(taxon != "")
  data_traits_taxon <- data_traits_taxon |>
                            filter(!(taxon %in% data_traits_species[,1]))

  # Species-level
  data_traits_completed_species <- data_traits_species |>
    left_join(info_fishbase, by = c("species" = "Species")) |>
    dplyr::select(-Genus)
  rownames(data_traits_completed_species) <- data_traits_completed_species[,"species"]
  colnames(data_traits_completed_species)[1] <- "taxon"

  # Taxon-level (genus)
  data_traits_completed_genus <- info_fishbase %>%
    semi_join(data_traits_taxon, by = c("Genus" = "taxon")) %>%
    group_by(Genus) %>%
    summarise(
      CommonLengthEstim = mean(CommonLengthEstim, na.rm = TRUE),
      TrophicLevel = mean(TrophicLevel, na.rm = TRUE),
      DemersPelag = names(which.max(table(DemersPelag)))
    ) %>%
    ungroup() |> as.data.frame()
  rownames(data_traits_completed_genus) <- data_traits_completed_genus$Genus
  colnames(data_traits_completed_genus)[1] <- "taxon"

  # Bind both
  data_traits_completed_species_genus <- rbind(data_traits_completed_species, data_traits_completed_genus)

  # Return
  return(data_traits_completed_species_genus)

  }











