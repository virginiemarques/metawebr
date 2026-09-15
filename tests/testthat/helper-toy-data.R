# Small datasets shared by the tests

toy_taxa <- c("Gadus_morhua", "Clupea_harengus", "Mallotus_villosus", "Reinhardtius_hippoglossoides",
              "Boreogadus_saida", "Somniosus_microcephalus", "Gasterosteus_aculeatus", "Lycodes_esmarkii")

toy_traits <- function() {
  data.frame(
    species = toy_taxa,
    CommonLengthEstim = c(100, 30, 15, 56, 25, 256.2, 5.1, 47.5),
    TrophicLevel = c(4.09, 3.38, 3.15, 4.38, 3.12, 4.22, 3.31, 3.91),
    DemersPelag = c("benthopelagic", "benthopelagic", "pelagic-oceanic", "benthopelagic",
                    "demersal", "benthopelagic", "benthopelagic", "bathydemersal"),
    row.names = toy_taxa
  )
}

toy_pars <- c(a0 = -0.5616698, a1 = 0.9611534, b0 = 0.1247726, b1 = 0.1359744)

# Corrected probability metaweb and its binary version (threshold 0.5)
toy_metaweb <- function(traits = toy_traits()) {
  suppressMessages(correct_metaweb_fish(apply_model_metaweb(traits, toy_pars), traits))
}

toy_sites <- function() {
  det <- data.frame(
    station = rep(c("s1", "s2", "s3"), c(5, 4, 3)),
    taxon = c(toy_taxa[1:5], toy_taxa[c(1, 3, 6, 8)], toy_taxa[c(2, 4, 7)]),
    count = 1
  )
  make_species_station_matrix(det, station_col = "station", species_col = "taxon", count_col = "count")
}
