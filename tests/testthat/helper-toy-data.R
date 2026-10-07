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

# Small FishBase tables, as returned by fb_table(), for offline tests
fake_fb_tables <- function() {
  taxa <- data.frame(
    SpecCode = 1:6,
    Species = c("Gadus morhua", "Chelon auratus", "Sarpa salpa", "Diplodus sargus",
                "Diplodus vulgaris", "Clupea harengus"),
    Genus = c("Gadus", "Chelon", "Sarpa", "Diplodus", "Diplodus", "Clupea"),
    Family = c("Gadidae", "Mugilidae", "Sparidae", "Sparidae", "Sparidae", "Clupeidae")
  )
  list(
    taxa = taxa,
    species = data.frame(SpecCode = 1:6, CommonLength = c(100, 40, 30, 25, 22, 30), Length = NA,
                         DemersPelag = "benthopelagic",
                         DepthRangeShallow = c(150, 0, 0, 0, 0, 0),
                         DepthRangeDeep = c(600, 20, 20, 50, 160, 300)),
    estimate = data.frame(SpecCode = 1:6, Troph = c(4.4, 2.4, 2.0, 3.1, 3.2, 3.2)),
    synonyms = data.frame(
      synonym = c("Liza aurata", "Diplodus sargus", "Sargus vulgaris", "Sargus vulgaris", "Gadus callarias"),
      Status = c("synonym", "misapplied name", "synonym", "synonym", "synonym"),
      SpecCode = c(2L, 5L, 4L, 5L, 1L)
    ),
    fooditems = data.frame(
      SpecCode = c(rep(2L, 4), rep(3L, 3), 1L, 1L, 1L, 6L),
      FoodI = c("detritus", "plants", "plants", "zoobenthos", "plants", "plants", "zoobenthos",
                "nekton", "nekton", "zoobenthos", "zooplankton"),
      PreySpecCode = c(rep(NA, 7), 6L, 4L, NA, NA)
    ),
    ecology = data.frame(SpecCode = c(1L, 4L, 5L),
                         FeedingType = c("hunting macrofauna (predator)", "grazing on aquatic plants",
                                         "variable"))
  )
}

local_fake_fishbase <- function(env = parent.frame()) {
  tables <- fake_fb_tables()
  local_mocked_bindings(fb_table = function(table) tables[[table]], .env = env)
}
