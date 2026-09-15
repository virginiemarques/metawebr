# Build `df_presence_obis` and `df_traits_obis`.
#
# Fish checklists from OBIS (records since 2000) for three Mediterranean
# rocky-reef sites and three tropical western Atlantic coral-reef sites, keeping
# the `n_top` species with the most records at each site, and their FishBase
# traits. Requires internet access.
#
# Run from the package root: source("data-raw/df_presence_obis.R")

devtools::load_all()

sites <- data.frame(
  station = c("Marseille", "Medes_Islands", "Port_Cros",
              "Florida_Keys", "Belize_Barrier_Reef", "Bonaire"),
  region = rep(c("Mediterranean", "Atlantic_reef"), each = 3),
  xmin = c(5.20, 3.15, 6.30, -81.60, -88.20, -68.45),
  xmax = c(5.50, 3.30, 6.50, -80.30, -87.70, -68.15),
  ymin = c(43.10, 41.95, 42.95, 24.40, 16.80, 12.00),
  ymax = c(43.30, 42.10, 43.05, 25.00, 17.50, 12.35)
)

fish_taxa <- c(Actinopterygii = 10194, Elasmobranchii = 10193)  # WoRMS AphiaIDs
n_top <- 40

bbox_wkt <- function(s) {
  sprintf("POLYGON((%s %s, %s %s, %s %s, %s %s, %s %s))",
          s$xmin, s$ymin, s$xmax, s$ymin, s$xmax, s$ymax, s$xmin, s$ymax, s$xmin, s$ymin)
}

site_checklist <- function(s) {
  cl <- do.call(rbind, lapply(fish_taxa, function(id) {
    x <- as.data.frame(robis::checklist(taxonid = id, geometry = bbox_wkt(s), startdate = "2000-01-01"))
    x[, c("species", "taxonRank", "is_marine", "records")]
  }))
  cl <- cl[cl$taxonRank == "Species" & !is.na(cl$species) & cl$is_marine %in% TRUE, ]
  cl <- cl[order(-cl$records, cl$species), ]
  cl <- head(cl, n_top)
  data.frame(station = s$station,
             region = s$region,
             lon = (s$xmin + s$xmax) / 2,
             lat = (s$ymin + s$ymax) / 2,
             species = gsub(" ", "_", cl$species),
             count = cl$records)
}

df_presence_obis <- do.call(rbind, lapply(split(sites, seq_len(nrow(sites))), site_checklist))
rownames(df_presence_obis) <- NULL

df_traits_obis <- get_traits(df_presence_obis, column_species = "species")

usethis::use_data(df_presence_obis, df_traits_obis, overwrite = TRUE, compress = "xz")
