# metawebr <img src="man/figures/logo.png" align="right" width="120" />

Build regional fish food webs (*metawebs*) from body size, FishBase traits and
diet data, and compute network indices for local communities such as eDNA
samples or surveys.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("virginiemarques/metawebr", build_vignettes = TRUE)
```

## What it does

**1. Build and curate a metaweb**
(`vignette("metaweb", package = "metawebr")`)

- calibrate the allometric niche model of Gravel et al. (2013) on
  predator–prey body sizes, and cross-validate it;
- harmonise taxon names (eDNA/NCBI, surveys) to FishBase accepted names, with
  WoRMS as a fallback;
- gather traits from FishBase, and fill gaps from congeners or with
  phylogenetic imputation;
- predict the probability of every trophic link, then correct it with
  ecological rules: herbivory (from trophic levels or FishBase diet data),
  water column, depth overlap, body size;
- record which rule and which data changed each link, and review the changes;
- binarise the metaweb and validate it against observed feeding links.

**2. Explore local food webs**
(`vignette("local-food-webs", package = "metawebr")`)

- extract the local food web of each site;
- compute about 45 network indices (connectance, trophic levels, omnivory,
  top/intermediate/basal taxa, modularity, path lengths...);
- compare them with null models;
- plot local webs, statically or interactively.

## Quick start

``` r
library(metawebr)

# Metaweb of the species in the example OBIS data
traits <- clean_traits(df_traits_obis)
pars <- c(a0 = -0.5616698, a1 = 0.9611534, b0 = 0.1247726, b1 = 0.1359744)
mw <- apply_model_metaweb(traits, pars)
mw <- correct_metaweb_fish(mw, traits)
mw_bin <- correct_metaweb_proba(mw, traits)

# Local food webs and their indices
pa <- make_species_station_matrix(df_presence_obis, station_col = "station",
                                  species_col = "species", count_col = "count")
pa <- pa[, colnames(pa) %in% rownames(mw_bin)]
indices <- get_indic_cells(pa, mw_bin)
```

FishBase tables are downloaded once per session; set
`options(metawebr.cache_dir = "<folder>")` to keep them between sessions.

## Reference

Gravel, D., Poisot, T., Albouy, C., Velez, L. & Mouillot, D. (2013). Inferring
food web structure from predator–prey body size relationships. *Methods in
Ecology and Evolution*, 4, 1083–1090.
