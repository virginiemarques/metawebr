# metawebr <img src="man/figures/logo.png" align="right" width="120" />

Tools to build a regional fish food web (a *metaweb*) from body size and
FishBase traits, and to compute network indices for local communities, such as
eDNA samples, optionally compared with null models.

## Install

``` r
# install.packages("remotes")
remotes::install_git("git@gitlab.ethz.ch:vmarques/metaweb.git", build_vignettes = TRUE)

library(metawebr)
```

## Quick start

``` r
data(df_presence_obis)  # fish at 3 Mediterranean and 3 Atlantic reef sites (OBIS)
data(df_traits_obis)    # their FishBase traits

# Niche-model parameters, estimated with metaweb_mod_parameters()
pars <- c(a0 = -0.5616698, a1 = 0.9611534, b0 = 0.1247726, b1 = 0.1359744)

traits <- clean_traits(df_traits_obis)
metaweb <- apply_model_metaweb(traits, pars) |>
  correct_metaweb_fish(traits) |>
  correct_metaweb_proba(traits)

pa <- make_species_station_matrix(df_presence_obis, station_col = "station",
                                  species_col = "species", count_col = "count")
med <- pa[c("Marseille", "Medes_Islands", "Port_Cros"), ]

# Local indices compared with an equiprobable null model
res <- get_indic_cells(med, metaweb, null_model = "equiprobable", n_null = 199)
res$indices
plot_indic_null(res$null_model, indices = c("Connectance", "TL_moy", "Modularity"))
```

The full workflow is described in the vignette:
`vignette("tutorial", package = "metawebr")`.

## Tests

``` r
devtools::test()
```
