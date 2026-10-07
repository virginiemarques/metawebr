# Example binary metaweb of the species of df_presence_obis, built as in
# vignette("metaweb"): niche model calibrated on df_interaction_fish, default
# corrections, threshold chosen from FishBase trophic levels.
devtools::load_all()

pars <- c(a0 = -0.5616698, a1 = 0.9611534, b0 = 0.1247726, b1 = 0.1359744)
traits <- clean_traits(df_traits_obis)
mw <- apply_model_metaweb(traits, pars)
mw_corr <- correct_metaweb_fish(mw, traits)
threshold <- get_proba(mw_corr, traits, plot_diag = FALSE)
metaweb_obis <- correct_metaweb_proba(mw_corr, threshold = threshold)
attr(metaweb_obis, "data_sources") <- NULL
attr(metaweb_obis, "threshold") <- threshold

usethis::use_data(metaweb_obis, overwrite = TRUE)
