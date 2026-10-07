#' metawebr: build fish metawebs and compute local food-web indices
#'
#' Infers a regional fish food web (metaweb) from body size with the
#' allometric niche model of Gravel et al. (2013), curates it with ecological
#' rules and FishBase diet data, validates it against observed feeding links,
#' and computes network indices for local communities, optionally compared with
#' null models.
#'
#' Two vignettes walk through the workflow:
#' `vignette("metaweb", package = "metawebr")` and
#' `vignette("local-food-webs", package = "metawebr")`.
#'
#' @section Building the metaweb:
#' * Size model: [metaweb_mod_parameters()], [get_fig_eval()],
#'   [apply_model_metaweb()].
#' * Taxon names: [harmonise_taxa()], [rename_taxa()].
#' * Traits: [get_traits()], [get_traits_higher_edna()], [infer_traits()],
#'   [clean_traits()], [impute_traits_phylo()], [add_depth_range()].
#' * Curation: [set_herbivores()], [get_herbivores_fishbase()],
#'   [correct_metaweb_fish()], [compare_metawebs()].
#' * Binarisation: [get_proba()], [correct_metaweb_proba()].
#' * Validation: [cv_size_model()], [validate_metaweb()], [get_diet_links()].
#' * FishBase cache: [clear_fishbase_cache()], [fishbase_cache_info()].
#'
#' @section Local food webs:
#' * Site data: [make_species_station_matrix()].
#' * Indices: [get_indic_cells()].
#' * Null models: [get_indic_null()], [plot_indic_null()].
#' * Plots: [plot_tree_network()], [plot_network_focus_sp()],
#'   [plot_network_interactive()].
#'
#' @section Data:
#' [df_interaction_fish], [df_depth_env], [df_presence_obis], [df_traits_obis],
#' [metaweb_obis].
#'
#' @keywords internal
#' @importFrom stats dnorm lm median runif sd setNames
#' @importFrom graphics abline hist points
#' @importFrom utils head read.csv read.table write.table
"_PACKAGE"

# Columns used through non-standard evaluation in ggplot2 / ggraph
utils::globalVariables(c("name", "Trophic_Level", "Degree_in", "ses", "site", "significant"))
