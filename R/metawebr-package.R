#' metawebr: build fish metawebs and compute local food-web indices
#'
#' Infers a regional fish food web (metaweb) from body size with the
#' allometric niche model of Gravel et al. (2013), corrects it with simple
#' ecological rules, and computes network indices for local communities,
#' optionally compared with null models.
#'
#' @keywords internal
#' @importFrom stats dnorm lm median runif sd setNames
#' @importFrom graphics abline hist points
#' @importFrom utils head read.csv read.table write.table
"_PACKAGE"

# Columns used through non-standard evaluation in ggplot2 / ggraph
utils::globalVariables(c("name", "Trophic_Level", "Degree_in", "ses", "site", "significant"))
