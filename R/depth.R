#' Add depth ranges to a trait table
#'
#' Adds the depth range of each taxon (`DepthMin`, `DepthMax`, in m), used by
#' the `depth_overlap` correction of [correct_metaweb_fish()], from FishBase
#' (`DepthRangeShallow`, `DepthRangeDeep`) or from [df_depth_env]. Genera and
#' families get the widest range of their species. A `Depth_source` column
#' records where each range comes from; it is reported for the links removed by
#' the depth correction.
#'
#' Existing `DepthMin` and `DepthMax` values are kept unless `overwrite = TRUE`,
#' so you can fill gaps from a second source:
#' `add_depth_range(add_depth_range(traits), "df_depth_env")`.
#'
#' @param data_traits Traits with taxon names as row names (`Genus_species`,
#'   genera or families).
#' @param source `"fishbase"` (default) or `"df_depth_env"`.
#' @param overwrite Replace existing depth values?
#'
#' @return `data_traits` with the columns `DepthMin`, `DepthMax` and
#'   `Depth_source`. Taxa without data keep `NA` and are not affected by the
#'   depth correction.
#'
#' @examples
#' traits <- data.frame(TrophicLevel = c(3.1, 4.1),
#'                      row.names = c("Acroteriobatus_annulatus", "Unknown_fish"))
#' add_depth_range(traits, source = "df_depth_env")
#'
#' @export
add_depth_range <- function(data_traits, source = c("fishbase", "df_depth_env"), overwrite = FALSE) {
  source <- match.arg(source)
  check_traits(data_traits)
  taxa <- rownames(data_traits)

  ref <- if (source == "fishbase") {
    info <- fishbase_traits_table()
    info[, c("Species", "Genus", "Family", "DepthMin", "DepthMax")]
  } else {
    d <- metawebr::df_depth_env
    sp <- gsub(" ", "_", d$Species)
    data.frame(Species = sp, Genus = sub("_.*$", "", sp), Family = NA_character_,
               DepthMin = d$Depth_min, DepthMax = d$Depth_max)
  }
  ranges <- t(vapply(taxa, function(x) {
    rows <- ref$Species == x
    if (!any(rows)) rows <- ref$Genus %in% x
    if (!any(rows)) rows <- ref$Family %in% x
    c(min_or_na(ref$DepthMin[rows]), max_or_na(ref$DepthMax[rows]))
  }, numeric(2)))

  for (col in c("DepthMin", "DepthMax", "Depth_source")) {
    if (!col %in% names(data_traits)) data_traits[[col]] <- if (col == "Depth_source") NA_character_ else NA_real_
  }
  fill <- (overwrite | is.na(data_traits$DepthMin) | is.na(data_traits$DepthMax)) &
    !is.na(ranges[, 1]) & !is.na(ranges[, 2])
  data_traits$DepthMin[fill] <- ranges[fill, 1]
  data_traits$DepthMax[fill] <- ranges[fill, 2]
  data_traits$Depth_source[fill] <- source

  bad <- which(data_traits$DepthMin > data_traits$DepthMax)
  if (length(bad) > 0) {
    warning(sprintf("%d taxa have DepthMin > DepthMax: %s.", length(bad), format_taxa(taxa[bad])), call. = FALSE)
  }
  n_missing <- sum(is.na(data_traits$DepthMin) | is.na(data_traits$DepthMax))
  message(sprintf("Depth ranges from %s for %d taxa; %d taxa without a depth range.",
                  source, sum(fill), n_missing))
  data_traits
}

min_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0) min(x) else NA_real_
}

max_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) > 0) max(x) else NA_real_
}
