#' Report and fill missing trait values
#'
#' Reports the number of missing values of each trait. Then, optionally, fills
#' missing values from other taxa of the same genus in the table (median length
#' and trophic level, most common water position) and removes the taxa that
#' still have missing traits. The genus is the part of the row name before the
#' first underscore, so a genus-level taxon (`Gadus`) is grouped with its
#' species (`Gadus_morhua`).
#'
#' @param data_traits Traits with taxon names as row names and the columns
#'   `CommonLengthEstim`, `TrophicLevel` and `DemersPelag`, e.g. from
#'   [get_traits()].
#' @param fill_by_genus Fill missing values from the same genus?
#' @param remove_incomplete Remove taxa with a missing trait (after filling)?
#'
#' @return `data_traits`, filled and/or filtered.
#'
#' @seealso [infer_traits()] to use all FishBase species of a genus or family.
#'
#' @examples
#' traits <- data.frame(
#'   CommonLengthEstim = c(100, NA, 30),
#'   TrophicLevel = c(4.1, 3.9, 3.4),
#'   DemersPelag = c("benthopelagic", NA, "benthopelagic"),
#'   row.names = c("Gadus_morhua", "Gadus_ogac", "Clupea_harengus")
#' )
#' clean_traits(traits)
#'
#' @export
clean_traits <- function(data_traits, fill_by_genus = TRUE, remove_incomplete = TRUE) {
  trait_cols <- c("CommonLengthEstim", "TrophicLevel", "DemersPelag")
  check_traits(data_traits, trait_cols)
  if (is.factor(data_traits$DemersPelag)) {
    data_traits$DemersPelag <- as.character(data_traits$DemersPelag)
  }

  for (col in trait_cols) {
    message(sprintf("There are %d missing values for %s out of %d taxa",
                    sum(is.na(data_traits[[col]])), col, nrow(data_traits)))
  }

  if (fill_by_genus) {
    genus <- sub("_.*$", "", rownames(data_traits))
    for (col in trait_cols) {
      na <- is.na(data_traits[[col]])
      if (!any(na)) next
      numeric_trait <- is.numeric(data_traits[[col]])
      by_genus <- vapply(split(data_traits[[col]], genus),
                         if (numeric_trait) median_or_na else mode_or_na,
                         if (numeric_trait) numeric(1) else character(1))
      filled <- by_genus[genus[na]]
      data_traits[[col]][na] <- filled
      if (any(!is.na(filled))) {
        message(sprintf("---------- %d missing values of %s filled from the same genus",
                        sum(!is.na(filled)), col))
      }
    }
  }

  if (remove_incomplete) {
    incomplete <- !stats::complete.cases(data_traits[, trait_cols])
    if (any(incomplete)) {
      message(sprintf("---------- %d taxa with missing traits removed: %s",
                      sum(incomplete), format_taxa(rownames(data_traits)[incomplete])))
      data_traits <- data_traits[!incomplete, , drop = FALSE]
    }
  }

  data_traits
}

#' Infer traits of species, genera and families from FishBase
#'
#' Gets traits for taxa identified at different taxonomic levels, e.g. eDNA
#' assignments. Species get their own FishBase traits. Genera and families get
#' the mean length and trophic level, and the most common water position, of
#' all FishBase species they contain (`method = "mean_level"`).
#'
#' With `method = "100matrix"`, trait values still missing afterwards (e.g.
#' species absent from FishBase) are then imputed from the phylogeny with
#' [impute_traits_phylo()], by default over the 100 complete fish trees of
#' fishtree.
#'
#' @param method `"mean_level"` (default) or `"100matrix"` to also impute
#'   missing values with [impute_traits_phylo()].
#' @param ... Arguments passed to [impute_traits_phylo()] when
#'   `method = "100matrix"`, e.g. `trees`, `n_runs` or `mc.cores`.
#' @param taxa Character vector of taxon names: species as `Genus_species`
#'   (spaces are converted to underscores), genera or families.
#' @param rank Taxonomic rank of each taxon (`"species"`, `"genus"` or
#'   `"family"`), recycled if of length 1. `NULL` (default) detects it by
#'   looking the names up in FishBase.
#' @param fishbase_table Optional table of FishBase species traits, to avoid
#'   downloading it again. Needs the columns `Species`, `Genus`, `Family`,
#'   `CommonLengthEstim`, `TrophicLevel` and `DemersPelag`.
#'
#' @return A data frame with taxon names as row names and the columns `taxon`,
#'   `rank`, `n_species` (number of FishBase species averaged),
#'   `CommonLengthEstim`, `TrophicLevel` and `DemersPelag`. Taxa not found in
#'   FishBase have missing traits, with a warning.
#'
#' @examples
#' \dontrun{
#' infer_traits(c("Gadus_morhua", "Sebastes", "Liparidae"))
#' }
#'
#' @export
infer_traits <- function(taxa, rank = NULL, fishbase_table = NULL,
                         method = c("mean_level", "100matrix"), ...) {
  method <- match.arg(method)
  out <- mean_level_traits(taxa, rank, fishbase_table)
  if (method == "100matrix") {
    out <- impute_traits_phylo(out, ...)
  }
  out
}

# FishBase traits of species, and means of genera and families
mean_level_traits <- function(taxa, rank = NULL, fishbase_table = NULL) {
  if (!is.character(taxa)) {
    stop("`taxa` must be a character vector.", call. = FALSE)
  }
  levels <- c(species = "Species", genus = "Genus", family = "Family")
  info <- if (is.null(fishbase_table)) fishbase_traits_table() else fishbase_table
  missing_cols <- setdiff(c(levels, "CommonLengthEstim", "TrophicLevel", "DemersPelag"), names(info))
  if (length(missing_cols) > 0) {
    stop(sprintf("`fishbase_table` is missing column(s): %s.", paste(missing_cols, collapse = ", ")),
         call. = FALSE)
  }

  taxa <- gsub(" ", "_", taxa)
  if (is.null(rank)) {
    rank <- ifelse(taxa %in% info$Species, "species",
                   ifelse(taxa %in% info$Genus, "genus",
                          ifelse(taxa %in% info$Family, "family", NA_character_)))
  } else {
    if (!length(rank) %in% c(1, length(taxa)) || !all(rank %in% names(levels))) {
      stop("`rank` must be \"species\", \"genus\" or \"family\", of length 1 or length(taxa).",
           call. = FALSE)
    }
    rank <- rep_len(rank, length(taxa))
  }
  keep <- !is.na(taxa) & taxa != "" & !duplicated(taxa)
  taxa <- taxa[keep]
  rank <- rank[keep]

  out <- data.frame(taxon = taxa, rank = rank, n_species = 0L, CommonLengthEstim = NA_real_,
                    TrophicLevel = NA_real_, DemersPelag = NA_character_, row.names = taxa)
  for (r in names(levels)) {
    sel <- which(rank == r)
    if (length(sel) == 0) next
    agg <- aggregate_traits(info, taxa[sel], levels[[r]])
    out[sel, names(agg)[-1]] <- agg[, -1]
  }
  warn_not_found(taxa[out$n_species == 0])
  out
}
