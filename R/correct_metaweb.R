#' Correct a fish metaweb with ecological rules
#'
#' The body-size model alone predicts some implausible links. This function
#' removes or rewires them using the traits of each taxon. Corrections are
#' applied in this order:
#'
#' * `herbivory`: herbivores eat no other taxa. Herbivores are taxa with a
#'   trophic level below 2.4, unless `data_traits` has a `Herbivore` column
#'   (see [set_herbivores()]). The message lists the taxa treated as
#'   herbivores so you can check them.
#' * `water_pos`: removes links between taxa living in incompatible parts of
#'   the water column (column `DemersPelag`, e.g. bathydemersal vs.
#'   pelagic-oceanic), and removes cannibalism.
#' * `small_fish`: taxa shorter than 10 cm eat no other taxa.
#' * `PP_SP_add`: adds `PrimaryProducer` and `SecondaryProducer` nodes.
#'   Herbivores eat primary producers; taxa with a trophic level in (2.4, 3.7]
#'   and small fish eat secondary producers; secondary producers eat primary
#'   producers and themselves.
#' * `exception_BigFish`: taxa longer than 100 cm with a trophic level in
#'   (2.4, 3.7), such as large planktivorous sharks, eat only secondary
#'   producers. Requires the producer nodes.
#'
#' Taxa with a missing trait are left untouched by the corrections that use
#' it, with a warning.
#'
#' The correction that changed each link is recorded in the result; see
#' [compare_metawebs()] to list the changed links for manual checks.
#'
#' @param data_meta Metaweb from [apply_model_metaweb()] (prey in rows,
#'   predators in columns).
#' @param data_traits Traits with taxon names as row names, containing every
#'   taxon of `data_meta`. Needs `TrophicLevel`, `CommonLengthEstim` (cm) and
#'   `DemersPelag`, or only those used by the selected corrections. An optional
#'   logical `Herbivore` column overrides the trophic-level rule for herbivores.
#' @param herbivory,water_pos,small_fish,PP_SP_add,exception_BigFish Logical:
#'   apply this correction?
#'
#' @return The corrected matrix. With `PP_SP_add = TRUE` it has two extra rows
#'   and columns, `PrimaryProducer` and `SecondaryProducer`.
#'
#' @seealso [set_herbivores()] to change which taxa are treated as herbivores,
#'   [compare_metawebs()] to check the changed links.
#'
#' @examples
#' traits <- data.frame(
#'   CommonLengthEstim = c(100, 30, 5),
#'   TrophicLevel = c(4.1, 3.4, 3.3),
#'   DemersPelag = c("benthopelagic", "benthopelagic", "benthopelagic"),
#'   row.names = c("Gadus_morhua", "Clupea_harengus", "Gasterosteus_aculeatus")
#' )
#' mw <- apply_model_metaweb(traits, c(-0.56, 0.96, 0.12, 0.14))
#' mw_corr <- correct_metaweb_fish(mw, traits)
#' compare_metawebs(mw, mw_corr, traits)
#'
#' @export
correct_metaweb_fish <- function(data_meta,
                                 data_traits,
                                 herbivory = TRUE,
                                 water_pos = TRUE,
                                 small_fish = TRUE,
                                 PP_SP_add = TRUE,
                                 exception_BigFish = TRUE) {
  data_meta <- check_metaweb(data_meta, "data_meta")
  attr(data_meta, "corrections") <- NULL
  has_herb_col <- is.data.frame(data_traits) && "Herbivore" %in% names(data_traits)
  use_tl <- (herbivory && !has_herb_col) || PP_SP_add || exception_BigFish
  use_length <- small_fish || PP_SP_add || exception_BigFish
  check_traits(data_traits, c(if (use_tl) "TrophicLevel",
                              if (use_length) "CommonLengthEstim",
                              if (water_pos) "DemersPelag"))
  if (PP_SP_add && any(producer_names %in% rownames(data_meta))) {
    stop("`data_meta` already contains producer nodes; use `PP_SP_add = FALSE`.", call. = FALSE)
  }

  taxa <- setdiff(rownames(data_meta), producer_names)
  check_taxa_in_traits(taxa, data_traits)
  traits <- data_traits[taxa, , drop = FALSE]
  if (use_tl) warn_missing_trait(traits, "TrophicLevel", "trophic-level based")
  if (use_length) warn_missing_trait(traits, "CommonLengthEstim", "body-length based")
  if (water_pos) warn_missing_trait(traits, "DemersPelag", "water position")

  tl <- if ("TrophicLevel" %in% names(traits)) traits$TrophicLevel
  len <- if (use_length) traits$CommonLengthEstim
  sp_herb <- herbivore_taxa(traits, taxa)
  sp_small <- taxa[which(len < 10)]

  # Last correction that changed each link
  input <- data_meta
  changed_by <- matrix(NA_character_, nrow(data_meta), ncol(data_meta), dimnames = dimnames(data_meta))

  # ******** HERBIVORY *********** #
  if (herbivory) {
    rule <- if (has_herb_col) "`Herbivore` column, or trophic level < 2.4 where it is NA" else "trophic level < 2.4"
    wrapped_message(sprintf("%d taxa treated as herbivores (%s), with no fish prey: %s",
                            length(sp_herb), rule,
                            if (length(sp_herb) > 0) paste(sp_herb, collapse = ", ") else "none"))
    wrapped_message(paste("If some are wrong, change them with",
                          "`data_traits <- set_herbivores(data_traits, add = , remove = )`",
                          "and run correct_metaweb_fish() again."))
    previous <- data_meta
    data_meta[, sp_herb] <- 0
    changed_by[previous != data_meta] <- "herbivory"
    message("---------- Herbivory corrected")
  }

  # ******** WATER POSITION *********** #
  if (water_pos) {
    previous <- data_meta
    data_meta[taxa, taxa] <- Distri_correction(traits$DemersPelag, data_meta[taxa, taxa, drop = FALSE])
    changed_by[previous != data_meta] <- "water_position"

    previous <- data_meta
    idx <- match(taxa, rownames(data_meta))
    data_meta[cbind(idx, idx)] <- 0
    changed_by[previous != data_meta] <- "cannibalism"
    message("---------- Water position corrected")
  }

  # ******** SMALL FISH *********** #
  if (small_fish) {
    message(sprintf("There are %d small taxa (< 10 cm) in the dataset", length(sp_small)))
    previous <- data_meta
    data_meta[, sp_small] <- 0
    changed_by[previous != data_meta] <- "small_fish"
    message("---------- Small fish corrected")
  }

  # ******** PP and SP addition *********** #
  if (PP_SP_add) {
    sp_eps <- taxa[which(tl > 2.4 & tl <= 3.7)]
    n <- nrow(data_meta)
    node_names <- c(rownames(data_meta), producer_names)
    out <- matrix(0, n + 2, n + 2, dimnames = list(node_names, node_names))
    out[seq_len(n), seq_len(n)] <- data_meta
    previous <- out
    out["PrimaryProducer", sp_herb] <- 1
    out["SecondaryProducer", unique(c(sp_eps, sp_small))] <- 1
    out[producer_names, "SecondaryProducer"] <- 1

    changed_out <- matrix(NA_character_, n + 2, n + 2, dimnames = dimnames(out))
    changed_out[seq_len(n), seq_len(n)] <- changed_by
    changed_out[previous != out] <- "producers"
    data_meta <- out
    changed_by <- changed_out
    message("---------- PP and SP corrected")
  }

  # ******** exception_BigFish *********** #
  if (exception_BigFish) {
    if (!"SecondaryProducer" %in% rownames(data_meta)) {
      warning("`exception_BigFish` needs the producer nodes (`PP_SP_add = TRUE`); correction skipped.",
              call. = FALSE)
    } else {
      big_basal_sp <- taxa[which(tl > 2.4 & tl < 3.7 & len > 100)]
      previous <- data_meta
      data_meta[, big_basal_sp] <- 0
      data_meta["SecondaryProducer", big_basal_sp] <- 1
      changed_by[previous != data_meta] <- "big_fish"
      message("---------- Exception - big fish corrected")
    }
  }

  attr(data_meta, "corrections") <- correction_log(input, data_meta, changed_by)
  data_meta
}

# Links changed by correct_metaweb_fish(): prey, predator, value before and
# after, and the last correction that changed them. Input taxa come first in
# the output, in the same order.
correction_log <- function(input, output, changed_by) {
  before <- matrix(0, nrow(output), ncol(output))
  before[seq_len(nrow(input)), seq_len(ncol(input))] <- input
  idx <- which(!is.na(changed_by) & before != output, arr.ind = TRUE)
  data.frame(prey = rownames(output)[idx[, 1]],
             predator = colnames(output)[idx[, 2]],
             before = before[idx],
             after = unclass(output)[idx],
             correction = changed_by[idx])
}

wrapped_message <- function(text, width = 80) {
  message(paste(strwrap(text, width = width, exdent = 2), collapse = "\n"))
}

# Remove links between taxa of incompatible water-column habitats.
# `data`: habitat of each taxon, in the order of the rows of `Lniche`.
# Taxa with a missing or unknown habitat keep all their links.
Distri_correction <- function(data, Lniche) {
  incompatible <- habitat_incompatibility()
  code <- match(as.character(data), rownames(incompatible))
  known <- !is.na(code)
  remove <- matrix(FALSE, length(code), length(code))
  remove[known, known] <- incompatible[code[known], code[known]]
  Lniche[remove] <- 0
  Lniche
}

# Symmetric logical matrix: TRUE when two FishBase habitats (DemersPelag) can't
# interact.
habitat_incompatibility <- function() {
  pairs <- rbind(
    c("bathydemersal", "pelagic-oceanic"),
    c("bathydemersal", "pelagic-neritic"),
    c("bathydemersal", "reef-associated"),
    c("bathydemersal", "pelagic"),
    c("bathypelagic", "pelagic-oceanic"),
    c("bathypelagic", "pelagic-neritic"),
    c("bathypelagic", "reef-associated"),
    c("pelagic-oceanic", "reef-associated"),
    c("pelagic-oceanic", "demersal"),
    c("pelagic-oceanic", "benthopelagic")
  )
  habitats <- sort(unique(c(pairs)))
  m <- matrix(FALSE, length(habitats), length(habitats), dimnames = list(habitats, habitats))
  m[pairs] <- TRUE
  m[pairs[, 2:1]] <- TRUE
  m
}
