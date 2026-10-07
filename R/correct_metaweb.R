#' Correct a fish metaweb with ecological rules
#'
#' The body-size model alone predicts some implausible links. This function
#' removes or rewires them using the traits of each taxon. Corrections are
#' applied in this order:
#'
#' * `herbivory`: herbivores eat no other taxa. Herbivores are taxa with a
#'   trophic level below 2.4, unless `data_traits` has a `Herbivore` column
#'   (see [set_herbivores()], which can also use FishBase diet data). The
#'   message lists the taxa treated as herbivores so you can check them.
#' * `water_pos`: removes links between taxa living in incompatible parts of
#'   the water column (column `DemersPelag`, e.g. bathydemersal vs.
#'   pelagic-oceanic), and removes cannibalism.
#' * `depth_overlap`: removes links between taxa whose depth ranges
#'   (`DepthMin`, `DepthMax`, in m, see [add_depth_range()]) don't overlap,
#'   allowing a gap of `depth_tolerance` m. Off by default.
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
#' Each link changed is recorded in the result with the correction that
#' changed it and the data it was based on (`source`): for herbivory, the
#' `Herbivore_source` of the predator (`"trophic_level"`, `"diet_items"`,
#' `"feeding_type"`, `"manual"`); otherwise the trait columns used, such as
#' `"DemersPelag"`, or `"rule"` for fixed rules. See [compare_metawebs()] to
#' list the changed links for manual checks.
#'
#' @param data_meta Metaweb from [apply_model_metaweb()] (prey in rows,
#'   predators in columns).
#' @param data_traits Traits with taxon names as row names, containing every
#'   taxon of `data_meta`. Needs `TrophicLevel`, `CommonLengthEstim` (cm) and
#'   `DemersPelag`, or only those used by the selected corrections, and
#'   `DepthMin` and `DepthMax` for `depth_overlap`. An optional logical
#'   `Herbivore` column overrides the trophic-level rule for herbivores.
#' @param herbivory,water_pos,depth_overlap,small_fish,PP_SP_add,exception_BigFish
#'   Logical: apply this correction?
#' @param depth_tolerance Depth ranges closer than this (m) count as
#'   overlapping, for `depth_overlap`.
#'
#' @return The corrected matrix. With `PP_SP_add = TRUE` it has two extra rows
#'   and columns, `PrimaryProducer` and `SecondaryProducer`.
#'
#' @seealso [set_herbivores()] to change which taxa are treated as herbivores,
#'   [add_depth_range()], [compare_metawebs()] to check the changed links.
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
                                 exception_BigFish = TRUE,
                                 depth_overlap = FALSE,
                                 depth_tolerance = 0) {
  data_meta <- check_metaweb(data_meta, "data_meta")
  attr(data_meta, "corrections") <- NULL
  has_herb_col <- is.data.frame(data_traits) && "Herbivore" %in% names(data_traits)
  use_tl <- (herbivory && !has_herb_col) || PP_SP_add || exception_BigFish
  use_length <- small_fish || PP_SP_add || exception_BigFish
  check_traits(data_traits, c(if (use_tl) "TrophicLevel",
                              if (use_length) "CommonLengthEstim",
                              if (water_pos) "DemersPelag",
                              if (depth_overlap) c("DepthMin", "DepthMax")))
  if (!is.numeric(depth_tolerance) || length(depth_tolerance) != 1 || is.na(depth_tolerance) ||
      depth_tolerance < 0) {
    stop("`depth_tolerance` must be a single number >= 0.", call. = FALSE)
  }
  if (PP_SP_add && any(producer_names %in% rownames(data_meta))) {
    stop("`data_meta` already contains producer nodes; use `PP_SP_add = FALSE`.", call. = FALSE)
  }

  taxa <- setdiff(rownames(data_meta), producer_names)
  check_taxa_in_traits(taxa, data_traits)
  traits <- data_traits[taxa, , drop = FALSE]
  if (use_tl) warn_missing_trait(traits, "TrophicLevel", "trophic-level based")
  if (use_length) warn_missing_trait(traits, "CommonLengthEstim", "body-length based")
  if (water_pos) warn_missing_trait(traits, "DemersPelag", "water position")
  if (depth_overlap) {
    warn_missing_trait(data.frame(DepthMin_DepthMax = ifelse(is.na(traits$DepthMax), NA, traits$DepthMin),
                                  row.names = taxa),
                       "DepthMin_DepthMax", "depth overlap")
  }

  tl <- if ("TrophicLevel" %in% names(traits)) traits$TrophicLevel
  len <- if (use_length) traits$CommonLengthEstim
  sp_herb <- herbivore_taxa(traits, taxa)
  herb_source <- herbivore_sources(traits, taxa)
  sp_small <- taxa[which(len < 10)]

  # Last correction that changed each link, and the data it was based on
  input <- data_meta
  changes <- list(by = matrix(NA_character_, nrow(data_meta), ncol(data_meta), dimnames = dimnames(data_meta)))
  changes$source <- changes$by

  # ******** HERBIVORY *********** #
  if (herbivory) {
    rule <- if (has_herb_col) "`Herbivore` column, or trophic level < 2.4 where it is NA" else "trophic level < 2.4"
    by_source <- table(herb_source[sp_herb])
    wrapped_message(sprintf("%d taxa treated as herbivores (%s%s), with no fish prey: %s",
                            length(sp_herb), rule,
                            if (length(by_source) > 0) paste0("; sources: ", paste(names(by_source), by_source,
                                                                                  sep = " ", collapse = ", ")) else "",
                            if (length(sp_herb) > 0) paste(sp_herb, collapse = ", ") else "none"))
    wrapped_message(paste("If some are wrong, change them with",
                          "`data_traits <- set_herbivores(data_traits, add = , remove = )`",
                          "and run correct_metaweb_fish() again."))
    previous <- data_meta
    data_meta[, sp_herb] <- 0
    changes <- record_changes(changes, previous, data_meta, "herbivory", herb_source)
    message("---------- Herbivory corrected")
  }

  # ******** WATER POSITION *********** #
  if (water_pos) {
    previous <- data_meta
    data_meta[taxa, taxa] <- Distri_correction(traits$DemersPelag, data_meta[taxa, taxa, drop = FALSE])
    changes <- record_changes(changes, previous, data_meta, "water_position", "DemersPelag")

    previous <- data_meta
    idx <- match(taxa, rownames(data_meta))
    data_meta[cbind(idx, idx)] <- 0
    changes <- record_changes(changes, previous, data_meta, "cannibalism", "rule")
    message("---------- Water position corrected")
  }

  # ******** DEPTH OVERLAP *********** #
  if (depth_overlap) {
    previous <- data_meta
    apart <- depth_apart(traits$DepthMin, traits$DepthMax, depth_tolerance)
    sub <- data_meta[taxa, taxa, drop = FALSE]
    sub[apart] <- 0
    data_meta[taxa, taxa] <- sub
    depth_source <- if ("Depth_source" %in% names(traits)) {
      setNames(ifelse(is.na(traits$Depth_source), "DepthMin/DepthMax",
                      paste0("DepthMin/DepthMax (", traits$Depth_source, ")")), taxa)
    } else {
      "DepthMin/DepthMax"
    }
    changes <- record_changes(changes, previous, data_meta, "depth", depth_source)
    message("---------- Depth overlap corrected")
  }

  # ******** SMALL FISH *********** #
  if (small_fish) {
    message(sprintf("There are %d small taxa (< 10 cm) in the dataset", length(sp_small)))
    previous <- data_meta
    data_meta[, sp_small] <- 0
    changes <- record_changes(changes, previous, data_meta, "small_fish", "CommonLengthEstim")
    message("---------- Small fish corrected")
  }

  # ******** PP and SP addition *********** #
  if (PP_SP_add) {
    sp_eps <- taxa[which(tl > 2.4 & tl <= 3.7)]
    n <- nrow(data_meta)
    node_names <- c(rownames(data_meta), producer_names)
    out <- matrix(0, n + 2, n + 2, dimnames = list(node_names, node_names))
    out[seq_len(n), seq_len(n)] <- data_meta
    changes <- lapply(changes, function(m) {
      grown <- matrix(NA_character_, n + 2, n + 2, dimnames = dimnames(out))
      grown[seq_len(n), seq_len(n)] <- m
      grown
    })

    previous <- out
    out["PrimaryProducer", sp_herb] <- 1
    changes <- record_changes(changes, previous, out, "producers", herb_source)
    previous <- out
    out["SecondaryProducer", unique(c(sp_eps, sp_small))] <- 1
    changes <- record_changes(changes, previous, out, "producers", "TrophicLevel/CommonLengthEstim")
    previous <- out
    out[producer_names, "SecondaryProducer"] <- 1
    changes <- record_changes(changes, previous, out, "producers", "rule")
    data_meta <- out
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
      changes <- record_changes(changes, previous, data_meta, "big_fish", "TrophicLevel/CommonLengthEstim")
      message("---------- Exception - big fish corrected")
    }
  }

  log <- correction_log(input, data_meta, changes)
  attr(data_meta, "corrections") <- log
  # Data the corrections relied on, kept after binarisation (see validate_metaweb())
  attr(data_meta, "data_sources") <- unique(c(if (herbivory) unname(herb_source[sp_herb]), log$source))
  data_meta
}

# Record, for the links that differ between `previous` and `current`, the
# correction and its source. `source` is one string, or a vector named by
# predator (taxa without a name get "rule").
record_changes <- function(changes, previous, current, correction, source) {
  idx <- which(previous != current, arr.ind = TRUE)
  if (nrow(idx) == 0) return(changes)
  changes$by[idx] <- correction
  src <- if (length(source) == 1 && is.null(names(source))) {
    rep(source, nrow(idx))
  } else {
    unname(source[colnames(current)[idx[, 2]]])
  }
  src[is.na(src)] <- "rule"
  changes$source[idx] <- src
  changes
}

# Links changed by correct_metaweb_fish(): prey, predator, value before and
# after, and the last correction that changed them, with its source. Input
# taxa come first in the output, in the same order.
correction_log <- function(input, output, changes) {
  before <- matrix(0, nrow(output), ncol(output))
  before[seq_len(nrow(input)), seq_len(ncol(input))] <- input
  idx <- which(!is.na(changes$by) & before != output, arr.ind = TRUE)
  data.frame(prey = rownames(output)[idx[, 1]],
             predator = colnames(output)[idx[, 2]],
             before = before[idx],
             after = unclass(output)[idx],
             correction = changes$by[idx],
             source = changes$source[idx])
}

# Logical matrix (prey in rows, predators in columns): TRUE when the depth
# ranges of two taxa are more than `tolerance` apart. Taxa with a missing
# depth are never apart.
depth_apart <- function(dmin, dmax, tolerance = 0) {
  gap <- pmax(outer(dmin, dmax, "-"), outer(dmax, dmin, "-") * -1)  # > 0 when the ranges don't overlap
  out <- gap > tolerance
  out[is.na(out)] <- FALSE
  out
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
