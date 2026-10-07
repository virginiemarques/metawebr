#' Harmonise taxon names to FishBase accepted names
#'
#' Resolves taxon names from eDNA assignments (e.g. NCBI taxonomy), surveys or
#' other sources to the names accepted by FishBase, so that they match the
#' traits and diets taken from FishBase and the names of other datasets. Each
#' name goes through these steps, and the first that succeeds is kept:
#'
#' 1. **Cleaning.** Underscores become spaces; `cf.` and `aff.` are dropped;
#'    `Genus sp.`, `Genus spp.` and `Genus sp. XYZ` become the genus;
#'    subspecies and author names are dropped; `Genus a/Genus b` and
#'    `Genus a/b` become the genus. Names such as `uncultured fish`,
#'    `unidentified` or hybrids (`A x B`) are marked `"not_a_taxon"`.
#' 2. **User table** (`lookup`): your own corrections, which override
#'    everything else.
#' 3. **FishBase** (offline once cached, see [clear_fishbase_cache()]):
#'    accepted species, genus and family names, then the FishBase synonym
#'    table. Misapplied names are ignored. A synonym pointing to several
#'    species is marked `"ambiguous"` and left unresolved.
#' 4. **WoRMS** (`worms = TRUE`, needs the \pkg{worrms} package and an
#'    internet connection): names still unresolved are matched with the WoRMS
#'    TAXAMATCH service, which also corrects spelling mistakes. The WoRMS
#'    accepted name is then looked up in FishBase again. Names WoRMS knows but
#'    FishBase doesn't, such as mammals or invertebrates in eDNA data, are
#'    marked `"not_in_fishbase"`.
#'
#' Use [rename_taxa()] to replace the names of a data frame by the harmonised
#' names.
#'
#' @param taxa Character vector of taxon names (species as `Genus species` or
#'   `Genus_species`, genera or families). Duplicates are allowed.
#' @param lookup Optional data frame with columns `input` and `accepted`: your
#'   own name corrections, applied first. `input` is matched against the
#'   original and the cleaned name (with spaces or underscores); `accepted` is
#'   used as is (`NA` marks the taxon as `"not_a_taxon"`).
#' @param worms Match names unresolved by FishBase with WoRMS?
#' @param aphia_id Also get the WoRMS AphiaID of every resolved name (one
#'   WoRMS request per 50 names)?
#'
#' @return A data frame with one row per unique name of `taxa`, in order:
#'   * `input`: the name as given;
#'   * `query`: the cleaned name;
#'   * `accepted`: the FishBase accepted name, as `Genus_species` (or the
#'     genus or family), `NA` when unresolved;
#'   * `rank`: `"species"`, `"genus"` or `"family"`;
#'   * `status`: `"accepted"` (already the FishBase name), `"synonym"`
#'     (FishBase synonym), `"user"` (from `lookup`), `"worms_exact"` or
#'     `"worms_fuzzy"` (matched through WoRMS, `"worms_fuzzy"` when the
#'     spelling was corrected), `"ambiguous"`, `"not_in_fishbase"`,
#'     `"not_a_taxon"` or `"unresolved"`;
#'   * `source`: `"user"`, `"fishbase"` or `"worms"`;
#'   * `SpecCode`: FishBase species code (species only);
#'   * `worms_name`, `AphiaID`: WoRMS accepted name and identifier, when WoRMS
#'     was queried;
#'   * `note`: what the cleaning step changed, or the candidate species of an
#'     ambiguous synonym.
#'
#'   Check the rows with `status` `"worms_fuzzy"`, `"ambiguous"` and
#'   `"unresolved"` by hand, and put your decisions in `lookup`.
#'
#' @seealso [rename_taxa()]
#'
#' @examples
#' \dontrun{
#' res <- harmonise_taxa(c("Liza aurata", "Diplodus sp.", "Gadus morua",
#'                         "Pomatoschistus minutus/lozanoi", "uncultured fish"))
#' res[, c("input", "accepted", "status")]
#' }
#'
#' @export
harmonise_taxa <- function(taxa, lookup = NULL, worms = TRUE, aphia_id = FALSE) {
  if (!is.character(taxa)) stop("`taxa` must be a character vector.", call. = FALSE)
  input <- unique(taxa[!is.na(taxa)])
  if (length(input) == 0) stop("`taxa` has no name.", call. = FALSE)
  if ((worms || aphia_id) && !requireNamespace("worrms", quietly = TRUE)) {
    stop("Package 'worrms' is needed to query WoRMS: install.packages(\"worrms\"), or use `worms = FALSE`.",
         call. = FALSE)
  }

  cleaned <- clean_taxon_names(input)
  out <- data.frame(input = input, query = cleaned$query, accepted = NA_character_,
                    rank = NA_character_, status = ifelse(cleaned$not_a_taxon, "not_a_taxon", "unresolved"),
                    source = NA_character_, SpecCode = NA_integer_, worms_name = NA_character_,
                    AphiaID = NA_integer_, note = cleaned$note)
  todo <- function() out$status == "unresolved"

  # User table
  if (!is.null(lookup)) {
    if (!is.data.frame(lookup) || !all(c("input", "accepted") %in% names(lookup))) {
      stop("`lookup` must be a data frame with columns `input` and `accepted`.", call. = FALSE)
    }
    key <- gsub("_", " ", as.character(lookup$input))
    i <- match(gsub("_", " ", out$input), key)
    i[is.na(i)] <- match(out$query[is.na(i)], key)
    hit <- !is.na(i)
    acc <- gsub(" ", "_", as.character(lookup$accepted[i[hit]]))
    out$accepted[hit] <- acc
    out$status[hit] <- ifelse(is.na(acc), "not_a_taxon", "user")
    out$source[hit] <- "user"
    out$rank[hit] <- guess_rank(acc)
  }

  # FishBase
  sel <- which(todo())
  fb <- fishbase_match(out$query[sel])
  out[sel, c("accepted", "rank", "status", "SpecCode")] <- fb[, c("accepted", "rank", "status", "SpecCode")]
  out$note[sel] <- ifelse(is.na(fb$note), out$note[sel], fb$note)
  out$source[out$status %in% c("accepted", "synonym")] <- "fishbase"

  # WoRMS, then FishBase again on the WoRMS accepted name
  # (ambiguous FishBase synonyms too: the WoRMS accepted name can settle them)
  sel <- which(out$status %in% c("unresolved", "ambiguous"))
  if (worms && length(sel) > 0) {
    wm <- worms_match(out$query[sel])
    found <- !is.na(wm$worms_name)
    out$worms_name[sel] <- wm$worms_name
    out$AphiaID[sel] <- wm$AphiaID
    fb <- fishbase_match(wm$worms_name[found])
    rows <- sel[found]
    resolved <- fb$status %in% c("accepted", "synonym")
    out[rows[resolved], c("accepted", "rank", "SpecCode")] <- fb[resolved, c("accepted", "rank", "SpecCode")]
    out$status[rows[resolved]] <- ifelse(wm$match_type[found][resolved] == "exact", "worms_exact", "worms_fuzzy")
    out$source[rows[resolved]] <- "worms"
    # still unresolved: unknown to FishBase, unless FishBase found it ambiguous
    still <- rows[!resolved]
    out$status[still] <- ifelse(out$status[still] == "ambiguous" | fb$status[!resolved] == "ambiguous",
                                "ambiguous", "not_in_fishbase")
  }

  if (aphia_id) {
    sel <- which(!is.na(out$accepted) & is.na(out$AphiaID))
    if (length(sel) > 0) {
      wm <- worms_match(gsub("_", " ", out$accepted[sel]))
      out$worms_name[sel] <- wm$worms_name
      out$AphiaID[sel] <- wm$AphiaID
    }
  }

  report_harmonisation(out)
  out
}

#' Replace taxon names by harmonised names
#'
#' Replaces the names of a column by the accepted names found by
#' [harmonise_taxa()], keeping the original names in a new column. Several
#' original names can become the same accepted name (synonyms, or species
#' merged into a genus); functions such as [make_species_station_matrix()]
#' then sum their counts.
#'
#' @param data Data frame, e.g. eDNA detections.
#' @param column Name of the column of taxon names.
#' @param harmonised Output of [harmonise_taxa()] for these names.
#' @param drop_unresolved Remove rows whose name has no accepted name
#'   (`"unresolved"`, `"ambiguous"`, `"not_in_fishbase"` or `"not_a_taxon"`)?
#'   Otherwise they keep their original name, with underscores.
#'
#' @return `data` with `column` replaced and a column `<column>_original`.
#'
#' @examples
#' det <- data.frame(taxon = c("Liza aurata", "Chelon_auratus", "Gobius sp."), reads = c(10, 5, 3))
#' harm <- data.frame(input = c("Liza aurata", "Chelon_auratus", "Gobius sp."),
#'                    accepted = c("Chelon_auratus", "Chelon_auratus", "Gobius"))
#' rename_taxa(det, "taxon", harm)
#'
#' @export
rename_taxa <- function(data, column, harmonised, drop_unresolved = FALSE) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  if (!column %in% names(data)) stop(sprintf("`data` has no column `%s`.", column), call. = FALSE)
  if (!is.data.frame(harmonised) || !all(c("input", "accepted") %in% names(harmonised))) {
    stop("`harmonised` must be the output of harmonise_taxa().", call. = FALSE)
  }
  original <- as.character(data[[column]])
  accepted <- harmonised$accepted[match(original, harmonised$input)]
  missing <- unique(original[!is.na(original) & !original %in% harmonised$input])
  if (length(missing) > 0) {
    warning(sprintf("%d names are not in `harmonised` and were kept: %s.",
                    length(missing), format_taxa(missing)), call. = FALSE)
  }

  data[[paste0(column, "_original")]] <- original
  data[[column]] <- ifelse(is.na(accepted), gsub(" ", "_", original), accepted)
  if (drop_unresolved) {
    drop <- is.na(accepted)
    if (any(drop)) {
      message(sprintf("%d rows without an accepted name removed (%s).", sum(drop),
                      format_taxa(unique(original[drop]))))
      data <- data[!drop, , drop = FALSE]
    }
  }
  n_merged <- sum(tapply(data[[paste0(column, "_original")]], data[[column]],
                         function(x) length(unique(x))) > 1)
  if (n_merged > 0) message(sprintf("%d accepted names gather several original names.", n_merged))
  data
}

# Cleaned query names (with spaces), whether they are not taxa at all, and a
# note on what was changed
clean_taxon_names <- function(x) {
  q <- gsub("\\s+", " ", trimws(gsub("_", " ", x)))
  q <- sub("^([a-z])", "\\U\\1", q, perl = TRUE)  # capitalise the genus
  note <- rep(NA_character_, length(q))
  set_note <- function(sel, text) note[sel] <<- ifelse(is.na(note[sel]), text, paste(note[sel], text, sep = "; "))

  not_a_taxon <- q == "" |
    grepl("^(uncultured|unidentified|unclassified|unknown|environmental|metagenome)\\b", q,
          ignore.case = TRUE, perl = TRUE) |
    grepl(" x |\u00d7", q)

  # Qualifiers
  sel <- grepl("\\b(cf|aff)\\.? ", q, perl = TRUE)
  q[sel] <- gsub("\\b(cf|aff)\\.? ", "", q[sel], perl = TRUE)
  set_note(sel, "qualifier cf./aff. dropped")

  # Alternatives "Genus a/Genus b" or "Genus a/b": the genus when shared
  sel <- grepl("/", q) & !not_a_taxon
  for (i in which(sel)) {
    parts <- trimws(strsplit(q[i], "/")[[1]])
    genera <- unique(sub(" .*$", "", parts[grepl("^[A-Z]", parts)]))
    if (length(genera) == 1) {
      q[i] <- genera
      note[i] <- "several species: genus kept"
    } else {
      not_a_taxon[i] <- TRUE
      note[i] <- "several genera"
    }
  }

  # Genus sp. / spp. / sp. XYZ
  sel <- grepl("^[A-Z][a-z]+ spp?\\.?( .*)?$", q) & !not_a_taxon
  q[sel] <- sub(" .*$", "", q[sel])
  set_note(sel, "sp. dropped: genus")

  # Subspecies, varieties and author names
  sel <- grepl("^[A-Z][a-z]+ [a-z-]+ .+$", q) & !not_a_taxon
  q[sel] <- sub("^([A-Z][a-z]+ [a-z-]+) .+$", "\\1", q[sel])
  set_note(sel, "reduced to species")

  list(query = q, not_a_taxon = not_a_taxon, note = note)
}

# Match cleaned names (with spaces) to FishBase: accepted species, genera and
# families, then synonyms. Returns accepted (with underscores), rank, status
# ("accepted", "synonym", "ambiguous" or "unresolved"), SpecCode and note.
fishbase_match <- function(q) {
  n <- length(q)
  out <- data.frame(accepted = rep(NA_character_, n), rank = rep(NA_character_, n),
                    status = rep("unresolved", n), SpecCode = rep(NA_integer_, n), note = rep(NA_character_, n))
  if (n == 0) return(out)
  taxa <- fb_table("taxa")
  syn <- fb_table("synonyms")

  i <- match(q, taxa$Species)
  hit <- !is.na(i)
  out[hit, c("accepted", "rank", "status")] <- list(taxa$Species[i[hit]], "species", "accepted")
  out$SpecCode[hit] <- taxa$SpecCode[i[hit]]

  for (level in c("Genus", "Family")) {
    hit <- out$status == "unresolved" & q %in% taxa[[level]]
    out[hit, c("accepted", "rank", "status")] <- list(q[hit], tolower(level), "accepted")
  }

  todo <- which(out$status == "unresolved")
  syn <- syn[syn$synonym %in% q[todo] & syn$Status != "misapplied name" &
               syn$SpecCode %in% taxa$SpecCode, c("synonym", "Status", "SpecCode")]
  for (k in todo) {
    s <- syn[syn$synonym == q[k], , drop = FALSE]
    if (nrow(s) == 0) next
    codes <- unique(s$SpecCode[s$Status == "accepted name"])
    if (length(codes) == 0) codes <- unique(s$SpecCode[s$Status %in% c("synonym", "provisionally accepted name")])
    if (length(codes) == 0) codes <- unique(s$SpecCode)
    if (length(codes) == 1) {
      out$accepted[k] <- taxa$Species[match(codes, taxa$SpecCode)]
      out[k, c("rank", "status")] <- list("species", "synonym")
      out$SpecCode[k] <- codes
    } else {
      out$status[k] <- "ambiguous"
      out$note[k] <- paste("synonym of", paste(taxa$Species[match(codes, taxa$SpecCode)], collapse = ", "))
    }
  }
  out$accepted <- gsub(" ", "_", out$accepted)
  out
}

# WoRMS TAXAMATCH, 50 names per request. Returns worms_name, AphiaID (of the
# accepted name) and match_type, NA when not found.
worms_match <- function(names) {
  out <- data.frame(worms_name = rep(NA_character_, length(names)), AphiaID = NA_integer_,
                    match_type = NA_character_)
  chunks <- split(seq_along(names), ceiling(seq_along(names) / 50))
  for (idx in chunks) {
    res <- tryCatch(worrms::wm_records_taxamatch(names[idx], marine_only = FALSE),
                    error = function(e) {
                      warning("WoRMS request failed: ", conditionMessage(e), call. = FALSE)
                      NULL
                    })
    if (is.null(res)) next
    for (j in seq_along(idx)) {
      d <- as.data.frame(res[[j]])
      if (nrow(d) == 0 || !"valid_name" %in% names(d)) next
      d <- d[!is.na(d$valid_name), , drop = FALSE]
      if (nrow(d) == 0) next
      best <- if (any(d$status == "accepted")) which(d$status == "accepted")[1] else 1
      out$worms_name[idx[j]] <- d$valid_name[best]
      out$AphiaID[idx[j]] <- as.integer(d$valid_AphiaID[best])
      out$match_type[idx[j]] <- d$match_type[best]
    }
  }
  out
}

# Rank from the form of a name: Genus_species, a family (-idae) or a genus
guess_rank <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(grepl("_", x), "species", ifelse(grepl("idae$", x), "family", "genus")))
}

report_harmonisation <- function(out) {
  counts <- table(factor(out$status, levels = c("accepted", "synonym", "user", "worms_exact", "worms_fuzzy",
                                                "ambiguous", "not_in_fishbase", "not_a_taxon", "unresolved")))
  counts <- counts[counts > 0]
  message(sprintf("%d names: %s", nrow(out), paste(names(counts), counts, sep = " ", collapse = ", ")))
  check <- out$input[out$status %in% c("worms_fuzzy", "ambiguous", "unresolved")]
  if (length(check) > 0) {
    wrapped_message(sprintf("Check by hand (status worms_fuzzy, ambiguous or unresolved): %s",
                            format_taxa(check, n = 10)))
  }
}
