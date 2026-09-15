# Internal input checks and small helpers --------------------------------------

producer_names <- c("PrimaryProducer", "SecondaryProducer")

# TRUE when row names hold taxon names, i.e. they are neither automatic
# (1, 2, 3, ...) nor purely numeric. Tibbles never have row names.
has_taxon_rownames <- function(x) {
  .row_names_info(x) > 0 && !all(grepl("^[0-9]+$", rownames(x)))
}

check_traits <- function(data_traits, cols = character(), arg = "data_traits") {
  if (!is.data.frame(data_traits)) {
    stop(sprintf("`%s` must be a data frame.", arg), call. = FALSE)
  }
  if (nrow(data_traits) == 0) {
    stop(sprintf("`%s` has no rows.", arg), call. = FALSE)
  }
  missing_cols <- setdiff(cols, names(data_traits))
  if (length(missing_cols) > 0) {
    stop(sprintf("`%s` is missing required column(s): %s.", arg,
                 paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  if (!has_taxon_rownames(data_traits)) {
    stop(sprintf(paste0("`%s` must have taxon names as row names, ",
                        "e.g. `rownames(%s) <- %s$species` (tibbles have no row names)."),
                 arg, arg, arg), call. = FALSE)
  }
  invisible(data_traits)
}

# Returns the metaweb as a numeric matrix, or stops with an informative error.
check_metaweb <- function(MW, arg = "MW") {
  if (is.data.frame(MW)) MW <- data.matrix(MW)
  if (!is.matrix(MW) || !is.numeric(MW)) {
    stop(sprintf("`%s` must be a numeric matrix.", arg), call. = FALSE)
  }
  if (nrow(MW) != ncol(MW)) {
    stop(sprintf("`%s` must be square (prey in rows, predators in columns); it is %d x %d.",
                 arg, nrow(MW), ncol(MW)), call. = FALSE)
  }
  if (is.null(rownames(MW)) || is.null(colnames(MW))) {
    stop(sprintf("`%s` must have taxon names as row and column names.", arg), call. = FALSE)
  }
  if (!identical(rownames(MW), colnames(MW))) {
    stop(sprintf("Row and column names of `%s` must be identical and in the same order.", arg),
         call. = FALSE)
  }
  if (anyDuplicated(rownames(MW)) > 0) {
    stop(sprintf("`%s` has duplicated taxon names.", arg), call. = FALSE)
  }
  if (anyNA(MW)) {
    stop(sprintf("`%s` contains missing values.", arg), call. = FALSE)
  }
  if (any(MW < 0 | MW > 1)) {
    stop(sprintf("Values of `%s` must be probabilities or 0/1 links (between 0 and 1).", arg),
         call. = FALSE)
  }
  MW
}

check_taxa_in_traits <- function(taxa, data_traits, arg = "data_traits") {
  missing_taxa <- setdiff(taxa, rownames(data_traits))
  if (length(missing_taxa) > 0) {
    stop(sprintf("%d taxa of the metaweb are not row names of `%s`: %s.",
                 length(missing_taxa), arg, format_taxa(missing_taxa)), call. = FALSE)
  }
}

warn_missing_trait <- function(traits, col, correction) {
  na_taxa <- rownames(traits)[is.na(traits[[col]])]
  if (length(na_taxa) > 0) {
    warning(sprintf("%d taxa have a missing `%s` and are ignored by the %s correction(s): %s.",
                    length(na_taxa), col, correction, format_taxa(na_taxa)), call. = FALSE)
  }
}

format_taxa <- function(x, n = 5) {
  out <- paste(head(x, n), collapse = ", ")
  if (length(x) > n) out <- paste0(out, sprintf(" and %d more", length(x) - n))
  out
}

check_count <- function(x, arg, min = 1) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x != round(x)) {
    stop(sprintf("`%s` must be a whole number >= %d.", arg, min), call. = FALSE)
  }
}

# Model parameters as a named numeric vector c(a0, a1, b0, b1). Accepts a file
# written by metaweb_mod_parameters(), a numeric vector, a one-column data
# frame (read.table() output) or the list returned by metaweb_mod_parameters().
as_pars <- function(pars) {
  if (is.list(pars) && !is.data.frame(pars) && "calibrated_data" %in% names(pars)) {
    pars <- pars$calibrated_data
  }
  if (is.character(pars) && length(pars) == 1) {
    if (!file.exists(pars)) {
      stop(sprintf("Parameter file '%s' does not exist.", pars), call. = FALSE)
    }
    pars <- read.table(pars)
  }
  if (is.data.frame(pars)) pars <- pars[[1]]
  if (!is.numeric(pars) || length(pars) != 4 || anyNA(pars)) {
    stop("Model parameters must be four numbers: a0, a1, b0, b1.", call. = FALSE)
  }
  setNames(as.numeric(pars), c("a0", "a1", "b0", "b1"))
}

# lapply() or forked parallel::mclapply(), with errors from workers re-raised
run_parallel <- function(X, FUN, mc.cores = 1) {
  check_count(mc.cores, "mc.cores")
  if (mc.cores > 1 && .Platform$OS.type == "windows") {
    warning("`mc.cores` > 1 is not supported on Windows; running on one core.", call. = FALSE)
    mc.cores <- 1
  }
  if (mc.cores == 1) return(lapply(X, FUN))

  res <- parallel::mclapply(X, FUN, mc.cores = mc.cores)
  failed <- vapply(res, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop("Parallel computation failed: ",
         conditionMessage(attr(res[[which(failed)[1]]], "condition")), call. = FALSE)
  }
  res
}
