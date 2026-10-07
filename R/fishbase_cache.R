#' Cache of FishBase tables
#'
#' Functions that read FishBase ([get_traits()], [infer_traits()],
#' [harmonise_taxa()], [get_herbivores_fishbase()], [get_diet_links()],
#' [add_depth_range()]) download each FishBase table once per R session and
#' keep it in memory. Set `options(metawebr.cache_dir = "<folder>")` to also
#' save the tables to disk, so later sessions reuse them without downloading.
#' Each saved table records the rfishbase version and the download date
#' ([fishbase_cache_info()]).
#'
#' `clear_fishbase_cache()` empties the cache, e.g. to download a newer
#' FishBase release.
#'
#' @param disk Also delete the tables saved in `getOption("metawebr.cache_dir")`?
#'
#' @return `clear_fishbase_cache()`: `NULL`, invisibly. `fishbase_cache_info()`:
#'   a data frame with one row per cached table: `table`, `rows`, `rfishbase`
#'   (version used), `retrieved` (download time) and `location` (`"memory"` or
#'   `"disk"`).
#'
#' @examples
#' \dontrun{
#' options(metawebr.cache_dir = "~/fishbase_cache")
#' traits <- get_traits(df_presence_obis)  # downloads and saves the tables
#' fishbase_cache_info()
#' clear_fishbase_cache(disk = TRUE)
#' }
#'
#' @export
clear_fishbase_cache <- function(disk = FALSE) {
  rm(list = ls(fishbase_cache), envir = fishbase_cache)
  dir <- cache_dir()
  if (isTRUE(disk) && !is.null(dir)) {
    unlink(file.path(dir, paste0("fishbase_", fishbase_tables, ".rds")))
  }
  invisible(NULL)
}

#' @rdname clear_fishbase_cache
#' @export
fishbase_cache_info <- function() {
  info <- function(x, table, location) {
    data.frame(table = table, rows = nrow(x), rfishbase = attr(x, "rfishbase") %||% NA_character_,
               retrieved = attr(x, "retrieved") %||% as.POSIXct(NA), location = location)
  }
  in_memory <- intersect(fishbase_tables, ls(fishbase_cache))
  out <- lapply(in_memory, function(t) info(get(t, envir = fishbase_cache), t, "memory"))
  dir <- cache_dir()
  if (!is.null(dir)) {
    on_disk <- setdiff(fishbase_tables[file.exists(cache_file(dir, fishbase_tables))], in_memory)
    out <- c(out, lapply(on_disk, function(t) info(readRDS(cache_file(dir, t)), t, "disk")))
  }
  if (length(out) == 0) {
    return(data.frame(table = character(), rows = integer(), rfishbase = character(),
                      retrieved = as.POSIXct(character()), location = character()))
  }
  do.call(rbind, out)
}

fishbase_cache <- new.env(parent = emptyenv())

# FishBase tables used by the package, and the rfishbase function reading each
fishbase_tables <- c("taxa", "species", "estimate", "ecology", "fooditems", "synonyms")

fishbase_download <- function(table) {
  fun <- switch(table,
                taxa = rfishbase::load_taxa,
                species = rfishbase::species,
                estimate = rfishbase::estimate,
                ecology = rfishbase::ecology,
                fooditems = rfishbase::fooditems,
                synonyms = rfishbase::synonyms)
  as.data.frame(fun())
}

cache_dir <- function() {
  dir <- getOption("metawebr.cache_dir")
  if (is.null(dir) || identical(dir, "")) NULL else path.expand(dir)
}

cache_file <- function(dir, table) file.path(dir, paste0("fishbase_", table, ".rds"))

# A FishBase table, from memory, then from disk, then downloaded
fb_table <- function(table) {
  table <- match.arg(table, fishbase_tables)
  if (exists(table, envir = fishbase_cache, inherits = FALSE)) {
    return(get(table, envir = fishbase_cache))
  }
  dir <- cache_dir()
  file <- if (!is.null(dir)) cache_file(dir, table)
  if (!is.null(file) && file.exists(file)) {
    x <- readRDS(file)
  } else {
    x <- fishbase_download(table)
    attr(x, "rfishbase") <- as.character(utils::packageVersion("rfishbase"))
    attr(x, "retrieved") <- Sys.time()
    if (!is.null(file)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(x, file)
    }
  }
  assign(table, x, envir = fishbase_cache)
  x
}

`%||%` <- function(x, y) if (is.null(x)) y else x
