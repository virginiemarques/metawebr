#' Build a site x taxon matrix from detection data
#'
#' Converts long detection data (one row per site, taxon and count, e.g. eDNA
#' reads) to a site x taxon matrix, summing counts, and adds the
#' `PrimaryProducer` and `SecondaryProducer` columns, present at every site.
#'
#' @param df Data frame of detections.
#' @param station_col,species_col,count_col Names of the site, taxon and
#'   (numeric) count columns.
#' @param convert_binary Convert counts to presence (1) / absence (0)?
#'
#' @return A numeric matrix with sites in rows and taxa, followed by the two
#'   producers, in columns. Rows with a missing site or taxon are removed with
#'   a warning.
#'
#' @examples
#' det <- data.frame(station = c("s1", "s1", "s2"),
#'                   taxon = c("Gadus_morhua", "Clupea_harengus", "Gadus_morhua"),
#'                   reads = c(120, 30, 55))
#' make_species_station_matrix(det, station_col = "station",
#'                             species_col = "taxon", count_col = "reads")
#'
#' @importFrom rlang .data
#' @export
make_species_station_matrix <- function(df, station_col = "ID.Number", species_col = "taxon",
                                        count_col = "count", convert_binary = TRUE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  required_cols <- c(station_col, species_col, count_col)
  if (!all(required_cols %in% names(df))) {
    stop("Data frame must contain columns: ", paste(required_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(df[[count_col]])) {
    stop(sprintf("Column `%s` must be numeric.", count_col), call. = FALSE)
  }
  bad <- is.na(df[[station_col]]) | df[[station_col]] == "" |
    is.na(df[[species_col]]) | df[[species_col]] == ""
  if (any(bad)) {
    warning(sprintf("%d rows with a missing site or taxon were removed.", sum(bad)), call. = FALSE)
    df <- df[!bad, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("`df` has no detection left.", call. = FALSE)
  }

  df_matrix <- df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(c(station_col, species_col)))) |>
    dplyr::summarise(count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = dplyr::all_of(species_col), values_from = "count",
                       values_fill = 0)

  mat <- as.matrix(df_matrix[, -1])
  rownames(mat) <- as.character(df_matrix[[station_col]])
  mat <- cbind(mat, PrimaryProducer = 1, SecondaryProducer = 1)

  if (convert_binary) {
    mat[mat > 0] <- 1
  }
  mat
}
