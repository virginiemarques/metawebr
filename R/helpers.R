#' Title: make_species_station_matrix
#'
#' Transform from dna to matrix fitted for funtion of p/a, adds primary and secondary producers
#'
#' @description
#'
#'
#' @param x Description of the first parameter.
#' @param y Description of the second parameter (if applicable).
#' @param ... Other optional parameters passed to methods.
#'
#' @details
#'
#' @return
#' Description of the object that the function returns.
#' If the function doesn't return anything meaningful, you can say `NULL`.
#'
#' @examples
#'
#'
#' @import dplyr
#' @import tidyr
#'
#' @export

make_species_station_matrix <- function(df, station_col = "ID.Number", species_col = "taxon", count_col = "count", convert_binary = TRUE){

  # Ensure the needed columns exist
  required_cols <- c(station_col, species_col, count_col)
  if (!all(required_cols %in% names(df))) {
    stop("Data frame must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Create the wide matrix
  df_matrix <- df %>%
    group_by(across(all_of(c(station_col, species_col)))) %>%
    summarise(count = sum(.data[[count_col]], na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = {{species_col}}, values_from = count, values_fill = 0) |>
    mutate(PrimaryProducer = 1,
           SecondaryProducer = 1)

  # Convert to matrix and set row names
  mat <- as.matrix(df_matrix[,-1])
  rownames(mat) <- df_matrix[[station_col]]

  if (convert_binary) {
    mat[mat > 0] <- 1
  }

  return(mat)

}
