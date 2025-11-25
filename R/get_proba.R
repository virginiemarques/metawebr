#' Title: correct_metaweb_proba
#'
#' correct the metaweb based on the proba
#'
#' @description
#'
#'
#' @param MW Adjency matrix (or metaweb)
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
#'
#' @export


correct_metaweb_proba <- function(MW, df_traits){

  t <- get_proba(MW, df_traits)
  message(sprintf("Threshold value chosen is %.1f", t))

  MW[MW<t] <- 0
  MW[MW>=t] <- 1

  return(MW)
}


#' Title: get_proba
#'
#' Find the best proba to validate for the adjency matrix
#'
#' @description
#'
#'
#' @param MW Adjency matrix (or metaweb)
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
#' @importFrom NetIndices TrophInd
#'
#' @export

get_proba <- function(MW, df_traits, plot_diag = TRUE) {

  # Helper function find local minimum
  find_local_min <- function(df, value_col = "TLm") {
    if (!value_col %in% names(df)) {
      stop(paste("Column", value_col, "not found in the data frame"))
    }
    min_index <- which.min(df[[value_col]])
    df[min_index, , drop = FALSE]
  }

  # Threshold exploration
  thres <- seq(0.1, 0.975, 0.025)
  thres <- c(thres, 0.99, 0.995)
  TLm <- c()

  for (t in thres) {
    bin_net <- round(MW, 4)
    bin_net[bin_net < t]  <- 0
    bin_net[bin_net >= t] <- 1

    Troph <- NetIndices::TrophInd(Flow = bin_net, Tij = t(bin_net))
    sp <- which(is.element(rownames(Troph), rownames(df_traits)))
    TL <- abs(Troph$TL[sp] - df_traits$TrophicLevel[sp])
    TLm <- c(TLm, sum(TL, na.rm = TRUE))
  }

  df <- data.frame(TLm, thres)
  df_min <- find_local_min(df)[,2]

  # Optional plot
  if (isTRUE(plot_diag)) {
    # Wrap in print to force evaluation in all contexts
    print(
      plot(df$thres, df$TLm, pch = 1, xlab = "Binary threshold",
           ylab = "Sum of (observed TL - inferred TL)")
    )
    points(df$thres[which.min(df$TLm)], min(df$TLm), col = "red", pch = 20)
    abline(h = min(df$TLm), col = "red", lty = 2, lwd = 0.5)
  }

  return(df_min)
}
