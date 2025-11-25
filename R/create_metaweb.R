#' Title: apply_model_metaweb
#'
#' Create the metaweb
#'
#' @description
#'
#' Clean trait data. Get common lenght if missing (0.6*TL)
#'
#' @param data_traits Dataframe containing trait data (columms...)
#' @param path_pars Path to parameters dataframe
#'
#' @details
#'
#' @return
#' Description of the object that the function returns.
#' If the function doesn't return anything meaningful, you can say `NULL`.
#'
#' @examples
#'
#' @export
#'
apply_model_metaweb <- function(data_traits,
                                path_pars,
                                length_column = "CommonLengthEstim"){

  # Check
  if(identical(rownames(data_traits),seq(1:nrow(data_traits)))){
    stop("check rownames, should be taxa name")
  }

  Size <- log10(data_traits[,length_column])
  names(Size) <- rownames(data_traits)

  # Expand all pairs of species
  M <- expand.grid(Size,Size)
  MPrey <- M[,1]; MPred <- M[,2]

  Names <- expand.grid(names(Size),names(Size))
  PreyNames <- Names[,1];PredNames <- Names[,2]

  rm(M);rm(Names)

  # Compute interaction probability
  Pars <- read.table(path_pars)
  pLM <- pLMFitted(MPrey,MPred,Pars) # function extracted from Model_GenSA.R code

  # Transform into an adgency matrix
  mw <- matrix(pLM,nr = length(Size), nc = length(Size), byrow = FALSE)
  colnames(mw)<-rownames(mw)<- names(Size)

  return(mw)
}

#' Title: pLMFitted
#'
#' Helper for above
#'
#' @description
#'
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

pLMFitted = function(MPrey,MPred,Pars) {
  with(Pars, {
    a0 = Pars[1,]; a1= Pars[2,]
    b0 = Pars[3,]; b1=Pars[4,]
    o = a0 + a1*MPred
    #else if(isTemp == 1) o = a0 +(a1 + a2*Temp)*MPred + a3*Temp
    r = b0 + b1*MPred

    # Compute the conditional pLM for each predator
    exp(-(o-MPrey)^2/2/r^2)
  })
}
