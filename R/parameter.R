#' Title: metaweb_mod_parameters
#'
#' Find the parameter for the metaweb
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
#' @import GenSA
#'
#' @export

metaweb_mod_parameters <- function(data_path = "inst/extdata/size_barnes2008_FF.csv", path_par_output=NULL, observation_resample=50){
  # start
  data_end <- read.csv2(data_path,dec=".",sep=",") # 34 932 in original (l) - here 5847
  MPred <- log10(data_end$standardised_predator_length)
  MPrey <- log10(data_end$si_prey_length)

  # data with unique observations
  L <- paste(data_end$predator, data_end$prey, data_end$standardised_predator_length, data_end$si_prey_length)
  L_u <- unique(L)
  ind<- c()
  for (i in 1:nrow(data_end)){
    a<-paste(data_end$predator[i], data_end$prey[i], data_end$standardised_predator_length[i],data_end$si_prey_length[i])
    if(is.element(a,L_u)){
      ind <- c(ind,i)
      L_u <- L_u[L_u != a]
    }
  }
  data_end_u <- data_end[ind,]
  MPred_u <- MPred[ind]
  MPrey_u <- MPrey[ind]

  #Reduce bias of observation : some pair species are over represented:
  data_end_u_t <- data_end[ind,]

  for (i in unique(data_end_u_t$predator)) {
    df <- data_end_u_t[data_end_u_t$predator == i, ]
    preys <- unique(df$prey)

    for (j in preys) {
      df_ij <- df[df$prey == j,]
      if (nrow(df_ij) > observation_resample) {
        r <- sample(df_ij$X,nrow(df[df$prey == j,])-observation_resample )
        data_end_u_t <- data_end_u_t[-which(is.element(data_end_u_t$X, r)),]
      }}}

  MPred_u_t <- log10(data_end_u_t$standardised_predator_length)
  MPrey_u_t <- log10(data_end_u_t$si_prey_length)
  lm_M_u_t <- lm(MPrey_u_t~MPred_u_t)

  # Calibration - parameters
  pars <- c(a0 = lm_M_u_t$coefficients[1],a1 = lm_M_u_t$coefficients[2],b0 = sd(lm_M_u_t$residuals),b1 = 0)

  ### Setting the boundaries for the algorithm of parameters estimation.
  par_lo <- c(a0 = -10, a1 = 0, b0 = -10, b1 = -10)
  par_hi <- c(a0 = 10, a1 = 10, b0 = 10, b1 = 10)

  ### Here we define the body size data that we will use to calibrate the model  ###############
  data <- data.frame(MPrey = MPrey_u_t, MPred = MPred_u_t)

  ### **Maximum likelihood estimation**

  # Model from the model_genSA script
  estim.pars = GenSA::GenSA(par = pars, fn = model, lower = par_lo, upper= par_hi, control = list(verbose =TRUE, max.time = 1000, smooth=FALSE), data = data) #Search for parameters maximizing the posteriori probability of these observed interactions

  # Save model parameter - unless output path is null
  if(!is.null(path_par_output)){
    write.table(estim.pars$par,file=path_par_output)
  }

  # Return list
  return(
    list(calibration_data=data.frame(pars),
         calibrated_data = estim.pars$par)
  )
}

#' Title: model
#'
#' Helper
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

model = function(pars,data) {
  MPred = data$MPred #list of log10 length of observed predator
  MPrey = data$MPrey
  a0 = pars[1]
  a1 = pars[2]
  b0 = pars[3]
  b1 = pars[4]
  meanMprey = mean(unique(MPrey))
  sdMprey = sd(unique(MPrey))

  # Optimum and range
  o = a0 + a1*MPred
  r = b0 + b1*MPred

  # Compute the conditional
  pLM = exp(-(o-MPrey)^2/2/r^2)

  # Compute the marginal
  pM = dnorm(x=MPrey,mean=meanMprey,sd=sdMprey)

  #  Integrate the denominator
  pL = r/(r^2+sdMprey^2)^0.5*exp(-(o-meanMprey)^2/2/(r^2+sdMprey^2))

  # Compute the posterior probability
  pML = pLM*pM/pL

  pML[pML<=0] = .Machine$double.xmin # Control to avoid computing issues

  return(-sum(log(pML)))
}

