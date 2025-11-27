#' Title: correct_metaweb_fish
#'
#' Correct metaweb based on general ecology
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
#' @export
#'
correct_metaweb_fish <- function(data_meta,
                            data_traits,
                            herbivory = TRUE,
                            water_pos = TRUE,
                            small_fish = TRUE,
                            PP_SP_add = TRUE,
                            exception_BigFish = TRUE){

  # ******** HERBIVORY *********** #
  if(herbivory == TRUE){

    ### link correction
    sp_herb <- rownames(data_traits[which(data_traits$TrophicLevel<2.4),])

    # Message
    message(paste0("There is ", length(sp_herb), " species that are herbivorous in the dataset"))

    # Correction
    if(is.null(dim(sp_herb))){
      data_meta <- data_meta
    } else {
      for(i in 1:length(sp_herb)){data_meta[which(data_meta[,sp_herb[i]]>0),sp_herb[i]] <- 0}
      data_meta <- data_meta

      message("---------- Herbivory corrected")
    } # end if null
  } # end herbivory

  # ******** WATER POSITION *********** #

  if(water_pos == TRUE){

    Trait <-data_traits$DemersPelag
    names(Trait) <- rownames(data_traits)
    Trait <- Trait[rownames(data_meta)]

    Mw_herb_distri <- Distri_correction(data=Trait,Lniche=data_meta,mc.cores=6)
    rownames(Mw_herb_distri) <- colnames(Mw_herb_distri)<- colnames(data_meta)
    diag(Mw_herb_distri) <- 0 #No cannibalism in general
    data_meta <- Mw_herb_distri
    message("---------- Water position corrected")

  }

  # ******** SMALL FISH *********** #

  if(small_fish == TRUE){

    SM_species <- rownames(data_traits[which(data_traits$CommonLengthEstim<10),])
    message(paste0("There is ", length(SM_species), " species of small fish"))

    data_meta[,SM_species] <- 0
    message("---------- Small fish corrected")

  }

  # ******** PP and SP addition *********** #

  if(PP_SP_add == TRUE){

    # Nb herb
    sp_herb <- rownames(data_traits[which(data_traits$TrophicLevel<2.4),])
    SM_species <- rownames(data_traits[which(data_traits$CommonLengthEstim<10),])

    ## Primary producer box
    PP_prey <- rep(0,nrow(data_meta))
    names(PP_prey) <- rownames(data_meta)
    PP_prey[sp_herb] <- 1
    PP_predator <-rep(0,nrow(data_meta)+2)

    ## Secondary producer box
    sp_eps <- rownames(data_traits[which(data_traits$TrophicLevel>2.4 &data_traits$TrophicLevel<=3.7),])
    Sec_prod <- unique(c(sp_eps, SM_species)) #small fishes + low trophic level species
    PS_prey <- rep(0,nrow(data_meta))
    names(PS_prey) <- rownames(data_meta)
    PS_prey[Sec_prod] <- 1
    PS_predator <- c(rep(0,nrow(data_meta)),1,1)

    data_meta_temp <- rbind(data_meta,PrimaryProducer=PP_prey,SecondaryProducer=PS_prey)
    data_meta <- cbind(data_meta_temp ,PrimaryProducer =PP_predator,SecondaryProducer=PS_predator)

    message("---------- PP and SP corrected")

  }

  # ******** exception_BigFish *********** #

  if(exception_BigFish == TRUE){

    sp_eps <- rownames(data_traits[which(data_traits$TrophicLevel>2.4 &data_traits$TrophicLevel<3.7),])
    big_basal_sp <- sp_eps[data_traits[sp_eps,"CommonLengthEstim"]>100]
    planktivores_interact <- c(rep(0,nrow(data_meta)-2),0,1)

    for(species in big_basal_sp){
      data_meta[, species] <- planktivores_interact
    }

    message("---------- Exception - big fish corrected")

  }


  # ******** OUTPUT *********** #
  return(data_meta)

}

#' Title: Distri_correction
#'
#' Correction for water column distribution
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
#' @importFrom parallel mclapply
#'
#' @examples
#'
#'


Distri_correction  <- function(data,Lniche,mc.cores=4){

  res<- parallel::mclapply(1:dim(Lniche)[1],mc.cores=mc.cores,function(i){

    f=vector()
    for (j in 1:dim(Lniche)[1]) {

      if((i%%200)==0 & (j%%1)==400)  cat("i=",i,"j=",j,"\n")

      a <- data[i]; b <-data[j]

      f[j] = Lniche[i,j]

      if(a=="bathydemersal"&b=="pelagic-oceanic"|a=="pelagic-oceanic"&b=="bathydemersal"){ f[j] = 0}
      if(a=="bathydemersal"&b=="pelagic-neritic"|a=="pelagic-neritic"&b=="bathydemersal"){ f[j] = 0}
      if(a=="bathydemersal"&b=="reef-associated"|a=="reef-associated"&b=="bathydemersal"){ f[j] = 0}
      if(a=="bathydemersal"&b=="pelagic"|a=="pelagic"&b=="bathydemersal"){f[j] = 0}

      if(a=="bathypelagic"&b=="pelagic-oceanic"|a=="pelagic-oceanic"&b=="bathypelagic"){ f[j] = 0}
      if(a=="bathypelagic"&b=="pelagic-neritic"|a=="pelagic-neritic"&b=="bathypelagic"){ f[j] = 0}
      if(a=="bathypelagic"&b=="reef-associated"|a=="reef-associated"&b=="bathypelagic"){ f[j] = 0}
      if(a=="bathypelagic"&b=="pelagic"|a=="pelagic"&b=="bathypelagic"){ f[j] = 0}
      if(a=="pelagic-oceanic"&b=="reef-associated"|a=="reef-associated"&b=="pelagic-oceanic"){f[j] = 0}
      if(a=="pelagic-oceanic"&b=="demersal"|a=="demersal"&b=="pelagic-oceanic"){ f[j] = 0}
      if(a=="pelagic-oceanic"&b=="benthopelagic"|a=="benthopelagic"&b=="pelagic-oceanic"){ f[j] = 0}

    }#end of j
    f
  })

  do.call(rbind,res)
}
