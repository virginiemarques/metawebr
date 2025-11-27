#' Title: clean_traits
#'
#' Clean trait data
#'
#' @description
#'
#' Clean trait data. Get common lenght if missing (0.6*TL)
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
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when distinct group_by summarise left_join arrange slice select ungroup
#'
#' @export
#'
clean_traits <- function(data_traits,
                         data_presence=NULL,
                         TL_clean=FALSE){

  message("There is ", sum(is.na(data_traits$CommonLengthEstim)), " missing values for CL out of ", nrow(data_traits))
  message("There is ", sum(is.na(data_traits$TrophicLevel)), " missing values for TL out of ", nrow(data_traits))
  message("There is ", sum(is.na(data_traits$DemersPelag)), " missing values for Trophic Level out of ", nrow(data_traits))

  # --------------------------------------- #
  # Clean the common length if absent
  data_traits <- data_traits %>%
    mutate(Common_length = case_when(
      is.na(Common_length) ~ 0.60 * Max_length,
      TRUE ~ Common_length
    ))

  # --------------------------------------- #
  # Clean the TL if missing data
  # to deal with later .... what does this do??

  if(TL_clean == TRUE){
    names_na_length <- rownames(data_traits[which(is.na(data_traits$Common_length)),])

    if(length(names_na_length) > 0){
      species_length <- list()

      for(i in 1:length(names_na_length)){
        cat("i",i,"\n")
        d <-data_traits[grep(do.call(rbind,strsplit(names_na_length,"_"))[i,1],rownames(data_traits)),"Common_length"]
        species_length[[i]] <- median(d,na.rm=TRUE)
      } # end of i

      names(species_length) <- names_na_length
      species_length <- do.call(rbind,species_length)

      data_traits[rownames(species_length),"Common_length"]  <- species_length[,1]
      data_traits <- data_traits[-which(is.na(data_traits[,"Common_length"])),]
    } # end check is 0 for TL completion
  } # end TL clean check

  # --------------------------------------- #
  # By co-occurence (not done yet)


  return(data_traits)
}

#' Title: infer_traits
#'
#' Infer trait data when missing
#'
#' @description
#'
#' Clean trait data. Get common lenght if missing (0.6*TL)
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
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate case_when distinct group_by summarise left_join arrange slice select ungroup
#'
#' @export
#'

infer_traits <- function(data_traits,
                         data_presence,
                         method = c("mean_level", "100matrix")){

  # --------------------------------------- #
  # If not at the species level, either make a mean by genus/family
  # Or generate 100 traits matrices

  data_presence_uniq <- data_presence %>% distinct(species_name_corrected, genus_name_corrected, family_name_corrected, rank_ncbi_corrected)
  # add the traits to this one

  if(method == "mean_level"){

    data_traits_list <- split(data_traits, data_traits$rank)

    # ------------------------------ #
    # by species
    sp_all <- data_traits_list$species %>%
      distinct(taxa_edna, Common_length, Trophic_Level, Env_2)

    # ------------------------------ #
    # by genus: get the mean CL for each genus

    genus_CL <- data_traits_list$genus %>%
      group_by(taxa_edna) %>%
      summarise(Common_length = mean(Common_length, na.rm = TRUE))

    genus_trophic <- data_traits_list$genus %>%
      group_by(taxa_edna) %>%
      summarise(Trophic_Level = mean(Trophic_Level, na.rm = TRUE))

    genus_Env <- data_traits_list$genus[,c("taxa_edna", "Env_2")] %>% # this line is not robust at all
      group_by(taxa_edna, Env_2) %>%
      summarize(n = n()) %>% arrange(desc(n)) %>% slice(1) %>% dplyr::select(-n) %>%  ungroup()

    genus_all <- genus_CL %>%
      left_join(., genus_trophic) %>%
      left_join(., genus_Env)

    # ------------------------------ #
    # by family: get the mean CL for each genus
    family_CL <- data_traits_list$family %>%
      group_by(taxa_edna) %>%
      summarise(Common_length = mean(Common_length, na.rm = TRUE))

    family_trophic <- data_traits_list$family %>%
      group_by(taxa_edna) %>%
      summarise(Trophic_Level = mean(Trophic_Level, na.rm = TRUE))

    family_Env <- data_traits_list$family[,c("taxa_edna", "Env_2")] %>%
      group_by(taxa_edna, Env_2) %>%
      summarize(n = n()) %>% arrange(desc(n)) %>% slice(1) %>% dplyr::select(-n) %>%  ungroup()

    family_all <- family_CL %>%
      left_join(., family_trophic) %>%
      left_join(., family_Env)

    # ------------------------------ #
    # Combine
    df <- rbind(sp_all, genus_all, family_all)
    rownames(df) <- df$taxa_edna
    return(df)
  } else { if(method == "100matrix"){

    # voir code romane TBI

  } # end of if 100matrix
  } # end of else
} # end of function
