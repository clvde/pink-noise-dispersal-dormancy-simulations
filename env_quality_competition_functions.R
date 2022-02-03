# --------------------------------------------------------- #
# - Environmental Quality and Patch Competition Functions - #
# --------------------------------------------------------- #
library(tidyverse)
library(statmod)


# ---- Patch Competition ---- 

## Description: Determines which individual seedlings live and die in a patch 
## with a carrying capacity of K. Individuals are randomized and thinned to K if n > K.

competitonOutcome <- function(all_indivs, patch_qual_df, timept, appx_mean_K = 25){
  
  if(dim(all_indivs)[[1]] == 0){
    
    return(all_indivs)
    
  } else{
  
  timept_name <- colnames(patch_qual_df)[3+timept]
  
  # count number of indivs in each patch
  n_indivs_by_patch <- 
    all_indivs %>% 
    dplyr::group_by(location_row, location_col, location) %>% 
    dplyr::summarize(n())
  
  # contains the Ks for each patch based on (-1,1) patch quality value from tuneR::noise()
  patch_K_df <- 
    patch_qual_df %>%
    dplyr::select(lattice_row, lattice_col, lattice_position, timept_name) %>%
    dplyr::rename(timepoint = timept_name) %>%
    dplyr::mutate(patch_K = round((timepoint - min(timepoint))*25, digits = 0)) %>%
    dplyr::select(lattice_position, patch_K) %>%
    dplyr::arrange(lattice_position)

  # randomize location of all individuals within each patch
  # then number of the randomized individuals and remove the ones that are over K
  filt_indivs_over_K_df <- 
    all_indivs %>%
    dplyr::group_by(location) %>%
    # assign each row a random number between 0 and 10 (real valued)
    dplyr::mutate(rand_place = runif(n = length(gen), min = 0, max = 10)) %>%
    # then sort within each patch by the random number value
    dplyr::arrange(location, rand_place) %>%
    # next assign each individual a number within each patch (these are NOT unique among patches)
    dplyr::mutate(randn = 1:n()) %>%
    # join with the patch K matrix
    dplyr::full_join(patch_K_df, by = c("location" = "lattice_position")) %>%
    # filter out those with a number over K
    dplyr::mutate(over_K = (patch_K - randn)) %>%
    dplyr::filter(over_K >= 0) %>%
    # remove this category since it's not important external to this function
    dplyr::select(-c(patch_K, rand_place, randn, over_K))
  
  return(filt_indivs_over_K_df)
  }
  
}

# ---- Dormancy Mortality ---- 

## Description: from dormant individuals, there is a small probability of mortality 
## each generation that the seed remains dormant.

dormancyMortality <- function(dormant_indivs, p_mort) {
  
  n_indivs <- dim(dormant_indivs)[1]
  
  if(n_indivs == 0){
    
    return(dormant_indivs)
  } else {
  
  dormant_indivs$dorm_mort <- rbinom(n = n_indivs, size = 1, prob = p_dorm_mort)
  
  dormant_indivs <- 
    dormant_indivs %>% 
    dplyr::filter(dorm_mort == 0) %>%
    dplyr::select(-dorm_mort)
  
  return(dormant_indivs)
  }
  
}


