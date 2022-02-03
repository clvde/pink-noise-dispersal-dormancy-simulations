# ---------------------------------------- #
# -------- Reproduction functions -------- #
# ---------------------------------------- #
library(tidyverse)
library(statmod)

# ---- Fecundity ---- 

## Description: Inputs: mean fecundity (parameterizes a lognormal distribution), 
##              sigma shape parameter of the lognormal distribution, and the 
##              total number of maternal plants.  Outputs: realized number of 
##              offspring for each maternal plant. This generates intrinsic 
##              demographic stochasticity. 

fecundity <- function(lam = 25, nmom){
  
  fec <- rpois(n = nmom, lambda = lam)

  return(fec)
  
}

# ---- Mate Choice ----

chooseMate <- function(all_indivs, 
                       mat_ind_id, 
                       patch_distances){
  
  # get focal individual details
  maternal_indiv <- dplyr::filter(all_indivs, individual_id == mat_ind_id)
  
  # get focal individual location
  maternal_indiv_location <- c(maternal_indiv$location_row, maternal_indiv$location_col)

  # get locations of potential mates (i.e. all other individuals)
  all_indiv_locs <-
    dplyr::filter(all_indivs, individual_id != mat_ind_id) %>%
    dplyr::select(location_row, location_col) %>% 
    dplyr::mutate(location = paste("(", location_row, ", ", location_col, ")", sep = ""))

  # get possible distances (i.e. distances of all patches from patch where focal individual is located)
  possible_dists <- 
    dplyr::filter(patch_distances, 
                  row_from == maternal_indiv_location[[1]] & col_from == maternal_indiv_location[[2]]) %>%
    dplyr::select(point_to, distance)

  # get distances of potential mates (i.e. all other individuals)
  all_indiv_dists <- dplyr::inner_join(x = all_indiv_locs, y = possible_dists, by = c("location" = "point_to"))

  ###### ----- might need to change something about the assignment of all_indiv_dists. 
  ###### ----- Currently adds a missing grouping variable (location). Wondering if it 
  ###### ----- is saving "point_to" as the joining variable and then adding location 
  ###### ----- back in?
  
  # choose mate using a single draw from a multinomial distribution - probability is proportional to distance
  choose_mate <- rmultinom(n = 1, size = 1, prob = all_indiv_dists$distance)

  # select the chosen mate and return its information.
  chosen_mate <- all_indivs[which(choose_mate == 1),]

  return(chosen_mate)
  #return(all_indiv_dists)
}


# ---- Offspring Dispersal Outcome ----

## Description: Inputs: maternal plant, patch distance matrix 
##              Outputs: dispersal patch outcome for one offspring. If offspring
##                       die during dispersal, the natal patch is returned and 
##                       the individual is removed within createOffspring().




## NEED TO ADD FUNCTIONALITY FOR DISPERSAL MORTALITY
dispersalOutcome <- function(mat_indiv, patch_distances, p_mort) {
  
  # Dispersal outcome of offspring - yes/no, based on maternal plant p_disp (probability of dispersal)
  yn <- rbinom(n = 1, size = 1, prob = mat_indiv$p_disp)
  
  mort <- rbinom(n = 1, size = 1, prob = p_mort)
  
  if (yn){ # if dispersing
    
    if (!mort) { # if it didn't die during dispersal
      
      # distance traveled
      drawn_dist <- 
        rinvgauss(n = 1, 
                  mean=mat_indiv$mean_disp_dist, 
                  shape=13)
      
      # distances of all other patches from focal indiv patch, add column of 
      # differences between the drawn travel distance and the distances of all patches. 
      # whichever patch(es) have the smallest distance discrepancy will be chosen 
      # among randomly for the disperse-to patch.
      patch_choice_df <- 
        patch_distances %>%
        dplyr::mutate(distance_diff = abs(distance - drawn_dist)) %>%
        dplyr::arrange(by = distance_diff)
      
      # find the minimum distance discrepancy
      min_diff <- patch_choice_df$distance_diff %>% min
      
      # filter to only the minimum distance discrepancy patches
      rand_patch_choice_df <- 
        patch_choice_df %>%
        dplyr::filter(distance_diff <= min_diff)
      
      # sample one of the patches from rand_patch_choice_df
      chosen_transition <- rand_patch_choice_df[sample(x = 1:(dim(rand_patch_choice_df)[1]), size = 1),]
      
      # find the row nad column numbers of that patch
      chosen_patch <- c(chosen_transition$row_to, chosen_transition$col_to)
      
      return(chosen_patch)
      
    } else if (mort) { # if it did die during dispersal
      
      return(c(NA,NA))
    }
  } else { # survives but doesn't disperse
    return(c(mat_indiv$location_row, mat_indiv$location_col))
  }
}

# ---- Offspring Dormancy Outcome ----

## Description: Inputs: offspring probability of dormancy (average of parental 
##              values), offspring mean number of dormant generations. 
##              Outputs: number of generations spent in dormancy.

dormancyOutcome <- function(pdorm_offsp, mean_time_offsp) {
  
  # Dormancy outcome - yes/no, based on offspring p_dorm (probability of dormancy)
  yn <- rbinom(n = 1, size = 1, prob = pdorm_offsp)
  
  
  if (yn){ # if dormant
    
    # drawn the number of generations to stay dormant (using mean dorm time value for offspring)
    dorm_ngens <- rpois(n = 1, lambda = mean_time_offsp)
    
    return(dorm_ngens)
    
  } else {
    
    # return 0 as zero generations dormant
    return(0)
    
  }
}

# ---- Create Offspring ----

## Description: Inputs: generation number, individual id, focal (maternal) 
##              individual, chosen mater (paternal individual), matrix with all 
##              patch distances.
##              Outputs:

createOffspring <- function(generation, individual_id, maternal_ind, chosen_mate, patch_distances, p_disp_mort) {
  
  # (1) set up offspring row
  offsp_values <- matrix(data = NA, nrow = 1, ncol = 9) %>% as.data.frame
  colnames(offsp_values) <- c("gen", "individual_id", "location_row", 
                              "location_col", "p_disp", "p_dorm", 
                              "mean_disp_dist", "mean_dorm_time", "location")
  
  # (2) draw dispersal outcome - returns either mortality marker or resultant patch
  disp_out <- dispersalOutcome(mat_indiv = maternal_ind, 
                               patch_distances = patch_distances,
                               p_mort = p_disp_mort)
  
  if ((is.na(sum(disp_out)))){ # the offspring died during dispersal

    return(offsp_values)
    
  } 
  else { # the offspring lived, fill in the rest of the data
    
    offsp_values[, c("location_row", "location_col")] <- disp_out
    offsp_values$location <- paste("(", offsp_values$location_row, ", ", offsp_values$location_col, ")", sep = "")
    
    # (2) assign individual ID
    offsp_values$individual_id <- individual_id
    
    # (4) calculate offspring dispersal phenotypes
    offsp_values[, c("p_disp", "mean_disp_dist")] <- c(mean(c(maternal_ind$p_disp, chosen_mate$p_disp)),
                                                       mean(c(maternal_ind$mean_disp_dist, chosen_mate$mean_disp_dist))) 
    
    # (5) calculate offspring dormancy phenotypes
    offsp_values[,"p_dorm"] <- mean(c(maternal_ind$p_dorm, chosen_mate$p_dorm))
    offsp_values[,"mean_dorm_time"] <- mean(c(maternal_ind$mean_dorm_time, chosen_mate$mean_dorm_time))
    
    # (6) draw dormancy outcome (yes/no) and, if yes, n generations spent in dormancy
    realized_dorm_time <- dormancyOutcome(pdorm_offsp = offsp_values$p_dorm, mean_time_offsp = offsp_values$mean_dorm_time)
    
    # (7) calculate the realized germination generation for this seed.
    offsp_values$gen <- generation + realized_dorm_time
    
    
    return(offsp_values)
  }
}
  
  
  
  
  