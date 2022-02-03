# ---------------------------------------- #
# ---- Lattice patch distance calculator --- #
# ---------------------------------------- #

# Author: Courtney Van Den Elzen
# Date: Jan 3rd, 2022

# Description: This script takes lattice dimensions (must be square lattice) and 
# calculates each pairwise distance among patches

# Libraries
library(tidyverse)

calcPatchDistances <- 
  
  function(lattice_length = 20){
    
    row_col_values <- data.frame(row_from = 1:lattice_length, 
                                 col_from = 1:lattice_length, 
                                 row_to = 1:lattice_length, 
                                 col_to = 1:lattice_length)
    
    pairwise_value_combos <- tidyr::expand(row_col_values, row_from, col_from, row_to, col_to)
    
    pairwise_value_combos$point_from <- paste("(", pairwise_value_combos$row_from, ", ", 
                                              pairwise_value_combos$col_from, ")", 
                                              sep = "")
    
    pairwise_value_combos$point_to <- paste("(", pairwise_value_combos$row_to, ", ", 
                                            pairwise_value_combos$col_to, ")", 
                                            sep = "")
    
    pairwise_value_combos <- 
      pairwise_value_combos %>% 
      mutate(distance = 
               case_when(point_from == point_to ~ 0, 
                         point_from != point_to ~ sqrt((row_to - row_from)^2 + (col_to - col_from)^2)))
    
    
    return(pairwise_value_combos)
    
  }



