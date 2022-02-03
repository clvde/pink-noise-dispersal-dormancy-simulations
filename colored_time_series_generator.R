# ---------------------------------------- #
# ---- Coloured time series generator! --- #
# ---------------------------------------- #

# Author: Courtney Van Den Elzen
# Date: Dec 3rd, 2021

# Description: This script simulates approximate colored noise time series

# Libraries

library(tidyverse)
library(tuneR)
library(psd)

createLatticeValues <- function(lattice_length = 20, alpha_val = 0, total_gens = 100){
  
  # initialize matrix
  lattice_vals <- matrix(nrow = (lattice_length^2), ncol = (total_gens+3))
  
  # input row and column numbers
  lattice_vals[,1] <- rep(1:lattice_length, times = lattice_length)
  lattice_vals[,2] <- rep(1:lattice_length, each = lattice_length)
  
  # input full lattice position label
  lattice_vals[,3] <- paste("(", 
                            rep(1:lattice_length, times = lattice_length), 
                            ", ", 
                            rep(1:lattice_length, each = lattice_length), 
                            ")", 
                            sep = "")
  
  # fill all rows with time series values for each lattice position
  for (i in 1:(lattice_length^2)) {
    
    lattice_vals[i,(1+3):(total_gens+3)] <- tuneR::noise(kind = "power", 
                                                     alpha = alpha_val, 
                                                     duration = total_gens, 
                                                     xunit = "samples")@left %>% round(digits = 4)
    
  }
  
  # convert to data frame, name columns, and convert some data types to numeric
  lattice_vals <- as.data.frame(lattice_vals)
  colnames(lattice_vals) <- c("lattice_row", "lattice_col", "lattice_position", paste("t", 1:total_gens, sep = ""))
  lattice_vals <- lattice_vals %>% mutate_at(paste("t", 1:total_gens, sep = ""), as.numeric)
  
  return(lattice_vals)
  
}

    