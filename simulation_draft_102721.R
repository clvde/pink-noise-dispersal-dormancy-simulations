# Simulation Model Draft

# Author: Courtney Van Den Elzen
# Date: October 2021

# Libraries ---------------------------------------------------------------

# Data wrangling
library(tidyverse)

# Generation of times series for spatial lattice
source("./colored_time_series_generator.R")

# Generation of lattice patch distances (distances from center points of each lattice patch)
source("./patch_distance_calculator.R")

# All functions required for reproduction
source("./reproduction_functions.R")

# All other helper functions 
source("./env_quality_competition_functions.R")

# Global Variables --------------------------------------------------------

# Environment

m <- 20 # rows and columns in lattice (lattice always square)
G <- 100 # total generations
alpha <- 0 # noise color
Kmax <- 100

var_type <- "spatial"
sim_name <- paste("a0_sim_", as.character(Sys.time()), sep = "")

# Patch distances

patch_distances_df <- calcPatchDistances(lattice_length = 20)

# Patch quality matrices

# Generate time series data
a0_lattice_vals <- createLatticeValues(alpha_val = 0,
                                       lattice_length = m,
                                       total_gens = G)

# rows are the patch, cols ar the time point
a0_lattice_vals

write_csv(a0_lattice_vals, "./alpha0_whitenoise_lattice_timeseries.csv")

# Convert time series data to lattice form where rows are the lattice row, cols 
# are the lattice column, and only a single time point is represented. 
# might be worth using purrr::map to generate a list of the lattices for each time point when this is used in simulation. 
sptl_lattice_t1 <- matrix(data = a0_lattice_vals[,1], nrow = m, ncol = m)


# Population
seed_pop_N <- 100

location_row <- sample(1:m, size = seed_pop_N, replace = TRUE)

all_indivs_df <- data.frame(gen = rep(0, times = seed_pop_N),
                            individual_id = 1:seed_pop_N, 
                            location_row = sample(1:m, size = seed_pop_N, replace = TRUE),
                            location_col = sample(1:m, size = seed_pop_N, replace = TRUE),
                            p_disp = runif(n = seed_pop_N, min = 0, max = 1),
                            p_dorm = runif(n = seed_pop_N, min = 0, max = 1),
                            mean_disp_dist = runif(n = seed_pop_N, min = 0, max = 10),
                            mean_dorm_time = sample(1:3, size = seed_pop_N, replace = TRUE))


# Initialize individual counter
ind_id_counter <- 1

for (t in 1:G) {
  
  ## (1) Filter to germinants
  # the full list (includes dormant seeds) to only those who will 
  # be germinating this generation
  germinants_df <- all_indivs_df %>% dplyr::filter(gen == (t - 1))
  
  ## Remove germinants from the individuals list - should contain only dormant 
  ## seeds at this stage.
  all_indivs_df <- dplyr::filter(gen != (t-1))
  
  ## (2) COMPETITON SHOULD HAPPEN HERE
  
  #adults_df <- (some competition function)

  # count the number of adults
  n_adults <- dim(adults_df)[1]
  
  # print the reproductive population size. 
  print(paste("Reproductive population size: ", n_adults, sep = ""))
  
  ## (3) Create offspring
  
  # choose the number of offspring each mom will produce
  n_offspring_by_mom <- fecundity(mu = 0.5, sigma = 3.75, nmom = n_adults)
  
  # make a data frame with spaces for all offspring that will be produced - this saves time
  new_offsp_df <- matrix(nrow = sum(n_offspring_by_mom), ncol = 8) %>% as.data.frame
  colnames(new_offsp_df) <- colnames(adults_df)
  
  # create all offspring for all moms
  for (a in 1:n_adults) {
    
    # choose mom to focus on in this loop
    focal_mom <- adults[a,]
    
    # find the number off offspring she has
    n_offspring <- n_offspring_by_mom[a,]

    # if the number of offspring is more than 0, we need to create those individuals.
    if (n_offspring > 0) {
      
      
      for (o in 1:n_offspring){
        
        papa <- chooseMate(focal_indiv_id = focal_mom$individual_id, 
                           all_indivs_df = adults_df, 
                           distances_df = patch_distances_df)
        
        offsp_details <- createOffspring(generation = t, 
                                         individual_id = ind_id_counter, 
                                         focal_ind, 
                                         chosen_mate = papa, patch_distance_df)
        
        offsp_ <- rbind(all_indivs_df, offsp_details)
        
      }
    }
  }
}

all_indivs_df[1,]$individual_id


belh <- data.frame(A = c(1,2,3), B = c(4,5,6), C = c(7,8,9))
rbind(belh, c(10,11,12))


all_indivs_df 

blerg <- all_indivs_df %>% dplyr::filter(gen == (1 - 1))

all_indivs_df <- all_indivs_df %>% dplyr::filter(gen != (1-1))

