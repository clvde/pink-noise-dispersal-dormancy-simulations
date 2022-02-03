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
G <- 25 # total generations
alpha <- 0 # noise color
Kmean <- 25
p_dorm_mort <- 0.05
p_disp_mort <- 0.05

var_type <- "spatial"
sim_name <- paste("a0_sim_", as.character(Sys.time()), sep = "")

# Patch distances

patch_distances_df <- calcPatchDistances(lattice_length = 20)
patch_distances_df 
# Patch quality matrices

# Generate time series data
a0_lattice_vals <- createLatticeValues(alpha_val = 0,
                                       lattice_length = m,
                                       total_gens = G) 

write_csv(a0_lattice_vals, "./alpha0_whitenoise_lattice_timeseries.csv")

# Convert time series data to lattice form where rows are the lattice row, cols 
# are the lattice column, and only a single time point is represented. 
# might be worth using purrr::map to generate a list of the lattices for each 
# time point when this is used in simulation. 
sptl_lattice_t1 <- matrix(data = a0_lattice_vals[,1], nrow = m, ncol = m)

# Population
seed_pop_N <- 1000

location_row <- sample(1:m, size = seed_pop_N, replace = TRUE)

init_indivs_df <- 
  data.frame(gen = rep(0, times = seed_pop_N),
             individual_id = 1:seed_pop_N, 
             location_row = sample(1:m, size = seed_pop_N, replace = TRUE),
             location_col = sample(1:m, size = seed_pop_N, replace = TRUE),
             p_disp = runif(n = seed_pop_N, min = 0, max = 1),
             p_dorm = runif(n = seed_pop_N, min = 0, max = 1),
             mean_disp_dist = runif(n = seed_pop_N, min = 0, max = 10),
             mean_dorm_time = sample(1:3, size = seed_pop_N, replace = TRUE)) %>%
  dplyr::mutate(location = paste("(", location_row, ", ", location_col, ")", sep = ""))


init_indivs_df

colnames(init_indivs_df)

# Initialize individual counter
ind_id_counter <- dim(init_indivs_df)[[1]]

sim_start_time <- Sys.time()

for (t in 1:10) {
  
  gen_start_time <- Sys.time()
  
  if (t == 1) {
    # 9 columns
    all_indivs_df <- init_indivs_df
  }
  
  ## (1) GERMINATION
  # the full list (includes dormant seeds) to only those who will 
  # be germinating this generation
  
  # 10 columns (+ rand_place)
  germinants_df <- 
    all_indivs_df %>% 
    dplyr::filter(gen == (t - 1)) %>%
    dplyr::mutate(rand_place = runif(n = length(gen), min = 0, max = 10)) %>%
    dplyr::arrange(rand_place) %>%
    dplyr::select(-rand_place)
  
  
  ## Remove germinants from the individuals list - should contain only dormant
  ## seeds at this stage.
  
  # 9 columns
  dorm_new_indivs_df <-
    all_indivs_df %>%
    dplyr::filter(gen != (t-1))

  ## (2) COMPETITON
  # 9 columns
  adults_df <-
    competitonOutcome(all_indivs = germinants_df,
                      patch_qual_df = a0_lattice_vals,
                      timept = t,
                      appx_mean_K = 25)

  ## (3) DORMANCY MORTALITY
  
  dorm_new_indivs_df <-
    dormancyMortality(dormant_indivs = dorm_new_indivs_df,
                      p_mort = p_dorm_mort)
  
  
  ## (4) REPRODUCTION

  # count the number of adults
  n_adults <- dim(adults_df)[1]

  # print the reproductive population size.
  print(paste("Reproductive population size: ", n_adults, sep = ""))

  # choose the number of offspring each mom will produce
  n_offspring_by_mom <- fecundity(lam = 25, nmom = n_adults)

  # make a data frame with spaces for all offspring that will be produced - this saves time
  new_offsp_df <- matrix(nrow = sum(n_offspring_by_mom), ncol = 9) %>% as.data.frame
  colnames(new_offsp_df) <- colnames(adults_df)

  offsp_counter <- 0

  # create all offspring for all moms
  for (a in 1:n_adults) {

    # choose mom to focus on in this loop
    focal_mom <- adults_df[a,]

    # find the number off offspring she has
    n_offspring <- n_offspring_by_mom[[a]]

    # if the number of offspring is more than 0, we need to create those individuals.
    if (n_offspring > 0) {

      for (o in 1:n_offspring){

        #counts offspring produced each generation for placement in the data frame
        offsp_counter <- offsp_counter + 1

        # gives unique individual ID - corresponds with timing of creation
        ind_id_counter <- ind_id_counter + 1
        
        print("Here is where the location added message is coming from")
        papa <- chooseMate(all_indivs = adults_df,
                           mat_ind_id = focal_mom$individual_id,
                           patch_distances = patch_distances_df)

        offsp_details <- createOffspring(generation = t,
                                         individual_id = ind_id_counter,
                                         maternal_ind = focal_mom,
                                         chosen_mate = papa,
                                         patch_distances = patch_distances_df,
                                         p_disp_mort = p_disp_mort)
        
        new_offsp_df[offsp_counter, ] <- offsp_details
      }
    } else {
      print("no offspring produced")

    }

  }

  # remove any seeds that died during dispersal
  new_offsp_df <-
    new_offsp_df %>%
    dplyr::filter(!is.na(individual_id))
  
  # add the surviving seeds onto the all_indivs_df list
  dorm_new_indivs_df <- dplyr::bind_rows(dorm_new_indivs_df, new_offsp_df)

  if(t%%1 == 0) {

    write_csv(dorm_new_indivs_df, paste("./sim_outputs/all_seeds_end_gen", t, ".csv", sep = ""))

  }
  
  print("Gen ", t, " Complete.", sep = "")
  
  gen_end_time <- Sys.time()
  
  print(gen_end_time - gen_start_time)
  
}
sim_end_time <- Sys.time()

sim_end_time - sim_start_time

