#===============================================================================
# Libraries
#===============================================================================
library(doParallel)
library(doRNG)
library(dplyr)
library(lubridate)
library(pomp)
library(readxl)
library(stringr)
library(tictoc)
#===============================================================================
# Data
#===============================================================================

source("./R_scripts/irish_data.R")
source("./R_scripts/apple.R")

irish_data <- get_irish_data()

apple_data <- get_apple_data(start_date = "2020-02-29",
                             end_date = "2020-05-17")

driving_data <- filter(apple_data, transportation_type == "driving")
first_val    <- driving_data[1, "index"] %>% as.numeric()

driving_data <- driving_data %>%
  mutate(index = index / first_val, time = row_number())

raw_data <- driving_data$index
imp      <- na_interpolation(raw_data)

obs_df <- irish_data %>% mutate(y2 = imp) %>% 
  rename(y1 = y)

#===============================================================================
# Model setup
#===============================================================================

source("./R_scripts/POMP_models.R")

par_obj  <- get_params("2")
params   <- par_obj$all
pomp_mdl <- pomp_SEI3R_GBM2(obs_df, params)

#===============================================================================
# Local search
#===============================================================================
source("./R_scripts/local_search.R")

fn     <- "./Saved_objects/Irish_data/SEI3R_GBM/delete.rds"
ptb    <- rw.sd(zeta = 0.02, P_0 = ivp(0.02), tau = 0.02)
seed   <- 482947940
ls_obj <- local_search(pomp_mdl, params, ptb, fn, seed)

# source("./R_scripts/likelihood_funs.R")
# 
# ptb <- rw.sd(P_0 = ivp(0.02), tau = 0.02)
# fn  <- "./Saved_objects/Irish_data/SEI3R_GBM/ifp2.rds"
# 
# ifp_obj <- iter_filt_profile(mf1 = mifs_local[[1]], 
#                              guesses = guesses,
#                              fixed_params = fixed_params,
#                              perturbations = ptb,
#                              filename = fn,
#                              seed = 860523001)