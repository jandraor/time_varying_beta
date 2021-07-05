library(deSolve)
library(doParallel)
library(doRNG)
library(dplyr)
library(imputeTS)
library(lubridate)
library(pomp)
library(purrr)
library(readr)
library(readxl)
library(stringr)
library(tictoc)
library(tidyr)

folder <- "./Saved_objects/Irish_data/SEI3R_GBM/mdl_3" # includes immigration

source("R_scripts/get_weekly_df.R")
wkl_df <- get_weekly_df()

#===============================================================================
# Model setup
#===============================================================================

source("./R_scripts/POMP_models.R")

par_obj      <- get_params("GBM_3")
params       <- par_obj$all
fixed_params <- par_obj$fixed
pomp_mdl     <- pomp_SEI3R_GBM_immi(wkl_df, params, 1 / 128)