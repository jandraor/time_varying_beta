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

folder <- "./Saved_objects/Irish_data/SEI3R_GBM/weekly/model_2"

source("R_scripts/get_weekly_df.R")
wkl_df <- get_weekly_df()

#===============================================================================
# Model setup
#===============================================================================

source("./R_scripts/POMP_models.R")

par_obj      <- get_params("GBM_2")
params       <- par_obj$all
fixed_params <- par_obj$fixed
pomp_mdl     <- pomp_SEI3R_GBM2(wkl_df, params, 1 / 128)

#===============================================================================
# Zeta profile log lik
#===============================================================================
fn_ifp                   <- file.path(folder, "ifp_zeta.rds" )
ifp_obj                  <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
zeta_profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, "ifp_zeta_ll_mdl2.rds")
profile_results <- mif_ll(zeta_profile_mif_results, seed = 744170621, 
                          n_cores = 7, filename = fn)
