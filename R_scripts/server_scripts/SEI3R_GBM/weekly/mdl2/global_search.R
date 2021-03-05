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

par_obj  <- get_params("GBM_2")
params   <- par_obj$all
pomp_mdl <- pomp_SEI3R_GBM2(wkl_df, params, 1 / 128)

#===============================================================================
# Global search
#===============================================================================

fn     <- file.path(folder, "local_search_mdl2.rds")
ls_obj <- readRDS(fn)

mf1          <- ls_obj$result[[1]] 
source("./R_scripts/POMP_models.R")
par_obj      <- get_params("GBM_2")
fixed_params <- par_obj$fixed      

set.seed(76393492)

runif_design(
  lower = c(zeta = 0.1, P_0 = 1,  tau = 0.05, alpha = 0.05),
  upper = c(zeta = 3,   P_0 = 30, tau = 3, alpha = 0.3),
  nseq  = 300
) -> guesses

fn     <- file.path(folder, "Global_search_mdl2.rds")
seed   <- 971112215
source("./R_scripts/global_search.R")
gs_obj <- global_search(guesses, fixed_params, mf1, fn, seed, 7)