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
# Global search likelihood
#===============================================================================
source("./R_scripts/helpers.R")
fn_gs     <- file.path(folder, "Global_search_mdl2.rds")
gs_obj <- readRDS(fn_gs)
mifs_global <- extract_mif_results(gs_obj)

source("./R_scripts/likelihood_funs.R")
fn          <- file.path(folder, "Global_search_mdl2_ll.rds" )
ll_obj      <- mif_ll(mifs_global, Np = 100000, 983707478, 
                      n_cores = 7, fn)
