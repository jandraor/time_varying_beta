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
# Zeta profile mif
#===============================================================================
fn_gs_ll <- file.path(folder, "Global_search_mdl2_ll.rds")
ll_obj   <- readRDS(fn_gs_ll)

source("./R_scripts/helpers.R")

loglik_df <- extract_ll_df(ll_obj)

loglik_df %>% 
  filter(loglik > max(loglik)- 10, loglik.se < 2) %>%
  sapply(range) -> box

set.seed(917477792)

profile_design(
  zeta  = seq(0.8, 1.40,length = 40),
  lower = box[1, c("P_0" , "tau", "alpha")],
  upper = box[2, c("P_0", "tau", "alpha")],
  nprof = 15, type = "runif"
) -> guesses

zeta_ptb <- rw.sd(P_0 = ivp(0.02), tau = 0.02, alpha = 0.02)



fn_ls  <- file.path(folder, "local_search_mdl2.rds")
ls_obj <- readRDS(fn_ls)
mf1    <- ls_obj$result[[1]]

fn_ifp  <- file.path(folder, "ifp_zeta.rds" )

ifp_obj <- iter_filt_profile(mf1 = mf1, 
                             guesses = guesses,
                             fixed_params = fixed_params,
                             perturbations = zeta_ptb,
                             filename = fn_ifp,
                             seed = 762115017)