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

folder <- "./Saved_objects/Irish_data/SEI3R_GBM/weekly"

source("./R_scripts/helpers.R")
source("./R_scripts/likelihood_funs.R")


ifp_file        <- file.path(folder, "ifp2_zeta.rds" )
ifp_obj         <- readRDS(ifp_file)

zeta_profile_mif_results <- extract_mif_results(ifp_obj)
fn              <- file.path(folder, "zeta_profile_ll_mdl2.rds")
profile_results <- mif_ll(zeta_profile_mif_results, seed = 744170621, n_cores = 7,
                          filename = fn)