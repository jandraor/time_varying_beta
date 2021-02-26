library(deSolve)
library(doParallel)
library(doRNG)
library(dplyr)
library(GGally)
library(ggplot2)
library(ggpubr)
library(imputeTS)
library(lubridate)
library(pomp)
library(purrr)
library(readr)
library(readsdr)
library(readxl)
library(scales)
library(stringr)
library(tictoc)
library(tidyr)

folder <- "./Saved_objects/Irish_data/SEI3R_CIR/weekly/"

source("./R_scripts/helpers.R")

global_file <- file.path(folder, "Global_search_mdl2.rds")
gs_obj      <- readRDS(global_file)
mifs_global <- extract_mif_results(gs_obj$mf_results)

source("./R_scripts/likelihood_funs.R")
fn          <- file.path(folder, "Global_search_mdl2_ll.rds" )
ll_obj      <- mif_ll(mifs_global, Np = 100000, 1270401374, 
                      n_cores = 7, fn)