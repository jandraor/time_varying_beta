library(doParallel)
library(doRNG)
library(dplyr)
library(pomp)
library(purrr)
library(readr)
library(stringr)
library(tictoc)


folder <- "./Saved_objects/Irish_data/SEI3R_GBM/mdl_1"
dir.create(folder, showWarnings = FALSE)

source("./R_scripts/irish_data.R")
irish_data   <- get_irish_data()

source("./R_scripts/apple.R")
drv_data_obj <- get_driving_data()

daily_df <- data.frame(time = 1:79, 
                       y1 = irish_data$y,
                       y2 = drv_data_obj$df$y2)

#===============================================================================
# Model setup
#===============================================================================

source("./R_scripts/POMP_models.R")

par_obj  <- get_params("GBM_1")
params   <- par_obj$all
pomp_mdl <- get_POMP_model(daily_df, params, 1 / 128)

source("./R_scripts/likelihood_funs.R")
source("./R_scripts/helpers.R")

path <- "./Saved_objects/Irish_data/SEI3R_GBM/weekly/mdl_2/top_10.csv"
test_pars <- read_csv(path)
pars_list <- transpose(test_pars)

imap_dfr(pars_list, function(pars_set, i) {
  n_particles <- c(5e3, 1e4, 2e4)
  fn <- file.path(folder, str_glue("pf_sensitivity_{i}.rds"))
  pars <- c(unlist(pars_set), par_obj$fixed)
  
  pomp_mdl <- pomp_mdl %>% pomp(params = pars)
  
  pf_sensitivity(n_particles = n_particles, n_cores = 7, 
                 seed = 873446145, pomp_mdl = pomp_mdl, fn = fn,
                 n_iter = 14)
}) -> sens_df