#===============================================================================
# Libraries
#===============================================================================
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
library(tidyr)
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

source("./R_scripts/likelihood_funs.R")

source("./R_scripts/helpers.R")

Np_tries <- c(5e3, 1e4, 5e4, 1e5, 2.5e5, 5e5, 1e6, 2e6, 5e6)

map_df(Np_tries, function(Np) {
  tic.clearlog()
  tic()
  params["zeta"] <- 1.032
  params["P_0"]  <- 4.854
  params["tau"]  <- 0.212
  
  pomp_mdl2 <- pomp_mdl %>% pomp(params = params)
  
  
  lik_est <- estimate_avg_lik(pomp_mdl2, 8647872, 7, Np = Np, 
                              n_iter = 14)
  toc(quiet = FALSE, log = TRUE)
  log.lst <- tic.log(format = FALSE)
  
  data.frame(Np = Np, loglik = lik_est[[1]], loglik.se = lik_est[[2]],
             time = calculate_time(log.lst))
}) -> comparison_ll_df

write_csv(comparison_ll_df, 
          "./Saved_objects/Irish_data/SEI3R_GBM/likelihood_model_2.csv")

