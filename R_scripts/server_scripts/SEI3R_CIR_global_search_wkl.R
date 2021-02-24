#===============================================================================
# Data
#===============================================================================

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

#===============================================================================
# Weekly data
#===============================================================================

wkl_inc <- irish_data %>% slice(1:77)%>%
  mutate(week = ((time - 1) %/% 7) + 1) %>% 
  group_by(week) %>% summarise(y = sum(y)) %>% 
  mutate(time = week * 7)

wkl_df <- wkl_inc %>% rename(y1 = y) %>% 
  mutate(y2 = imp[1:length(imp) %% 7 == 0])

#===============================================================================
# Model setup
#===============================================================================

source("./R_scripts/POMP_models.R")

par_obj  <- get_params("CIR")
params   <- par_obj$all
pomp_mdl <- pomp_SEI3R_CIR(wkl_df, params)

#===============================================================================
# Local search
#===============================================================================
source("./R_scripts/local_search.R")
ptb    <- rw.sd(P_0 = ivp(0.02), tau = 0.02, nu = 0.02, upsilon = 0.02, 
                zeta = 0.02, alpha = 0.02)

fn     <- "./Saved_objects/Irish_data/SEI3R_CIR/local_search_wkl.rds"
ls_obj <- local_search(pomp_mdl, params, ptb, fn, 292718669, 7)

#===============================================================================
# Global search
#===============================================================================

mf1          <- ls_obj$result[[1]] 
fixed_params <- par_obj$fixed      

set.seed(76393492)

runif_design(
  lower = c(zeta = 0.1, P_0 = 1,  tau = 0.05, alpha = 0.05, nu = 0, 
            upsilon = 0.01),
  upper = c(zeta = 3,   P_0 = 30, tau = 3, alpha = 0.3, nu = 0.5, 
            upsilon = 0.25 ),
  nseq  = 300
) -> guesses

fn     <- "./Saved_objects/Irish_data/SEI3R_CIR/test_Global_search_wkl.rds"
seed   <- 435367905
source("./R_scripts/global_search.R")
gs_obj <- global_search(guesses, fixed_params, mf1, fn, seed, 7)
