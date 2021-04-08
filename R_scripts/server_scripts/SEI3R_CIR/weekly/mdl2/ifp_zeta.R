source("./R_scripts/model_setup_CIR.R")

fn         <- file.path(folder, "local_search_mdl2.rds")
ls_obj     <- readRDS(fn)
mifs_local <- ls_obj$result

source("./R_scripts/helpers.R")

fn        <- file.path(folder, "Global_search_mdl2_ll.rds" )
ll_obj    <- readRDS(fn) 
loglik_df <- extract_ll_df(ll_obj)

source("./R_scripts/likelihood_funs.R")

loglik_df %>% 
  filter(loglik > max(loglik)- 20, loglik.se < 2) %>%
  sapply(range) -> box

set.seed(917477792)

profile_design(
  zeta  = seq(0.5, 1.75,length = 40),
  lower = box[1, c("P_0" , "tau", "alpha", "nu", "upsilon")],
  upper = box[2, c("P_0", "tau", "alpha", "nu", "upsilon")],
  nprof = 15, type = "runif"
) -> guesses

ptb <- rw.sd(P_0 = ivp(0.02), tau = 0.02, alpha = 0.02, nu = 0.02,
             upsilon = 0.02)

fn      <- file.path(folder, "ifp_zeta.rds")

ifp_obj <- iter_filt_profile(mf1           = mifs_local[[1]], 
                             guesses       = guesses,
                             fixed_params  = fixed_params,
                             perturbations = ptb,
                             filename      = fn,
                             seed          = 553239040)