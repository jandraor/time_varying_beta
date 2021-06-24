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

set.seed(358559280)

profile_design(
  tau   = seq(0, 0.50,length = 50),
  lower = box[1, c("zeta" , "alpha", "P_0", "nu", "upsilon")],
  upper = box[2, c("zeta", "alpha", "P_0", "nu", "upsilon")],
  nprof = 30, type = "sobol"
) -> guesses_tau

ptb <- rw.sd(P_0 = ivp(0.02), alpha = 0.02, zeta = 0.02, nu = 0.02,
             upsilon = 0.02)

fn      <- file.path(folder, "ifp_tau.rds")

ifp_obj <- iter_filt_profile(mf1 = mifs_local[[1]], 
                             guesses = guesses_tau,
                             fixed_params = fixed_params,
                             perturbations = ptb,
                             filename = fn,
                             seed = 139567891)