#===============================================================================
# Global search
#===============================================================================
source("./R_scripts/server_scripts/global_search_mdl2.R")

#===============================================================================
# Beta profile
#===============================================================================

filter(ll_obj$ll_results) %>% 
  filter(loglik > max(loglik)- 100) %>%
  sapply(range) -> box

set.seed(7936910)

profile_design(
  zeta  = seq(0.75, 1.25,length = 40),
  lower = box[1, c("P_0" , "tau")],
  upper = box[2, c("P_0", "tau")],
  nprof = 15, type = "runif"
) -> guesses

ptb <- rw.sd(P_0 = ivp(0.02), tau = 0.02)
fn  <- "./Saved_objects/Irish_data/SEI3R_GBM/ifp_mdl2_beta_2.rds"

source("./R_scripts/likelihood_funs.R")

ifp_obj <- iter_filt_profile(mf1 = mf1, 
                             guesses = guesses,
                             fixed_params = fixed_params,
                             perturbations = ptb,
                             filename = fn,
                             seed = 95591631)
