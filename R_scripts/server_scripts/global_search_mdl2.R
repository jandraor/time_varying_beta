#===============================================================================
# Local search
#===============================================================================
source("./R_scripts/server_scripts/local_search_mdl2.R")

#===============================================================================
# Global search
#===============================================================================
mf1          <- ls_obj$result # Comes from previous script
fixed_params <- par_obj$fixed # Comes from previous script

set.seed(58741564)

runif_design(
  lower = c(zeta = 0.1, P_0 = 1,  tau = 0.05),
  upper = c(zeta = 3,   P_0 = 30, tau = 3),
  nseq  = 300
) -> guesses

fn     <- "./Saved_objects/Irish_data/SEI3R_GBM/test_Global_search2.rds"
seed   <- 1270401374
source("./R_scripts/global_search.R")
gs_obj <- global_search(guesses[1:14], fixed_params, mf1, fn, seed, 7)