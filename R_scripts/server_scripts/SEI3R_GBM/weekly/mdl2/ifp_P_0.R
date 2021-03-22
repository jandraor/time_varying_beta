source("./R_scripts/model_setup_GBM.R")

# Local search results----------------------------------------------------------
fn_ls  <- file.path(folder, "local_search_mdl2.rds")
ls_obj <- readRDS(fn_ls)
mf1    <- ls_obj$result[[1]]
#-------------------------------------------------------------------------------

# Global search results --------------------------------------------------------

fn_gs_ll <- file.path(folder, "Global_search_mdl2_ll.rds")
ll_obj   <- readRDS(fn_gs_ll)

source("./R_scripts/helpers.R")

loglik_df <- extract_ll_df(ll_obj)
#-------------------------------------------------------------------------------

# Guesses-----------------------------------------------------------------------
loglik_df %>% 
  filter(loglik > max(loglik)- 20, loglik.se < 2) %>%
  sapply(range) -> box

set.seed(416831519)

profile_design(
  P_0  = seq(0.1, 15,length = 40),
  lower = box[1, c("zeta" , "tau", "alpha")],
  upper = box[2, c("zeta", "tau", "alpha")],
  nprof = 15, type = "runif"
) -> guesses
#-------------------------------------------------------------------------------
source("./R_scripts/likelihood_funs.R")
P_0_ptb <- rw.sd(zeta = 0.02, tau = 0.02, alpha = 0.02)

fn_ifp     <- file.path(folder, "ifp_P_0.rds" )

ifp_obj <- iter_filt_profile(mf1 = mf1, 
                             guesses = guesses,
                             fixed_params = fixed_params,
                             perturbations = P_0_ptb,
                             filename = fn_ifp,
                             seed = 373492667)