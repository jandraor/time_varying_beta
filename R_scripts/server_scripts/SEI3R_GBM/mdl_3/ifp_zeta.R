source("./R_scripts/model_setup_GBM_3.R")

#===============================================================================
var_name <- "zeta"
#===============================================================================

# Local search results----------------------------------------------------------
fn_ls  <- file.path(folder, "local_search.rds")
ls_obj <- readRDS(fn_ls)
mf1    <- ls_obj$result[[1]]
#-------------------------------------------------------------------------------

# Global search results --------------------------------------------------------

fn_gs_ll <- file.path(folder, "Global_search_ll.rds")
ll_obj   <- readRDS(fn_gs_ll)

source("./R_scripts/helpers.R")

loglik_df <- extract_ll_df(ll_obj)
#-------------------------------------------------------------------------------

# Guesses-----------------------------------------------------------------------
loglik_df %>% 
  filter(loglik > max(loglik)- 20, loglik.se < 2) %>%
  sapply(range) -> box

set.seed(167826283)

other_pars <- unk_names[unk_names != var_name]

profile_design(
  zeta  = seq(0.4, 2,length = 50),
  lower = box[1, other_pars],
  upper = box[2, other_pars],
  nprof = 30, type = "sobol"
) -> guesses
#-------------------------------------------------------------------------------
source("./R_scripts/likelihood_funs.R")
zeta_ptb <- rw.sd(P_0 = ivp(0.02), tau = 0.02, alpha = 0.02, iota = 0.02)

fn_ifp     <- file.path(folder, "ifp_zeta.rds" )

ifp_obj <- iter_filt_profile(mf1 = mf1, 
                             guesses = guesses,
                             fixed_params = fixed_params,
                             perturbations = zeta_ptb,
                             filename = fn_ifp,
                             seed = 105449832)