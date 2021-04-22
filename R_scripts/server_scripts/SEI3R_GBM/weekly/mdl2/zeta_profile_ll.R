source("./R_scripts/model_setup_GBM.R")

var_name <- "zeta"

#===============================================================================
# Zeta profile log lik
#===============================================================================
fn_ifp                   <- file.path(folder, "ifp_zeta.rds" )
ifp_obj                  <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
zeta_profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, "ifp_zeta_ll_mdl2.rds")
profile_results <- mif_ll(zeta_profile_mif_results, seed = 744170621, 
                          n_cores = 7, filename = fn)
