source("./R_scripts/model_setup_CIR.R")

#===============================================================================
# Zeta profile 
#===============================================================================
fn_ifp  <- file.path(folder, "ifp_zeta.rds" )
ifp_obj <- readRDS(fn_ifp)

source("./R_scripts/helpers.R")
profile_mif_results <- extract_mif_results(ifp_obj)

source("./R_scripts/likelihood_funs.R")
fn              <- file.path(folder, "ifp_zeta_ll_mdl2.rds")
profile_results <- mif_ll(profile_mif_results, seed = 426820244, 
                          n_cores = 7, filename = fn)