source("./R_scripts/model_setup_GBM.R")

#===============================================================================
# Global search likelihood
#===============================================================================
source("./R_scripts/helpers.R")
fn_gs       <- file.path(folder, "Global_search_mdl2.rds")
gs_obj      <- readRDS(fn_gs)
mifs_global <- extract_mif_results(gs_obj)

source("./R_scripts/likelihood_funs.R")
fn          <- file.path(folder, "Global_search_mdl2_ll.rds" )
ll_obj      <- mif_ll(mifs_global, Np = 100000, 983707478, 
                      n_cores = 7, fn)
