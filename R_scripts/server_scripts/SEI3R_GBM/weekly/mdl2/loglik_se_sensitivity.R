source("./R_scripts/model_setup_GBM.R")

source("./R_scripts/likelihood_funs.R")
source("./R_scripts/helpers.R")

fn <- file.path(folder, "pf_sensitivity_mdl2.rds")
n_particles <- c(5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6)
sens_results2 <- pf_sensitivity(n_particles = n_particles, n_cores = 7, 
                                seed = 949549845, pomp_mdl = pomp_mdl, fn = fn,
                                n_iter = 14)