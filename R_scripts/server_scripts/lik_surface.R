source("./R_scripts/irish_data.R")

irish_data <- get_irish_data()

N_val       <- 4999970

zeta_guess  <- 3 / 5
P_guess     <- 1

gamma_val   <- 1 / 3 # Clinical
sigma_val   <- 1 / 3
rho_val     <- 1
alpha_val   <- 0.1
mu_val      <- 1 / 2 # Relative infectiousness
omega_val   <- 0.7 # Clinical fraction
eta_val     <- 1 / 2.1 # Preclinical
kappa_val   <- 1 / 5 # Subclinical

fixed_pars <- c(N = N_val, mu = mu_val, omega = omega_val, gamma = gamma_val,
                sigma = sigma_val, alpha = alpha_val, eta = eta_val, 
                kappa = kappa_val, rho = rho_val)

unknown_pars <- c(zeta = zeta_guess, P_0 = P_guess)

params <- c(fixed_pars, unknown_pars)

source("./R_scripts/POMP_models.R")

obs_df    <- irish_data %>% select(-date)
SEI3R_GBM <- get_POMP_model(obs_df, params)

params_bounds <- list(list(name = "P_0", min = 1, max = 30),
                      list(name = "zeta", min = 0.1, max = 3))

combinations  <- combn(1:length(params_bounds), 2, simplify = FALSE)
#-------------------------------------------------------------------------------
# This script is equivalent to chunk 'lik_surface' in Ireland_GBM.Rmd

source("./R_scripts/likelihood_slice.R")
source("./R_scripts/helpers.R")

params_bounds <- list(list(name = "P_0", min = 1, max = 30),
                      list(name = "zeta", min = 0.1, max = 3))

combinations  <- combn(1:length(params_bounds), 2, simplify = FALSE)

folder <- "./Saved_objects/Irish_data/SEI3R_GBM/"

p_list <- likelihood_slice(combinations, folder)
#-------------------------------------------------------------------------------