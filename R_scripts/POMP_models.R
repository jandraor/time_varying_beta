get_params <- function(mdl_id) {
  
  zeta_guess  <- 3 / 5
  B_guess     <- zeta_guess
  P_guess     <- 1
  
  N_val       <- 4999970
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
  
  if(mdl_id == "1") {
    unknown_pars <- c(B_0 = B_guess, P_0 = P_guess)
  }
  
  if(mdl_id == "2") {
    tau_guess    <- 1
    unknown_pars <- c(zeta = zeta_guess, P_0 = P_guess, tau = tau_guess)
  }
  
  params       <- c(fixed_pars, unknown_pars)
  
  list(fixed = fixed_pars,
       unknown = unknown_pars,
       all = params)
}

get_POMP_model <- function(obs_df, params) {
  
  Csnippet("
    y = rpois(C);  
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y,C,give_log) + 1e-5;
  ") -> dmeas
  
  Csnippet("
    S = 4999970 - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    B = B_0;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = B * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    B+= alpha*B*dW;
    C+= (rho*eta*P)*dt;
  ") -> SEI3R_GBM_step
  
  par_names <- names(params)
  
  obs_df %>%  
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step, delta.t = 0.01),
      statenames = c("S", "E", "P","I", "R", "A", "C", "B"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("B_0", "P_0")),
      cdir = ".", 
      cfile= "SEI3R_GBM"
    ) -> SEI3R_GBM
  
}

pomp_SEI3R_GBM2 <- function(obs_df, params) {
  Csnippet("
    y1  = rpois(C);  
    y2 = rnorm(B, tau);
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y1,C,give_log) + dnorm(y2, B, tau, give_log);
  ") -> dmeas
  
  Csnippet("
    S = 4999970 - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    B = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * B * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    B+= alpha*B*dW;
    C+= (rho*eta*P)*dt;
  ") -> SEI3R_GBM_step2
  

  par_names <- names(params)
  
  obs_df %>% 
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step2, delta.t = 0.01),
      statenames = c("S", "E", "P","I", "R", "A", "C", "B"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("zeta", "P_0", "tau")),
      obsnames = c("y1", "y2"),
      cdir = ".", 
      cfile= "SEI3R_GBM2"
    ) -> SEI3R_GBM2
}



