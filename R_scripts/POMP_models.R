get_fixed_params <- function() {
  sigma_val   <- 1 / 3   # Incubation period
  eta_val     <- 1 / 2.1 # Duration of preclinical infectiousness
  gamma_val   <- 1 / 2.9 # Duration of clinical infectiousness
  kappa_val   <- 1 / 5 # Duration of subclinical infectiousness
  N_val       <- 4937796
  mu_val      <- 1 / 2 # Relative infectiousness
  omega_val   <- 0.73 # Clinical fraction
  
  c(N = N_val, mu = mu_val, omega = omega_val, gamma = gamma_val,
    sigma = sigma_val, eta = eta_val, kappa = kappa_val)
}

get_params <- function(mdl_id) {
  
  zeta_guess    <- 3 / 5
  tau_guess     <- 1
  B_guess       <- zeta_guess
  P_guess       <- 1
  alpha_guess   <- 0.1
  nu_guess      <- 0.04
  upsilon_guess <- 0.1
  iota_guess    <- 1
  
  if(mdl_id == "GBM_1") {
    unknown_pars <- c(zeta = zeta_guess, P_0 = P_guess, tau = tau_guess,
                      alpha = alpha_guess)
  }
  
  if(mdl_id == "GBM_2") {
    unknown_pars <- c(zeta = zeta_guess, P_0 = P_guess, tau = tau_guess,
                      alpha = alpha_guess)
  }
  
  if(mdl_id == "GBM_3") {
    unknown_pars <- c(zeta = zeta_guess, P_0 = P_guess, tau = tau_guess,
                      alpha = alpha_guess, iota = iota_guess)
  }
  
  if(mdl_id == "CIR") {
    unknown_pars  <- c(zeta = zeta_guess, P_0 = P_guess, tau = tau_guess,
                       nu = nu_guess, upsilon = upsilon_guess,
                       alpha = alpha_guess)
  }
  
  fixed_pars   <- get_fixed_params()
  
  params       <- c(fixed_pars, unknown_pars)
  
  list(fixed = fixed_pars,
       unknown = unknown_pars,
       all = params)
}

get_POMP_model <- function(obs_df, params, dt = 0.01) {
  
  Csnippet("
    y1  = rpois(C);  
    y2 = rnorm(Z, tau); 
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y1,C,give_log) + dnorm(y2, Z, tau, give_log);
  ") -> dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * Z * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= alpha*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_GBM_step
  
  par_names <- names(params)
  
  obs_df %>%  
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("zeta", "P_0", "tau", "alpha")),
      cdir = ".", 
      cfile= "SEI3R_GBM"
    ) -> SEI3R_GBM
  
}

pomp_SEI3R_GBM2 <- function(obs_df, params, dt = 0.01) {
  
  Csnippet("
    y1  = rpois(C);  
    y2 = rnorm(Z, tau);
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y1,C,give_log) + dnorm(y2, Z, tau, give_log);
  ") -> dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * Z * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= alpha*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_GBM_step2
  

  par_names <- names(params)
  
  obs_df %>% 
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_step2, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("zeta", "P_0", "tau", "alpha")),
      obsnames = c("y1", "y2"),
      cdir = ".", 
      cfile= "SEI3R_GBM2"
    ) -> SEI3R_GBM2
}

pomp_SEI3R_CIR <- function(obs_df, pars, dt = 0.01) {
  Csnippet("
    y1  = rpois(C);  
    y2 = rnorm(Z, tau);
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y1,C,give_log) + dnorm(y2, Z, tau, give_log);
  ") -> dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW     = rnorm(0,sqrt(dt));
    double lambda = zeta * Z * (I + P + mu * A) / N;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= (nu * (upsilon - Z)) * dt + sqrt(alpha)*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_CIR_step
  
  par_names <- names(pars)
  
  obs_df %>% 
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_CIR_step, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("zeta", "P_0", "tau", "alpha"),
                                 logit = c("upsilon", "nu")),
      obsnames = c("y1", "y2"),
      cdir = ".", 
      cfile= "SEI3R_CIR"
    ) -> SEI3R_CIR
}

pomp_SEI3R_GBM_immi <- function(obs_df, params, dt = 0.01) {
  
  Csnippet("
    y1  = rpois(C);  
    y2  = rnorm(Z, tau);
  ") -> rmeas
  
  Csnippet("
    lik = dpois(y1,C,give_log) + dnorm(y2, Z, tau, give_log);
  ") -> dmeas
  
  Csnippet("
    S = N - P_0;
    E = 0;
    P = P_0;
    I = 0;
    A = 0;
    R = 0;
    Z = 1;
    C = 0;
  ") -> rinit
  
  Csnippet("
    double dW         = rnorm(0,sqrt(dt));
    double beta       = zeta * Z;
    double lambda_end = beta * (I + P + mu * A) / N;
    double lambda_exo = (beta * iota) / N;
    double lambda     = lambda_end + lambda_exo;
    S-= (S * lambda)*dt;
    E+= (S * lambda - sigma*E)*dt;
    P+= (omega * sigma * E - eta * P) * dt;
    I+= (eta * P - gamma*I)*dt;
    R+= (gamma*I + kappa * A)*dt;
    A+= ((1-omega) * sigma * E - kappa * A) * dt;
    Z+= alpha*Z*dW;
    C+= (eta*P)*dt;
  ") -> SEI3R_GBM_immi_step
  
  par_names <- names(params)
  
  obs_df %>% 
    pomp(
      times = "time", t0 = 0,
      rinit = rinit,
      rprocess = pomp::euler(SEI3R_GBM_immi_step, delta.t = dt),
      statenames = c("S", "E", "P","I", "R", "A", "C", "Z"),
      paramnames = par_names,
      params = params,
      accumvars = "C",
      rmeasure = rmeas,
      dmeasure = dmeas,
      partrans = parameter_trans(log   = c("zeta", "P_0", "tau", "alpha", 
                                           "iota")),
      obsnames = c("y1", "y2"),
      cdir = ".", 
      cfile= "SEI3R_GBM_immi"
    ) -> SEI3R_GBM_immi
}