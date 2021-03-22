estimate_r <- function(beta_t, omega = 0.73, gamma = 1/2.9, eta = 1 / 2.1, 
                       mu = 0.5, kappa = 1 / 5) {
  
  # Symptomatic contribution to R
  symp_cont  <- omega * beta_t / gamma
  
  # Pre-symptomatic contribution to R
  psymp_cont <- omega * beta_t / eta
  
  # Asymptomatic contribution to R
  asymp_cont <- (1 - omega) * beta_t * mu / kappa
  
  symp_cont + psymp_cont + asymp_cont
}