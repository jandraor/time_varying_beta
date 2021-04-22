estimate_r <- function(beta_t, omega = 0.73, gamma = 1/2.9, eta = 1 / 2.1, 
                       mu = 0.5, kappa = 1 / 5) {
  
  # Symptomatic contribution to R
  symp_cont  <- omega * beta_t / gamma
  
  # Pre-symptomatic contribution to R
  psymp_cont <- beta_t / eta
  
  # Asymptomatic contribution to R
  asymp_cont <- (1 - omega) * beta_t * mu / kappa
  
  symp_cont + psymp_cont + asymp_cont
}

estimate_r_via_matrix <- function(sigma = 1 /3,omega = 0.73, gamma = 1/2.9, 
                                  eta = 1 / 2.1, mu = 0.5, kappa = 1 / 5) {
  
  matrix_F <- matrix(c(0, beta_t, beta_t, mu * beta_t,
                rep(0, 12)), nrow = 4, byrow = TRUE)
  
  matrix_V <- matrix(c(sigma, 0, 0, 0,
                       -sigma, eta, 0, 0,
                       0, -omega * eta, gamma, 0,
                       0, -(1 - omega) * eta, 0, kappa),
                       nrow = 4, byrow = TRUE)
  
  inverse_V <- solve(matrix_V)
  
  FV_inv <- matrix_F %*% inverse_V
  
  eigensystem <- eigen(FV_inv)
  
  max(eigensystem$values)
}