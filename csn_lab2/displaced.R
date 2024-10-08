require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

calculate_zeta_right_truncated <- function(x) {
  #basic analysis
  N <- length(x)
  M <- sum(x)
  M_prim <- sum(log(x))

  #function calculating minus log likelihood
  minus_log_likelihood <- function(gamma, k_max) {
    gamma * M_prim + N * log(zeta(gamma, shift=k_max))
  }
  
  #start parameter value
  gamma_0 = 2
  k_max_0 = max(x)
  
  #calculating optimal parameter within bounds
  mle_calc <- mle(minus_log_likelihood,
                  start = list(gamma = gamma_0, k_max = k_max_0),
                  method = "L-BFGS-B",
                  lower = c(1.000001, 1))
  
  #calculating result values
  gamma <- attributes(summary(mle_calc))$coef[1]
  L <- minus_log_likelihood(gamma)
  AIC <- get_AIC(2*L, 2, N)
  return(list("ll"=-L, "param"=gamma, "param_0"=gamma_0, "AIC"=AIC))
}

calculate_zeta_gamma_2 <- function(x) {
  #basic analysis
  N <- length(x)
  M <- sum(x)
  M_prim <- sum(log(x))

  #calculating result values (parameters are constant)
  L <- 2 * M_prim + N* log(zeta(2))
  AIC <- get_AIC(2*L, 0, N)
  return(list("ll"=-L, "param"=2, "param_0"=2, "AIC"=AIC))
}


calculate_zeta <- function(x) {
  #basic analysis
  N <- length(x)
  M <- sum(x)
  M_prim <- sum(log(x))
  
  #function calculating minus log likelihood
  minus_log_likelihood <- function(gamma) {
    gamma * M_prim + N * log(zeta(gamma))
  }
  
  #start parameter value
  gamma_0 = 2
  
  #calculating optimal parameter within bounds
  mle_calc <- mle(minus_log_likelihood,
                  start = list(gamma = gamma_0),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  
  #calculating result values
  gamma <- attributes(summary(mle_calc))$coef[1]
  L <- minus_log_likelihood(gamma)
  AIC <- get_AIC(2*L, 1, N)
  return(list("ll"=-L, "param"=gamma, "param_0"=gamma_0, "AIC"=AIC))
}

calculate_poisson <- function(x) {
  #basic analysis
  N <- length(x)
  M <- sum(x)
  
  c_inner_sum <- function(i) {sum(log(2:x[i]))}
  C <- sum(sapply(1:N, FUN=c_inner_sum))
  
  #function calculating minus log likelihood
  minus_log_likelihood <- function(lambda) {
    -(M * log(lambda) - N * (lambda + log(1 - exp(-lambda))) - C)
  }
  
  #start parameter value
  lambda_0 = M/N
  
  #calculating optimal parameter within bounds
  mle_calc <- mle(minus_log_likelihood,
                  start = list(lambda = lambda_0),
                  method = "L-BFGS-B",
                  lower = c(0.0000001))
  
  #calculating result values
  lambda <- attributes(summary(mle_calc))$coef[1]
  L <- minus_log_likelihood(lambda)
  AIC <- get_AIC(2*L, 1, N)
  return(list("ll"=-L, "param"=lambda, "param_0"=lambda_0, "AIC"=AIC))
}

calculate_geometric <- function(x) {
  #basic analysis
  N <- length(x)
  M <- sum(x) 

  #function calculating minus log likelihood
  minus_log_likelihood <- function(q) {
    -((M - N) * log(1-q) + N*log(q))
  }
  
  #start parameter value
  q_0 <- N/M
  
  #calculating optimal parameter within bounds
  mle_geom <- mle(minus_log_likelihood,
                  start = list(q = q_0),
                  method = "L-BFGS-B",
                  lower = c(0.0000001),
                  upper = c(0.9999999))
  
  #calculating result values
  q <- attributes(summary(mle_geom))$coef[1]
  L <- minus_log_likelihood(q)
  AIC <- get_AIC(2*L, 1, N)
  return(list("ll"=-L, "param"=q, "param_0"=q_0, "AIC"=AIC))
}
