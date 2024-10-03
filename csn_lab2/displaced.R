get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

displaced_poisson <- function(x) {
  N <- length(x)
  M <- sum(x)

  c_inner_sum <- function(i) {sum(log(2:x[i]))}
  C <- sum(sapply(1:N, FUN=c_inner_sum))
  
  lambda <- M/N
  
  L <- (M * log(lambda) - N * (lambda + log(1 - exp(-lambda))) - C)
  K <- 1
  
  return(list("L"=L, "param"=lambda, "AIC"=get_AIC(-2*L, 1, N)))
}

displaced_geometric <- function(x) {
  N <- length(x)
  M <- sum(x)
  q <- N/M
  
  L <- ((M - N) * log((1-q), 2)) + (N*log(q, 2))
  K <- 1
  
  return(list("L"=L, "param"=q, "AIC"=get_AIC(-2*L, 1, N)))
}
