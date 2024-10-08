# Load Required Libraries
library(stats4)  # For MLE
library(VGAM)    # For zeta function

# Load Degree Sequence
degree_sequence <- read.table("./data/English_in-degree_sequence.txt", header = FALSE)

# Check Input Data
if (any(degree_sequence$V1 < 0)) stop("Degree sequence contains negative values.")
if (any(is.na(degree_sequence$V1))) stop("Degree sequence contains NA values.")

# Summary Statistics Calculation
N <- nrow(degree_sequence)
max_degree <- max(degree_sequence$V1)
sum_degrees <- sum(degree_sequence$V1)
mean_degree <- sum_degrees / N
inverse_mean_degree <- N / sum_degrees

# Create Summary Table for Properties
summary_table <- data.frame(Language = "English", 
                            N = N, 
                            Max_Degree = max_degree, 
                            Mean_Degree = mean_degree, 
                            Inverse_Mean_Degree = inverse_mean_degree)

# Print Summary Table
print("Summary Table of Properties:")
print(summary_table)

# Log-Likelihood Functions for Each Model

# Displaced Poisson
minus_log_likelihood_poisson <- function(lambda) {
  if (lambda <= 0) return(Inf)  # Ensure lambda is positive
  M <- sum(degree_sequence$V1)
  -N * (lambda + log(1 - exp(-lambda))) - sum(log(factorial(degree_sequence$V1)))
}

# Displaced Geometric
minus_log_likelihood_geometric <- function(q) {
  if (q <= 0 || q >= 1) return(Inf)  # Ensure q is in (0,1)
  M <- sum(degree_sequence$V1)
  -((M - N) * log(1 - q) + N * log(q))
}

# Zeta Distribution
minus_log_likelihood_zeta <- function(gamma) {
  if (gamma <= 1) return(Inf)  # Ensure gamma > 1
  -length(degree_sequence$V1) * log(zeta(gamma)) - gamma * sum(log(degree_sequence$V1))
}

# Right-Truncated Zeta Distribution
minus_log_likelihood_right_truncated_zeta <- function(gamma, kmax) {
  if (gamma <= 1 || kmax < 1) return(Inf)  # Ensure gamma > 1 and kmax >= 1
  M0 <- sum(log(degree_sequence$V1))
  C <- sum(sapply(degree_sequence$V1, function(k) sum(log(1:k))))
  -(-gamma * M0 - N * log(zeta(gamma, kmax)) + N * log(gamma) - C)
}

# Altmann Function
minus_log_likelihood_altmann <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) return(Inf)  # Ensure alpha, beta > 0
  c <- 1 / sum(degree_sequence$V1^(-alpha) * exp(-beta * degree_sequence$V1))
  L <- sum(log(c * degree_sequence$V1^(-alpha) * exp(-beta * degree_sequence$V1)))
  return(-L)
}

# Parameter Estimation and AIC Calculation

# 1. Displaced Poisson
mle_poisson <- tryCatch(
  mle(minus_log_likelihood_poisson, start = list(lambda = mean_degree * 1.5), method = "L-BFGS-B"),
  error = function(e) { cat("Error in Displaced Poisson MLE:", e$message, "\n"); return(NULL) }
)

# 2. Displaced Geometric
initial_q <- min(max(1 - (mean_degree / N), 0.01), 0.99)  # Ensure q is between 0 and 1
mle_geometric <- tryCatch(
  mle(minus_log_likelihood_geometric, start = list(q = initial_q), method = "L-BFGS-B"),
  error = function(e) { cat("Error in Displaced Geometric MLE:", e$message, "\n"); return(NULL) }
)

# 3. Zeta Distribution
mle_zeta <- tryCatch(
  mle(minus_log_likelihood_zeta, start = list(gamma = 2.5), method = "L-BFGS-B", lower = c(1.0000001)),
  error = function(e) { cat("Error in Zeta MLE:", e$message, "\n"); return(NULL) }
)

# 4. Right-Truncated Zeta Distribution
mle_right_truncated_zeta <- tryCatch(
  mle(minus_log_likelihood_right_truncated_zeta, 
      start = list(gamma = 2.5, kmax = max_degree), 
      method = "L-BFGS-B", 
      lower = c(1.0000001, 1)),
  error = function(e) { cat("Error in Right-Truncated Zeta MLE:", e$message, "\n"); return(NULL) }
)

# 5. Right_truncated Zeta
mle_altmann <- tryCatch(
  mle(minus_log_likelihood_altmann, 
      start = list(alpha = 1, beta = 1), 
      method = "L-BFGS-B", 
      lower = c(0.0001, 0.0001)),
  error = function(e) { cat("Error in Altmann MLE:", e$message, "\n"); return(NULL) }
)

# AIC Calculations
get_AIC <- function(mle_obj, K) {
  if (is.null(mle_obj)) return(NA)
  return(-2 * attributes(summary(mle_obj))$m2logL + 2 * K)
}

# AIC Values
AIC_poisson <- get_AIC(mle_poisson, 1)
AIC_geometric <- get_AIC(mle_geometric, 1)
AIC_zeta <- get_AIC(mle_zeta, 1)
AIC_right_truncated_zeta <- get_AIC(mle_right_truncated_zeta, 2)
AIC_altmann <- get_AIC(mle_altmann, 2)

# Summary Tables

# Summary Table for AIC
model_names <- c("Displaced Poisson", "Displaced Geometric", "Zeta", "Right-Truncated Zeta", "Altmann")
AIC_values <- c(AIC_poisson, AIC_geometric, AIC_zeta, AIC_right_truncated_zeta, AIC_altmann)

summary_AIC_table <- data.frame(Model = model_names, AIC = AIC_values)
print("Summary Table of AIC Values:")
print(summary_AIC_table)

# Summary Table for Parameters
parameters_table <- data.frame(
  Model = model_names,
  Lambda = c(ifelse(!is.null(mle_poisson), coef(mle_poisson), NA)),
  Q = c(ifelse(!is.null(mle_geometric), coef(mle_geometric), NA)),
  Gamma = c(ifelse(!is.null(mle_zeta), coef(mle_zeta), NA), 
            ifelse(!is.null(mle_right_truncated_zeta), coef(mle_right_truncated_zeta)[1], NA), 
            ifelse(!is.null(mle_altmann), coef(mle_altmann)[1], NA))),
Kmax = c(NA, NA, NA, 
         ifelse(!is.null(mle_right_truncated_zeta), coef(mle_right_truncated_zeta)[2], NA), NA),
Alpha = c(NA, NA, NA, NA, 
          ifelse(!is.null(mle_altmann), coef(mle_altmann)[1], NA)),
Beta = c(NA, NA, NA, NA, 
         ifelse(!is.null(mle_altmann), coef(mle_altmann)[2], NA))
)

print("Summary Table of Estimated Parameters:")
print(parameters_table)

# Identify the Best Model
best_model <- summary_AIC_table[which.min(AIC_values),]
print(paste("Best Model is:", best_model$Model))

