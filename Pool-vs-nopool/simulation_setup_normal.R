# Simulation setup for comparison of pooled reference null vs individual test

# Total number of tests
m <- 5000
# number of samples in each group
n <- 5
# Proportion of null
pi_0 <- 0.9

# Index of true null
true_null_idx <- 1:(m*pi_0)
non_null_idx <- (1:m)[-true_null_idx]

# Mean of each test in control group and treatment group
mu_control <- numeric(m)
mu_treat <- numeric(m)

mu_control <- rnorm(m, 0, 1)
mu_treat[true_null_idx] <- mu_control[true_null_idx]
mu_treat[non_null_idx] <- (abs(mu_control[non_null_idx]) + runif(length(non_null_idx), 0.5, 4))*sign(mu_control[non_null_idx])

# Standard deviation of each test in both groups
sigma_control <- rep(1, m)
sigma_treat <- sigma_control
sigma_treat[non_null_idx] <- rexp(length(non_null_idx), rate = 1)

# Generate data

control <- replicate(n, rnorm(m, mu_control, sigma_control))
treat <- replicate(n, rnorm(m, mu_treat, sigma_treat))

dat <- cbind(control, treat)
