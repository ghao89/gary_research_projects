# Simulation setup for multiple Fisher's Exact Test
# For details of the simulation study, check the book chapter

# Total number of tests
m <- 1000
# number of samples in each group
n <- 4
# Proportion of null
pi_0 <- 0.1

# Index of true null
true_null_idx <- sample(m, m*pi_0, replace = FALSE)
non_null_idx <- ifelse(length(true_null_idx) == 0, 1:m, (1:m)[-true_null_idx])

# Mean of each test in control group and treatment group
mu_control <- numeric(m)
mu_treat <- numeric(m)

mu_control <- rnorm(m, 0, 1)
mu_treat[true_null_idx] <- mu_control[true_null_idx]
mu_treat[non_null_idx] <- (abs(mu_control[non_null_idx]) + runif(length(non_null_idx), 5, 10))*sign(mu_control[non_null_idx])

# Standard deviation of each test in both groups
sigma_control <- 1
sigma_treat <- 1

# Generate data

control <- replicate(n, rnorm(m, mu_control, sigma_control))
treat <- replicate(n, rnorm(m, mu_treat, sigma_treat))

dat <- cbind(control, treat)
