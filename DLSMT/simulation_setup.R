# Simulation setup for multiple Fisher's Exact Test
# For details of the simulation study, check the book chapter

# Total number of tests
m <- 1000
# number of samples in each group
n <- 2
# Proportion of null
pi_0 <- 0.6
# Poisson mean
mu_pois <- 20
# Index of true null
true_null_idx <- sample(m, m*pi_0, replace = FALSE)
non_null_idx <- (1:m)[-true_null_idx]

# Binomial mean
p_1 <- rep(0, m)
p_2 <- rep(0, m)

p_1[true_null_idx] <- runif(length(true_null_idx))
p_2[true_null_idx] <- p_1[true_null_idx]
p_1[non_null_idx] <- runif(length(non_null_idx), 0, 0.5)
p_2[non_null_idx] <- p_1[non_null_idx] + runif(length(non_null_idx), 0.2, 0.5)

# Number of trials
n_trial_1 <- replicate(n, rpois(m, mu_pois))
n_trial_2 <- replicate(n, rpois(m, mu_pois))
n_trial <- cbind(n_trial_1, n_trial_2)

# Binomial draws
bin_1 <- apply(n_trial_1, 2, FUN = function(x, m, p) rbinom(m, x, p), m = m, p = p_1)
bin_2 <- apply(n_trial_2, 2, FUN = function(x, m, p) rbinom(m, x, p), m = m, p = p_2)
binom <- cbind(bin_1, bin_2)

# The constructed data can be used directly for fdrDiscreteNull package "GeneralizedFDREstimators"
# function
dat <- cbind(apply(bin_1, 1, sum), apply(bin_2, 1, sum), apply(n_trial_1, 1, sum), apply(n_trial_2, 1, sum))

