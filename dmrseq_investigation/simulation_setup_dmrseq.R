# Simulation setup for investigating dmrseq method

# Cut m tests (sites) evenly into r regions so that each region has same number of tests (sites) (prototype)
# Each region follows a Beta-binomial setup: 
# 1. the mean methylation level for each site is drawn from a Beta distribution; 
# 2. the number of trials for each site each sample is drawn uniformly from [5, 100]
# 3. the mean of Beta distribution for the ith region in control group is calculated from a Sin-function:
#    mu_i = sin(pi*(i/r))
# 4. the mean of Beta distribution for the ith region in treatment group:
#    i.  equals that of control group if that region is null
#    ii. (mu_i + delta) %% 1, if that region is non-null. delta is uniformly drawn from [0.1, 0.5]
# The proportion of null regions is pre-set by pi_0

# Total number of tests
m <- 5000

# Number of regions
r <- 100

# number of samples in each group
n <- 2

# Proportion of null
pi_0 <- 0.6

# Beta mean
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

