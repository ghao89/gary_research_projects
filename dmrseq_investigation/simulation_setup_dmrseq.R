# Simulation setup for investigating dmrseq method

###############################################################
##############        Simulation guidence       ###############
###############################################################

# Cut m tests (sites) evenly into r regions so that each region has same number of tests (sites) (prototype)
# Each region follows a Beta-binomial setup: 
# 1. the mean of Beta distribution for the ith region in control group is calculated from a Sin-function:
#    mu_i = sin(0.99*pi*(i/r)), use 0.99 to avoid mu = 1
# 2. the mean of Beta distribution for the ith region in treatment group:
#    i.  equals that of control group if that region is null
#    ii. (mu_i + delta) %% 1, if that region is non-null. delta is uniformly drawn from [0.05, 0.5]
# 3. the mean of Beta distribution for each site in the i-th region is drawn uniformly from [mu_i - 0.05, mu_i + 0.05];
# 4. the mean methylation level for each site across different samples (within group) is drawn from a Beta distribution; 
# 4. the number of trials for each site each sample is drawn uniformly from [5, 100]
# The proportion of null regions is pre-set by pi_0
# The number of samples in each group is n


###############################################################
##############        Simulation Setup        #################
###############################################################

# Total number of tests
m <- 50000

# Number of regions
r <- 500

# number of samples in each group
n <- 3

# Proportion of null
pi_0 <- 0.9

# Index of true null and true non-null
true_null_idx <- sample(r, r*pi_0, replace = FALSE)
non_null_idx <- (1:r)[-true_null_idx]

# Mean of Beta means on each null region (periodic seq(0.1, 1, by = 0.1))
mu_control <- (0.2*(1:r - 1) + 0.1) %% 1

# Mean of Beta means on each non-null region
mu_treat <- numeric(r)
mu_treat[true_null_idx] <- mu_control[true_null_idx]
mu_treat[non_null_idx] <- (mu_control[non_null_idx] + runif(length(non_null_idx), 0.05, 0.5)) %% 1

# Beta mean for each site in null region
beta_mean_control <- as.vector(t(replicate(m/r, runif(r, mu_control - 0.05, mu_control + 0.05))))
# Beta mean for each site in non-null region
beta_mean_treat <- as.vector(t(replicate(m/r, runif(r, ifelse(mu_treat - 0.05 > 0, mu_treat - 0.05, 0.005), ifelse(mu_treat + 0.05 >= 1, 0.995, mu_treat + 0.05)))))

# a (shape1) and b (shape2) parameterization for each Beta distribution
nv <- 100
a_control <- beta_mean_control*nv
b_control <- (1 - beta_mean_control)*nv
a_treat <- beta_mean_treat*nv
b_treat <- (1 - beta_mean_treat)*nv

# Generate the data
# Coverage in each group
coverage_control <- matrix(0, nrow = m, ncol = n)
coverage_treat <- matrix(0, nrow = m, ncol = n)

# Methylation in each group
methy_control <- matrix(0, nrow = m, ncol = n)
methy_treat <- matrix(0, nrow = m, ncol = n)

# Methylation level in each group
methy_level_control <- matrix(0, nrow = m, ncol = n)
methy_level_treat <- matrix(0, nrow = m, ncol = n)

for (i in 1:n) {
  coverage_control[, i] <- sample(5:100, m, replace = T)
  coverage_treat[, i] <- sample(5:100, m, replace = T)
  methy_level_control[, i] <- rbeta(m, a_control, b_control)
  methy_level_treat[, i] <- rbeta(m, a_treat, b_treat)
  methy_control[, i] <- rbinom(m, coverage_control[, i], methy_level_control[, i])
  methy_treat[, i] <- rbinom(m, coverage_treat[, i], methy_level_treat[, i])
}

# Construct BSseq data
sampleNames <- c(paste0("Control_", 1:n), paste0("Treat_", 1:n))
bs_dat <- BSseq(chr = rep("chr24", m),
                pos = 1:m,
                M = cbind(methy_control, methy_treat),
                Cov = cbind(coverage_control, coverage_treat),
                sampleNames = sampleNames)

save(bs_dat, file = "~/Dropbox/Private/Git_Projects/DLSMT/dmrseq_investigation/bs_dat.Rdata")

non_null_range <- IRanges(start = non_null_idx*100 - 99, end = non_null_idx*100)
