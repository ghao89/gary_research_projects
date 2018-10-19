# Demo code for examining performance of Liang's method

library(fdrDiscreteNull)
library(ggplot2)

# Load the simulated data
# dat: contains the simulated data used to perform multiple tests
#      the number of rows is the total number of tests
#      the 1st and 3rd columns contain the observed count and total number of trials of one binomial distribution
#      the 2nd and 4th columns contain the observed count and total number of trials of the other binomial distribtion
# non_null_idx: indices indicating non-null (alternative) hypotheses
# true_null_idx: indices indicating true null hypotheses

load(file = "dat.Rdata")
source("Helper/get_power.R")
source("Helper/get_fdr.R")

# Number of tests
m <- nrow(dat)

# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.04)

# Prepare the result storage
liang_power <- numeric(length(alpha))
liang_fdr <- numeric(length(alpha))

##########################################################
###################      Liang       #####################
########################################################## 

# Generate all the permutations
perms <- combn(seq(1:ncol(n_trial)), ncol(n_trial)/2)
perms <- perms[, 1:(ncol(perms)/2)]

# Storage of test statistic for each permutation
stat_liang <- matrix(0, nrow = m, ncol = ncol(perms))
for (j in 1:ncol(perms)) {
  b_1 <- apply(binom[, perms[,j]], 1, sum)
  b_2 <- apply(binom[, -perms[,j]], 1, sum)
  n_1 <- apply(n_trial[, perms[,j]], 1, sum)
  n_2 <- apply(n_trial[, -perms[,j]], 1, sum)
  stat_liang[, j] <- abs(b_1/n_1 - b_2/n_2)
}

# Calculate p-values based on permutations
rank <- apply(stat_liang, 1, FUN = function(x) rank(x, ties.method = "first")[1])
p_val_liang <- (ncol(perms) + 1 - rank)/(ncol(perms))

# Calculate the support of p-values 
p_cand <- seq(1/ncol(perms), 1, by = 1/ncol(perms))

# Calculate Liang's estimator of pi0
c <- c(0, p_cand[which(p_cand < 0.5)], 0.5)
pi_0_hat <- unlist(lapply(c, FUN = function(x, p_val) sum(p_val > x)/(1 - x)/length(p_val), p_val = p_val_liang))
idx <- which(pi_0_hat[-1] - pi_0_hat[-length(pi_0_hat)] > 0)
pi0_liang <- min(ifelse(length(idx) == 0, pi_0_hat[length(pi_0_hat)], pi_0_hat[min(idx) + 1]), 1)

# Storey's q-value method with Liang's pi0 estimator
# Use the GeneralizedFDREstimators function in fdrDiscreteNull package to get result
res <- GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.05)
q_liang <- pi0_liang*m*res$pvalues/unlist(lapply(res$pvalues, FUN = function(x, pval) max(sum(pval < x),1), pval = res$pvalues))
liang_detection <- lapply(alpha, FUN = function(x) which(q_liang <= x))

# Prepare results for plotting
liang_power <- unlist(lapply(liang_detection, FUN = function(x) get_power(x, non_null_idx)))
liang_fdr <- unlist(lapply(liang_detection, FUN = function(x) get_fdr(x, true_null_idx)))

summ <- data.frame(matrix(0, ncol = 4, nrow = length(alpha)))
summ[, 1] <- alpha
summ[, 2] <- liang_power
summ[, 3] <- liang_fdr
summ[, 4] <- rep("Liang", length(alpha))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot

