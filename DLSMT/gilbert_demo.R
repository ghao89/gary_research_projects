# Demo code for examining performance of Gilbert method

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
source("Helper/find_Gilbert_threshold.R")

# Number of tests
m <- nrow(dat)

# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.04)

# Prepare the result storage
gilbert_power <- numeric(length(alpha))
gilbert_fdr <- numeric(length(alpha))

##########################################################
###################   Gilbert 2005   #####################
##########################################################

# Use the GeneralizedFDREstimators function in fdrDiscreteNull package to get result
res <- GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.05)
pval_supp_min <- unlist(lapply(res$pvalSupp, min))
# Find the p-value threshold from the subset of tests
thres_gilbert <- unlist(lapply(alpha, FUN = function(x, p_val, pval_supp_min) find_Gilbert_threshold(p_val, pval_supp_min, x), p_val = res$pvalues, pval_supp_min = pval_supp_min))
gilbert_detection <- lapply(thres_gilbert, FUN = function(x, p_val) which(p_val <= x), p_val = res$pvalues)

# Prepare results for plotting
gilbert_power <- unlist(lapply(gilbert_detection, FUN = function(x) get_power(x, non_null_idx)))
gilbert_fdr <- unlist(lapply(gilbert_detection, FUN = function(x) get_fdr(x, true_null_idx)))

summ <- data.frame(matrix(0, ncol = 4, nrow = length(alpha)))
summ[, 1] <- alpha
summ[, 2] <- gilbert_power
summ[, 3] <- gilbert_fdr
summ[, 4] <- rep("Gilbert", length(alpha))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot

