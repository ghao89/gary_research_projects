# Demo code for examining performance of Storey's q-value method

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
storey_power <- numeric(length(alpha))
storey_fdr <- numeric(length(alpha))

##########################################################
###################      Storey      #####################
##########################################################  

# Use the GeneralizedFDREstimators function in fdrDiscreteNull package to get result
res <- GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.05)

lambda_storey <- 0.5
pi0_storey <- min((1 + sum(res$pvalues > lambda_storey))/((1 - lambda_storey)*m), 1)
q_storey <- pi0_storey*m*res$pvalues/unlist(lapply(res$pvalues, FUN = function(x, pval) max(sum(pval < x),1), pval = res$pvalues))
storey_detection <- lapply(alpha, FUN = function(x) which(q_storey <= x))

# Prepare results for plotting
storey_power <- unlist(lapply(storey_detection, FUN = function(x) get_power(x, non_null_idx)))
storey_fdr <- unlist(lapply(storey_detection, FUN = function(x) get_fdr(x, true_null_idx)))

summ <- data.frame(matrix(0, ncol = 4, nrow = length(alpha)))
summ[, 1] <- alpha
summ[, 2] <- storey_power
summ[, 3] <- storey_fdr
summ[, 4] <- rep("Storey", length(alpha))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot
