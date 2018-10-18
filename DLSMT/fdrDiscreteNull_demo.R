# Demo code for examining performance of FDR control procedures including:
# Benjamini-Hochberg procedure (BH)
# Benjamini-Hochberg-Heyse procedure (BHH)
# adaptive Benjamini-Hochberg procedure (aBH)
# adaptive Benjamini-Hochberg-Heyse procedure (aBHH)
# Habiger's method
# as all methods are implemented in package "fdrDiscreteNull"

library(fdrDiscreteNull)
library(ggplot2)

# Load the simulated data
# dat: contains the simulated data used to perform multiple tests
#      the number of rows is the total number of tests
#      the 1st and 3rd columns contain the observed count and total number of trials of one binomial distribution
#      the 2nd and 4th columns contain the observed count and total number of trials of the other binomial distribtion
# non_null_idx: indices indicating non-null (alternative) hypotheses
# true_null_idx: indices indicating true null hypotheses

load(file = "DLSMT/dat.Rdata")
source("DLSMT/Helper/get_power.R")
source("DLSMT/Helper/get_fdr.R")

# Number of replications for Habiger's method
rep <- 5

# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.02)

# Prepare the result storage
BH_power <- numeric(length(alpha))
aBH_power <- numeric(length(alpha))
BHH_power <- numeric(length(alpha))
aBHH_power <- numeric(length(alpha))
habiger_power <- matrix(0, nrow = length(alpha), ncol = rep)

BH_fdr <- numeric(length(alpha))
aBH_fdr <- numeric(length(alpha))
BHH_fdr <- numeric(length(alpha))
aBHH_fdr <- numeric(length(alpha))
habiger_fdr <- matrix(0, nrow = length(alpha), ncol = rep)

methods <- c("BH", "BHH", "aBH", "aBHH", "Habiger")


# Use the GeneralizedFDREstimators function in fdrDiscreteNull package to get result
res_alpha <- lapply(alpha, FUN = function(x) GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = x))

# BH procedure
BH_detection <- lapply(res_alpha, FUN = function(x) x$BH$IndicesOfDiscoveries)
BH_power <- unlist(lapply(BH_detection, FUN = function(x) get_power(x, non_null_idx)))
BH_fdr <- unlist(lapply(BH_detection, FUN = function(x) get_fdr(x, true_null_idx)))

# BHH method
BHH_detection <- lapply(res_alpha, FUN = function(x) x$BHH$IndicesOfDiscoveries)
BHH_power <- unlist(lapply(BHH_detection, FUN = function(x) get_power(x, non_null_idx)))
BHH_fdr <- unlist(lapply(BHH_detection, FUN = function(x) get_fdr(x, true_null_idx)))

# aBH method
aBH_detection <- lapply(res_alpha, FUN = function(x) x$aBH$IndicesOfDiscoveries)
aBH_power <- unlist(lapply(aBH_detection, FUN = function(x)  get_power(x, non_null_idx)))
aBH_fdr <- unlist(lapply(aBH_detection, FUN = function(x) get_fdr(x, true_null_idx)))

# aBHH method
aBHH_detection <- lapply(res_alpha, FUN = function(x) x$aBHH$IndicesOfDiscoveries)
aBHH_power <- unlist(lapply(aBHH_detection, FUN = function(x) get_power(x, non_null_idx)))
aBHH_fdr <- unlist(lapply(aBHH_detection, FUN = function(x) get_fdr(x, true_null_idx)))


for (i in 1:rep) { 
  res_alpha <- lapply(alpha, FUN = function(x) GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = x))
  # Habiger method
  habiger_detection <- lapply(res_alpha, FUN = function(x) x$SARP$IndicesOfDiscoveries)
  habiger_power[, i] <- unlist(lapply(habiger_detection, FUN = function(x) get_power(x, non_null_idx)))
  habiger_fdr[, i] <- unlist(lapply(habiger_detection, FUN = function(x) get_fdr(x, true_null_idx)))
  
  # Monitor the process
  print(paste0("The ", i, "th round is completed!"))
}

habiger_power_mean <- apply(habiger_power, 1, mean)
habiger_fdr_mean <- apply(habiger_fdr, 1, mean)
# habiger_power_sd <- apply(habiger_power, 1, sd)
# habiger_fdr_sd <- apply(habiger_fdr, 1, sd)

# Prepare results for plotting
summ <- data.frame(matrix(0, ncol = 4, nrow = length(alpha)*length(methods)))
summ[, 1] <- rep(alpha, length(methods))
summ[, 2] <- c(BH_power, BHH_power, aBH_power, aBHH_power, habiger_power_mean)
summ[, 3] <- c(BH_fdr, BHH_fdr, aBH_fdr, aBHH_fdr, habiger_fdr_mean)
summ[, 4] <- rep(methods, each = length(alpha))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot




