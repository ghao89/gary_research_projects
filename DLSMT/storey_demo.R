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
source("Helper/Generalized_FDR_Estimators.R")
source("Helper/MultipleTesting_HeteroDist_ModFuncs.R")

# Number of tests
m <- nrow(dat)

# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.04)

# Prepare the result storage
storey_power <- numeric(length(alpha))
storey_fdr <- numeric(length(alpha))
storey_rnd_power <- numeric(length(alpha))
storey_rnd_fdr <- numeric(length(alpha))

##########################################################
###################      Storey      #####################
##########################################################  


cellcountsmarginals <- getcellcountsandmarginals_DE(dat)
simallcellcounts <- cellcountsmarginals[[1]]

pvsupp <- lapply(simallcellcounts, FUN = function(x) pvalFETSupport(x))
rawpvalues <- unlist(lapply(pvsupp, FUN = function(x) x$rawpvalues))
rndpvalues <- unlist(lapply(pvsupp, FUN = function(x) x$rndpvalues))

# Based on raw p values
pi0_storey <- storeyPi0Est(0.5,rawpvalues)
storey_detection <- lapply(alpha, 
                           FUN = function(x, pval, pi_0)
                             StoreyFDREst(pval, x, pi_0)[[1]],
                           pval = rawpvalues,
                           pi_0 = pi0_storey)

# Based on randomized p values
pi0_rnd <- storeyPi0Est(0.5, rndpvalues)
storey_rnd_detection <- lapply(alpha, 
                               FUN = function(x, pval, pi_0)
                                 StoreyFDREst(pval, x, pi_0)[[1]],
                               pval = rndpvalues,
                               pi_0 = pi0_rnd)

# Prepare results for plotting
storey_power <- unlist(lapply(storey_detection, FUN = function(x) get_power(x, non_null_idx)))
storey_fdr <- unlist(lapply(storey_detection, FUN = function(x) get_fdr(x, true_null_idx)))
storey_rnd_power <- unlist(lapply(storey_rnd_detection, FUN = function(x) get_power(x, non_null_idx)))
storey_rnd_fdr <- unlist(lapply(storey_rnd_detection, FUN = function(x) get_fdr(x, true_null_idx)))

summ <- data.frame(matrix(0, ncol = 4, nrow = 2*length(alpha)))
summ[, 1] <- c(alpha,alpha)
summ[, 2] <- c(storey_power, storey_rnd_power)
summ[, 3] <- c(storey_fdr, storey_rnd_fdr)
summ[, 4] <- c(rep("Storey-Raw", length(alpha)), rep("Storey-Randomized", length(alpha)))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot
