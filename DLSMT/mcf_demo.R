# Demo code for examining performance of MCF method

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
source("Helper/get_prev.R")
source("Helper/mcf_detect.R")

# Number of tests
m <- nrow(dat)

# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.04)

# Prepare the result storage
mcf_power <- numeric(length(alpha))
mcf_fdr <- numeric(length(alpha))
# mcf_power <- matrix(0, nrow = length(alpha), ncol = rep)
# mcf_fdr <- matrix(0, nrow = length(alpha), ncol = rep)

# Use the GeneralizedFDREstimators function in fdrDiscreteNull package to get result
res_alpha <- lapply(alpha, FUN = function(x) GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = x))

##########################################################
###################        MCF       #####################
##########################################################    

p_org <- res_alpha[[1]]$pvalues
p_prev <- unlist(lapply(1:m, FUN = function(x, pval, pval_supp) get_prev(pval_supp[[x]], pval[x]), pval = p_org, pval_supp = res_alpha[[1]]$pvalSupp))

# Number of replications used to simulate expected proportion of rejections by Habiger's method
B <- 200
#B_ecdf <- B*m
# Simulate the randomized p-values as in Habiger's method
replicated_randp <- replicate(B, runif(m, p_prev, p_org))
randp_ecdf <- as.vector(replicated_randp)
lambda_storey <- 0.5
pi0 <- mean(apply(replicated_randp, 2, FUN = function(x) (1 + sum(x > lambda_storey))/((1 - lambda_storey)*m)))

# Thresholding p-values from Habiger's method
lambda_star <- unlist(lapply(res_alpha, FUN = function(x) x$SARP$Threshold))
# Get detections (rejections)
mcf_detection <- lapply(lambda_star, FUN = function(x, p_prev, p_org, randp_ecdf) mcf_detect(x, p_prev, p_org, randp_ecdf), p_prev = p_prev, p_org = p_org, randp_ecdf = randp_ecdf)

# Prepare results for plotting
mcf_power <- unlist(lapply(mcf_detection, FUN = function(x) get_power(x, non_null_idx)))
mcf_fdr <- unlist(lapply(mcf_detection, FUN = function(x) get_fdr(x, true_null_idx)))

summ <- data.frame(matrix(0, ncol = 4, nrow = length(alpha)))
summ[, 1] <- alpha
summ[, 2] <- mcf_power
summ[, 3] <- mcf_fdr
summ[, 4] <- rep("MCF", length(alpha))

colnames(summ) <- c("alpha", "power", "fdr", "method")

power_plot <- ggplot(data = summ, aes(x = alpha, y = power)) + geom_line(aes(color = method, linetype = method)) + xlab("Nominal FDR level") + ylab("True power")
power_plot
fdr_plot <- ggplot(data = summ, aes(x = alpha, y = fdr)) + geom_line(aes(color = method, linetype = method)) + geom_abline(slope = 1, intercept = 0) + xlab("Nominal FDR level") + ylab("True FDR")
fdr_plot





