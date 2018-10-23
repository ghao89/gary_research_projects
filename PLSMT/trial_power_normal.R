### Trials for comparing statistical power between t-test and test through balanced permutation
### This simulation assumes that both groups follow normal distribution with sd = 1

# Number of replications
rep <- 1000
# Number of samples drawn to get empirical distribution
m <- 10000
# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.02)
# Number of samples in each group
n <- 4
# Normal mean for each group
mu1 <- 0
mu2 <- 1
# Standard deviation for the normal distribution for each group
sigma <- 1

# Simulate control group and treatment group and repeat 'rep' times
control <- replicate(rep, rnorm(n, mu1, sigma))
treat <- replicate(rep, rnorm(n, mu2, sigma))
dat <- rbind(control, treat)

# Construct all the permutations (half when each group has equal number of samples based on symmetry)
# Randomly pick one from all the balanced permutations
perms <- combn(2*n, n)
perms <- perms[, seq_len(ncol(perms)/2)]
idx <- which(apply(perms, 2, FUN = function(x) sum(x %in% (1:n)) == floor(n/2)))
perm <- perms[, sample(idx, 1)]

##########################################################
###############     Two-sample t-test     ################    
##########################################################

stat_t_test <- unlist(lapply(1:rep, FUN = function(x, control, treat) 
  t.test(control[, x], treat[, x])$statistic, control = control, treat = treat))
p_val_t_test <- unlist(lapply(1:rep, FUN = function(x, control, treat) 
  t.test(control[, x], treat[, x])$p.value, control = control, treat = treat))
power_t_test <- unlist(lapply(alpha, FUN = function(x, p_val) sum(p_val <= x)/rep, p_val = p_val_t_test))

##########################################################
#######    Heuristic balanced permutation test    ########    
##########################################################
mu2_heuristic <- 0.05
J <- rpois(m, n/4*(mu2_heuristic - mu1)^2)
t_sam <- rt(m, df = 2*n - 2 + 2*J)
# Draw samples from the supposed distribution for the balanced permuted statistic given heuristic alt mean
emp_alt <- sqrt((2*n - 2)/(2*n - 2 + 2*J))*t_sam
alt_ecdf <- ecdf(abs(emp_alt))

# # Control and treatment group after a balanced permutation
# control_perm <- dat[perm, ]
# treat_perm <- dat[-perm, ]
# 
# stat_perm <- unlist(lapply(1:rep, FUN = function(x, control, treat) 
#   t.test(control[, x], treat[, x])$stat, control = control_perm, treat = treat_perm))

p_val_perm <- 1 - alt_ecdf(abs(stat_t_test))
power_perm <-  unlist(lapply(alpha, FUN = function(x, p_val) sum(p_val <= x)/rep, p_val = p_val_perm))

##########################################################
#######      Compare the power of two methods     ########    
##########################################################
library(ggplot2)
plot_data <- data.frame(matrix(0, ncol = 3, nrow = 2*length(alpha)))
plot_data[, 1] <- rep(alpha, 2)
plot_data[, 2] <- c(power_t_test, power_perm)
plot_data[, 3] <- c(rep("t_test", length(alpha)), rep("perm", length(alpha)))

colnames(plot_data) <- c("alpha", "power", "method")

ggplot(data = plot_data, aes(x = alpha, y = power, color = method)) + geom_line()


