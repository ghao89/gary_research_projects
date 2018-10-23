### Trials for comparing statistical power between t-test and test through balanced permutation
### This simulation assumes that both groups follow beta distribution

# Number of replications
rep <- 1000
# Number of samples drawn to get empirical distribution
m <- 10000
# Nominal significance level
alpha <- seq(0.01, 0.25, by = 0.02)
# Number of samples in each group
n <- 4
# Shape 1 for both groups
a <- 10
# Shape 2 for each group
b_1 <- 10
b_2 <- 10

# Simulate control group and treatment group and repeat 'rep' times
control <- replicate(rep, rbeta(n, a, b_1))
treat <- replicate(rep, rbeta(n, a, b_2))
control <- apply(control, 2, FUN = function(x) x/sd(x))
treat <- apply(treat, 2, FUN = function(x) x/sd(x))
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
# False rejections
fr_t_test <- unlist(lapply(alpha, FUN = function(x, p_val) sum(p_val <= x)/rep, p_val = p_val_t_test))

##########################################################
#######    Heuristic balanced permutation test    ########    
##########################################################
mu2_heuristic <- 0.5
J <- rpois(m, n/4*(mu2_heuristic)^2)
t_sam <- rt(m, df = 2*n - 2 + 2*J)
# Draw samples from the supposed distribution for the balanced permuted statistic given heuristic alt mean
emp_alt <- sqrt((2*n - 2)/(2*n - 2 + 2*J))*t_sam
alt_ecdf <- ecdf(abs(emp_alt))

p_val_perm <- 1 - alt_ecdf(abs(stat_t_test))
# False rejections
fr_perm <-  unlist(lapply(alpha, FUN = function(x, p_val) sum(p_val <= x)/rep, p_val = p_val_perm))

##########################################################
#######      Compare the power of two methods     ########    
##########################################################
library(ggplot2)
plot_data <- data.frame(matrix(0, ncol = 3, nrow = 2*length(alpha)))
plot_data[, 1] <- rep(alpha, 2)
plot_data[, 2] <- c(fr_t_test, fr_perm)
plot_data[, 3] <- c(rep("t_test", length(alpha)), rep("perm", length(alpha)))

colnames(plot_data) <- c("alpha", "false rejections", "method")

ggplot(data = plot_data, aes(x = alpha, y = `false rejections`, color = method)) + geom_line()
