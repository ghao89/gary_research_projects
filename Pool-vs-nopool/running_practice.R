library(ggplot2)
library(gridExtra)

rep <- 2
alpha <- seq(0, 0.5, by = 0.01)

power_t <- matrix(0, nrow = rep, ncol = length(alpha))
fdr_t <- matrix(0, nrow = rep, ncol = length(alpha))
power_pool <- matrix(0, nrow = rep, ncol = length(alpha))
fdr_pool <- matrix(0, nrow = rep, ncol = length(alpha))
power_pool_reinf <- matrix(0, nrow = rep, ncol = length(alpha))
fdr_pool_reinf <- matrix(0, nrow = rep, ncol = length(alpha))

for (j in 1:rep) {
  
    # Total number of tests
    m <- 5000
    # number of samples in each group
    n <- 5
    # Proportion of null
    pi_0 <- 0.75
    
    #----------Simulation setup (Normal)----------#
    true_null_idx <- 1:(m*pi_0)
    non_null_idx <- (1:m)[-true_null_idx]

    # Mean of each test in control group and treatment group
    mu_control <- numeric(m)
    mu_treat <- numeric(m)

    mu_treat[true_null_idx] <- mu_control[true_null_idx]
    mu_treat[non_null_idx] <- runif(length(non_null_idx), 1, 4)

    # Standard deviation of each test in both groups
    sigma_control <- rep(1, m)
    sigma_treat <- sigma_control
    sigma_treat[non_null_idx] <- rexp(length(non_null_idx), rate = 1)

    # Generate data
    control <- replicate(n, rnorm(m, mu_control, sigma_control))
    treat <- replicate(n, rnorm(m, mu_treat, sigma_treat))

    dat <- cbind(control, treat)
    
    # #----------Simulation setup (Cauchy)----------#
    # # Index of true null
    # true_null_idx <- 1:(m*pi_0)
    # non_null_idx <- (1:m)[-true_null_idx]
    # 
    # # Mean of each test in control group and treatment group
    # loc_control <- numeric(m)
    # loc_treat <- numeric(m)
    # 
    # loc_treat[true_null_idx] <- loc_control[true_null_idx]
    # loc_treat[non_null_idx] <- runif(length(non_null_idx), 5, 10)
    # 
    # # Standard deviation of each test in both groups
    # #scale_control <- rep(1, m)
    # scale_control <- rep(1, m)
    # scale_treat <- scale_control
    # scale_treat[non_null_idx] <- rexp(length(non_null_idx), rate = 1)
    # 
    # # Generate data
    # control <- replicate(n, rcauchy(m, loc_control, scale_control))
    # treat <- replicate(n, rnorm(m, loc_treat, scale_treat))
    # 
    # dat <- cbind(control, treat)
    
    #----------Conduct multiple testing----------#
    perms <- combn(2*n, n)
    perms <- perms[, 1:(ncol(perms)/2)]
    stat <- matrix(0, nrow = m, ncol = ncol(perms))
    
    for (i in 1:ncol(perms)) {
      perm <- perms[, i]
      control <- dat[, perm]
      treat <- dat[, -perm]
      
      stat[, i] <- unlist(lapply(1:nrow(treat),
                                 FUN = function(x, treat, control)
                                   t.test(control[x, ], treat[x, ])$statistic,
                                 treat = treat,
                                 control = control))
    }
    
    #------------------Individual t-test------------------#
    p_val_t <- apply(dat, 
                     1, 
                     FUN = function(x)
                       t.test(x[1:n], x[n + (1:n)])$p.value)
    q_val_t <- p.adjust(p_val_t, method = "BH")
    power_t[j, ] <- unlist(lapply(alpha, 
                             FUN = function(x, q_val, non_null_idx)
                               sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
                             q_val = q_val_t,
                             non_null_idx = non_null_idx))
    fdr_t[j, ] <- unlist(lapply(alpha, 
                           FUN = function(x, q_val, non_null_idx)
                             sum(!(which(q_val <= x) %in% non_null_idx))/max(sum(q_val <= x), 1),
                           q_val = q_val_t,
                           non_null_idx = non_null_idx))
    
    # #------------------Pooled permutation test------------------#
    # ref_cdf <- ecdf(abs(as.vector(stat[, 2:ncol(stat)])))
    # p_val_pool <- 1 - ref_cdf(abs(stat[, 1]))
    # q_val_pool <- p.adjust(p_val_pool, method = "BH")
    # 
    # power_pool[j, ] <- unlist(lapply(alpha, 
    #                             FUN = function(x, q_val, non_null_idx)
    #                               sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
    #                             q_val = q_val_pool,
    #                             non_null_idx = non_null_idx))
    # fdr_pool[j, ] <- unlist(lapply(alpha, 
    #                           FUN = function(x, q_val, non_null_idx)
    #                             sum(!(which(q_val <= x) %in% non_null_idx))/max(sum(q_val <= x), 1),
    #                           q_val = q_val_pool,
    #                           non_null_idx = non_null_idx))
    
    #------------------Pooled permutation test with transformation------------------#
    ref_cdf <- ecdf(qnorm(pt(as.vector(stat[true_null_idx, 2:ncol(stat)]), df = 2*n - 2)))
    p_val_pool <- 2*(0.5 - abs(0.5 - ref_cdf(qnorm(pt(stat[, 1], df = 2*n - 2)))))
    q_val_pool <- p.adjust(p_val_pool, method = "BH")

    power_pool[j, ] <- unlist(lapply(alpha,
                                FUN = function(x, q_val, non_null_idx)
                                  sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
                                q_val = q_val_pool,
                                non_null_idx = non_null_idx))
    fdr_pool[j, ] <- unlist(lapply(alpha,
                              FUN = function(x, q_val, non_null_idx)
                                sum(!(which(q_val <= x) %in% non_null_idx))/max(sum(q_val <= x),1),
                              q_val = q_val_pool,
                              non_null_idx = non_null_idx))
    
    detected_idx <- non_null_idx
    reinf_ref_cdf <- ecdf(abs(as.vector(stat[-detected_idx, 2:ncol(stat)])))
    p_val_pool_reinf <- 1 - reinf_ref_cdf(abs(stat[, 1]))
    q_val_pool_reinf <- p.adjust(p_val_pool_reinf, method = "BH")
    
    power_pool_reinf[j, ] <- unlist(lapply(alpha,
                                      FUN = function(x, q_val, non_null_idx)
                                        sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
                                      q_val = q_val_pool_reinf,
                                      non_null_idx = non_null_idx))
    
    fdr_pool_reinf[j, ] <- unlist(lapply(alpha,
                                    FUN = function(x, q_val, non_null_idx)
                                      sum(!(which(q_val <= x) %in% non_null_idx))/max(sum(q_val <= x), 1),
                                    q_val = q_val_pool_reinf,
                                    non_null_idx = non_null_idx))
    print(paste0("The ", j, "th iteration is done!"))
    
}

# power_pool_storey <- unlist(lapply(alpha,
#                               FUN = function(x, q_val, non_null_idx)
#                                 sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
#                               q_val = q_val_pool_storey,
#                               non_null_idx = non_null_idx))
# fdr_pool_storey <- unlist(lapply(alpha,
#                             FUN = function(x, q_val, non_null_idx)
#                               sum(!(which(q_val <= x) %in% non_null_idx))/sum(q_val <= x),
#                             q_val = q_val_pool_storey,
#                             non_null_idx = non_null_idx))

power_pool_mean <- apply(power_pool, 2, mean)
fdr_pool_mean <- apply(fdr_pool, 2, mean)
power_t_mean <- apply(power_t, 2, mean)
fdr_t_mean <- apply(fdr_t, 2, mean)
power_pool_reinf_mean <- apply(power_pool_reinf, 2, mean)
fdr_pool_reinf_mean <- apply(fdr_pool_reinf, 2, mean)

#----------Plot----------#
dat_plot <- cbind.data.frame(c(power_pool_mean, power_pool_reinf_mean, power_t_mean),
                             c(fdr_pool_mean, fdr_pool_reinf_mean, fdr_t_mean),
                             c(alpha, alpha, alpha),
                             c(rep("Pooled", length(alpha)), 
                               rep("Pooled perfect", length(alpha)), 
                               rep("Individual t", length(alpha))))
colnames(dat_plot) <- c("power", "fdr", "alpha", "method")
power_plot <- ggplot(data = dat_plot) + geom_line(aes(x = alpha, y = power, color = method)) + ggtitle("Power")
fdr_plot <- ggplot(data = dat_plot) + geom_line(aes(x = alpha, y = fdr, color = method)) + geom_abline(slope = 1, intercept = 0) + ggtitle("FDR")
grid.arrange(power_plot, fdr_plot, nrow = 1, top = paste0("Comparison when n = ", n, ", pi_0 = ", pi_0))
print(paste0("The step for n = ", n, ", pi_0 = ", pi_0, " has completed!"))

