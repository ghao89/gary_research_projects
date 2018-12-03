# Comparison under normal setup

n <- ncol(dat)/2
m <- nrow(dat)

perms <- combn(2*n, n)
perms <- perms[, 1:(ncol(perms)/2)]
stat <- matrix(0, nrow = nrow(dat), ncol = ncol(perms))

for (i in 1:ncol(perms)) {
  perm <- perms[, i]
  control <- dat[, perm]
  treat <- dat[, -perm]
  
  stat[, i] <- abs(apply(control, 1, mean) - apply(treat, 1, mean))
}

alpha <- seq(0.01, 0.5, by = 0.01)

#------------------Individual permutation test------------------#
p_val <- (ncol(perms) + 1 - rowRanks(stat, ties.method = "min")[,1])/ncol(perms)

p_cand <- (1:ncol(perms))/ncol(perms)

p_count <- unlist(lapply(p_cand, FUN = function(x, p_val) sum(p_val == x), p_val = p_val))
avg_right <- unlist(lapply(p_cand, 
                           FUN = function(x, p_val, p_cand) 
                             sum(p_val > x)/sum(p_cand > x),
                           p_val = p_val,
                           p_cand = p_cand))
lambda <- p_cand[which(avg_right > p_count)[1]]

pi_0_hat <- sum(p_val > lambda)/(1 - lambda)/m

q_val <- unlist(lapply(p_val,
                       FUN = function(x, pi_0, p_val)
                         m*pi_0*x/max(sum(p_val <= x), 1),
                       pi_0 = pi_0_hat,
                       p_val = p_val))

#------------------Individual t-test------------------#
p_val_t <- apply(dat, 
                 1, 
                 FUN = function(x)
                   t.test(x[1:n], x[n + (1:n)])$p.value)
q_val_t <- p.adjust(p_val_t, method = "BH")

power_t <- unlist(lapply(alpha, 
                        FUN = function(x, q_val, non_null_idx)
                          sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
                        q_val = q_val_t,
                        non_null_idx = non_null_idx))
fdr_t <- unlist(lapply(alpha, 
                        FUN = function(x, q_val, non_null_idx)
                          sum(!(which(q_val <= x) %in% non_null_idx))/sum(q_val <= x),
                        q_val = q_val_t,
                        non_null_idx = non_null_idx))

#------------------Pooled permutation test------------------#
ref_cdf <- ecdf(as.vector(stat[, 2:ncol(stat)]))
p_val_pool <- 1 - ref_cdf(stat[, 1])
q_val_pool <- p.adjust(p_val_pool, method = "BH")

power_pool <- unlist(lapply(alpha, 
                    FUN = function(x, q_val, non_null_idx)
                      sum(which(q_val <= x) %in% non_null_idx)/length(non_null_idx),
                    q_val = q_val_pool,
                    non_null_idx = non_null_idx))
fdr_pool <- unlist(lapply(alpha, 
                      FUN = function(x, q_val, non_null_idx)
                        sum(!(which(q_val <= x) %in% non_null_idx))/sum(q_val <= x),
                      q_val = q_val_pool,
                      non_null_idx = non_null_idx))

plot(alpha, power_pool, col = "red", type = "l")
lines(alpha, power_t, col = "blue")
plot(alpha, fdr_pool, col = "red", type = "l")
lines(alpha, fdr_t, col = "blue")
abline(a = 0, b = 0.9)
