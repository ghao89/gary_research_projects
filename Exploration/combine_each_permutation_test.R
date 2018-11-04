dat <- binom/n_trial

B <- 10

calc <- function(y, idx) {
  control <- y[, idx]
  treat <- y[, -idx]
  n <- ncol(control)
  t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
    ifelse(var(control[x,]) == 0 && var(treat[x,]) == 0, 0, (mean(control[x, ]) - mean(treat[x, ]))/sqrt(var(control[x, ])/n + var(treat[x, ])/n)),
    control = control, treat = treat))
  z_stat <- qnorm(pt(t_stat, df = 2*n - 2))  
  return(z_stat)
}

perms <- combn(ncol(dat), ncol(dat)/2)
B_perms <- perms[, sample(1:ncol(perms), B)]

z <- unlist(lapply(1:B, FUN = function(x, y, perms) 
  #calc(y, sample(1:ncol(y), ncol(y)/2)), 
  calc(y, perms[, x]),
  y = dat, perms = B_perms
))

plot(density(z))
lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = "red")
