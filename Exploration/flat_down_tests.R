# Flat the m tests into 1 test

dat <- binom/n_trial

n <- ncol(dat)/2
left_mean <- apply(dat[, 1:n], 1, mean)
dat <- dat - left_mean

control <- as.vector(dat[sample(1:nrow(control), 50), 1:n])
treat <- as.vector(dat[sample(1:nrow(treat), 50), (n + 1):(2*n)])

B <- 10000
y <- c(control, treat)

calc <- function(y, idx) {
  control <- y[idx]
  treat <- y[-idx]
  n <- length(control)
  t_stat <- (mean(control) - mean(treat))/sqrt(var(control)/n + var(treat)/n)
  z_stat <- qnorm(pt(t_stat, df = 2*n - 2))  
  return(z_stat)
}

z <- unlist(lapply(1:B, FUN = function(x, y) 
  calc(y, sample(1:(2*length(control)), length(control))), y = y
))

plot(density(z))
lines(seq(-4, 4, by = 0.1), dnorm(seq(-4, 4, by = 0.1)), col = "red")
