permutation_null <- function(dat, 
                             B = NULL) {
  if (ncol(dat) %% 4 != 0) {
    stop("The number of columns contained in the data should be a multiple of 4!")
  }
  
  n <- ncol(dat)/2
  if (n <= 10) {
    
    perms <- combn(2*n, n)
    if (is.null(B) || B > ncol(perms)) {
      if (B > ncol(perms)) {
        warning("B is larger than the number of all possible permutations.")
      }
      warning("Use all permutations.")
      z <- apply(perms, 2, FUN = function(x, y) 
                                        calc(y, x),y = dat)
      return(z)
      
    } else {
      B_perms <- perms[, sample(1:ncol(perms), B)]
      z <- unlist(lapply(1:B, FUN = function(x, y, perms) 
                                        calc(y, perms[, x]),
                          y = dat, perms = B_perms
                        )
                  )
      return(z)
    }
  } else {
    if (is.null(B)) {
      warning("Number of samples in each group is too large (> 10). Use default B = 100 to randomly sample from all possible permutations.")
      B <- 100
    }

    z <- unlist(lapply(1:B, FUN = function(x, y) 
                                        calc(y, sample(1:(2*n), n)), y = dat
                      )
                )
    return(z)
  }

}

calc <- function(y, 
                 idx) {
  
  control <- y[, idx]
  treat <- y[, -idx]
  n <- ncol(control)
  t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
    ifelse(var(control[x,]) == 0 && var(treat[x,]) == 0, 0, (mean(control[x, ]) - mean(treat[x, ]))/sqrt(var(control[x, ])/n + var(treat[x, ])/n)),
    control = control, treat = treat))
  z_stat <- qnorm(pt(t_stat, df = 2*n - 2))  
  return(z_stat)
  
}