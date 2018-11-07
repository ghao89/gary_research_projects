permutation_null <- function(dat, 
                             B = NULL,
                             stat = "z_stat",
                             perfect = FALSE) {
  calc_z <- function(y, 
                   idx,
                   perfect = FALSE) {
    control <- y[, idx]
    treat <- y[, -idx]
    n <- ncol(control)
    
    if (perfect) {
      t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
        ifelse(var(control[x,]) == 0 && var(treat[x,]) == 0, 0, (mean(control[x, ]) - mean(treat[x, ]))/sqrt(var(control[x, ])/n + var(treat[x, ])/n)),
        control = control, treat = treat))      
    } else {
      t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
        t.test(control[x, ], treat[x, ])$statistic,
        control = control, treat = treat))
    }

    z_stat <- qnorm(pt(t_stat, df = 2*n - 2))  
    return(z_stat)
  }
  
  calc_t <- function(y,
                     idx,
                     perfect = FALSE) {
    control <- y[, idx]
    treat <- y[, -idx]
    n <- ncol(control)
    
    if (perfect) {
      t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
        ifelse(var(control[x,]) == 0 && var(treat[x,]) == 0, 0, (mean(control[x, ]) - mean(treat[x, ]))/sqrt(var(control[x, ])/n + var(treat[x, ])/n)),
        control = control, treat = treat))      
    } else {
      t_stat <- unlist(lapply(1:nrow(control), FUN = function(x, control, treat) 
        t.test(control[x, ], treat[x, ])$statistic,
        control = control, treat = treat))      
    }

    return(t_stat)    
  }
  
  calc <- ifelse(stat == "z_stat", calc_z, calc_t)
  
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
                                        calc(y, x, perfect), y = dat)
      return(z)
      
    } else {
      B_perms <- perms[, sample(1:ncol(perms), B)]
      z <- unlist(lapply(1:B, FUN = function(x, y, perms) 
                                        calc(y, perms[, x], perfect),
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
                                        calc(y, sample(1:(2*n), n), perfect), y = dat
                      )
                )
    return(z)
  }

}
