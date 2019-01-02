library(fdrDiscreteNull)
library(ggplot2)

m <- 5000
alpha <- seq(0.01, 0.25, by = 0.02)
  
all_mm <- read.table("Randomization-test-and-Empirical-null/all_mm", header = TRUE)
all_mm <- all_mm[, 2:7]

dat <- all_mm[sample(nrow(all_mm), m, replace = TRUE), c(3,5,4,6)]

summ <-  GeneralizedFDREstimators(data = dat, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.01)

p_val <- summ$pvalues
p_val_supp <- summ$pvalSupp
  
  
  
trans_then_randomize <- function(p_val, p_val_supp) {
  
  p_minus <- ifelse(min(p_val_supp) == p_val, 0, p_val_supp[which(p_val_supp == p_val) - 1])
  
  if (p_val == 1) {
    
    if (p_minus != 0) {
      z_left <- qnorm(p_minus/2)
      z_right <- -qnorm(p_minus/2)
      return(runif(1, z_left, z_right))   
    } else {
      return(rnorm(1, 0, 100))
    }

  } else {
    
    if (p_minus != 0) {
      if (runif(1) < 0.5) {
        z_left <- qnorm(p_minus/2)
        z_right <- qnorm(p_val/2)
        return(runif(1, z_left, z_right))
      } else {
        z_left <- -qnorm(p_val/2)
        z_right <- -qnorm(p_minus/2)
        return(runif(1, z_left, z_right))
      }
    } else {
      draw <- rnorm(1, 0, 100)
      return(draw - sign(draw)*qnorm(p_val/2))
    }
  }
}

zz <- unlist(lapply(1:length(p_val), FUN = function(x, p_val, p_val_supp)
              trans_then_randomize(p_val[x], p_val_supp[[x]]),
              p_val = p_val,
              p_val_supp = p_val_supp))