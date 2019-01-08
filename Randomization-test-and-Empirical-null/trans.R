trans_to_zval <- function(rnd_p_val) {
  return(((runif(length(rnd_p_val)) < 0.5)*2 - 1)*qnorm(rnd_p_val/2))
}