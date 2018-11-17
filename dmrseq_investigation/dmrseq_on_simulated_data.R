###############################################################
##############            Try dmrseq          #################
###############################################################

library(dmrseq)

alpha <- 0.05
cutoff <- seq(0.01, 0.1, by = 0.01)
fp_bef_merge <- numeric(length(cutoff))
tp_merge <- numeric(length(cutoff))
num_detect <- numeric(length(cutoff))


for (c in cutoff) {
  idx <- which(cutoff == c)
  load(file = "~/Dropbox/Private/Git_Projects/DLSMT/dmrseq_investigation/bs_dat.Rdata")
  
  n <- length(sampleNames(bs_dat))/2
  
  perms <- combn(2*n, n)
  i <- 1
  perm <- c(perms[, i], perms[, ncol(perms) - i + 1])
  p_data <- data.frame(cbind(c(rep("Control",n), rep("Treat", n))[perm],
                             sampleNames(bs_dat)[perm]))
  colnames(p_data) <- c("Group", "Rep")
  pData(bs_dat) <- p_data
  
  t <- proc.time()
  dmrseq_dmr <- dmrseq(bs = bs_dat, 
                       testCovariate = "Group",
                       cutoff = c,
                       smooth = FALSE)
  proc.time() - t
  
  detected <- dmrseq_dmr$index[dmrseq_dmr$qval <= alpha]
  num_detect[idx] <- length(detected)
  tp_bef_merge <- intersect(dmrseq_dmr$index[dmrseq_dmr$qval <= alpha], non_null_range)
  tp_merge[idx] <- sum(unique(floor(start(tp_bef_merge)/100)) %in% (non_null_idx - 1))
  fp_bef_merge[idx] <- length(detected) - length(tp_bef_merge)
  
}
