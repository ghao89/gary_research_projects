###############################################################
##############            Try dmrseq          #################
###############################################################

library(dmrseq)

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
                     cutoff = 0.05,
                     smooth = FALSE)
proc.time() - t

dmrseq_dmr$index[dmrseq_dmr$qval <= 0.05]

intersect(dmrseq_dmr$index[dmrseq_dmr$qval <= 0.05], non_null_range)
