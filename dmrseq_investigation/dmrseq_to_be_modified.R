function (bs, testCovariate, adjustCovariate = NULL, cutoff = 0.1, 
          minNumRegion = 5, smooth = TRUE, bpSpan = 1000, minInSpan = 30, 
          maxGapSmooth = 2500, maxGap = 1000, verbose = TRUE, maxPerms = 10, 
          matchCovariate = NULL, BPPARAM = bpparam(), stat = "stat", 
          block = FALSE, blockSize = 5000, chrsPerChunk = 1) 
{
  stopifnot(is(bs, "BSseq"))
  if (!(is.null(cutoff) || length(cutoff) %in% seq_len(2))) 
    stop("'cutoff' has to be either NULL or a vector of length 1 or 2")
  if (length(cutoff) == 2) 
    cutoff <- sort(cutoff)
  if (is.null(cutoff) | abs(cutoff) > 1 | abs(cutoff) == 0) 
    stop("Must specify a value for cutoff between 0 and 1")
  subverbose <- max(as.integer(verbose) - 1L, 0)
  if (minNumRegion < 3) {
    stop("minNumRegion must be at least 3")
  }
  if (!(stat %in% c("L", "area", "beta", "stat", "avg"))) {
    stop("Specified '", stat, "' as the test statistic which is not ", 
         "in the results. Please specify a valid name from one of ", 
         "L, area, beta, stat, or avg")
  }
  if (block) {
    message("Searching for large scale blocks with at least ", 
            blockSize, " basepairs.")
    if (minInSpan < 100 && bpSpan < 2000 && maxGapSmooth < 
        1e+05) {
      warning("When block=TRUE, it is recommended to increase the values ", 
              "of minInSpan, bpSpan, and maxGapSmooth in order to widen ", 
              "the smoothing window")
    }
  }
  if (is.character(testCovariate)) {
    if (length(testCovariate) > 1) 
      stop("Only one testCovariate can be specified")
    if (is.character(adjustCovariate)) {
      if (sum(testCovariate %in% adjustCovariate) > 0) 
        stop("adjustCovariate can't contain testCovariate")
    }
    if (is.character(matchCovariate)) {
      if (sum(testCovariate %in% matchCovariate)) 
        stop("matchCovariate can't contain testCovariate")
    }
    testCovariate <- which(colnames(pData(bs)) == testCovariate)
    if (length(testCovariate) == 0) {
      stop("testCovariate not found in pData(). ", "Please specify a valid testCovariate")
    }
  }
  if (is.character(adjustCovariate)) {
    if (is.character(matchCovariate)) {
      if (matchCovariate == adjustCovariate) 
        stop("matchCovariate can't be identical to adjustCovariate")
    }
    adjustCovariate <- which(colnames(pData(bs)) %in% adjustCovariate)
    if (length(adjustCovariate) == 0) {
      stop("adjustCovariate not found in pData(). ", "Please specify a valid adjustCovariate")
    }
  }
  if (chrsPerChunk != 1) {
    if (chrsPerChunk%%1 != 0) {
      stop("chrsPerChunk must be an integer")
    }
    else if (chrsPerChunk < 1) {
      stop("chrsPerChunk must be strictly positive")
    }
    else if (chrsPerChunk > length(unique(seqnames(bs)))) {
      stop("chrsPerChunk can't be larger than the total", 
           " number of chromosomes")
    }
    else if (!identical(as.character(seqnames(bs)@values), 
                        seqlevels(bs))) {
      stop("BSseq object must be ordered if breaking computation ", 
           "into multiple chromosomes per chunk (see bsseq::orderBSseq())")
    }
  }
  if (ncol(pData(bs)) < max(testCovariate, adjustCovariate)) {
    stop("Error: pData(bs) has too few columns.  ", "\n              Please specify valid ", 
         "covariates to use in the analysis")
  }
  coeff <- seq(2, (2 + length(testCovariate) - 1))
  testCov <- pData(bs)[, testCovariate]
  fact <- TRUE
  if (length(unique(testCov)) == 1) {
    message("Warning: only one unique value of the specified ", 
            "covariate of interest.  Assuming null comparison and ", 
            "splitting sample group into two equal groups")
    testCov <- rep(1, length(testCov))
    testCov[seq_len(round(length(testCov)/2))] <- 0
  }
  else if (length(unique(testCov)) > 2 && !is.numeric(testCov)) {
    message("Performing a global test of H0: no difference among ", 
            length(unique(testCov)), " groups (assuming the test ", 
            "covariate ", colnames(pData(bs))[testCovariate], 
            " is a factor).")
    coeff <- c(coeff, coeff + length(unique(testCov)) - 2)
  }
  else if (length(unique(testCov)) > 2 && is.numeric(testCov)) {
    message("Assuming the test ", "covariate ", colnames(pData(bs))[testCovariate], 
            " is continuous.")
    fact <- FALSE
  }
  else {
    message("Assuming the test ", "covariate ", colnames(pData(bs))[testCovariate], 
            " is a factor.")
    if (min(table(testCov)) < 2) 
      stop("At least one group has only one sample! ", 
           "Replicates are required to run dmrseq.")
    testCov <- as.factor(testCov)
  }
  sampleSize <- table(testCov)
  if (!is.null(adjustCovariate)) {
    mmdat <- data.frame(testCov = testCov)
    adjustCov <- pData(bs)[, adjustCovariate, drop = FALSE]
    mmdat <- cbind(mmdat, adjustCov)
    frm <- paste0("~", paste0(colnames(mmdat), collapse = " + "))
    design <- model.matrix(as.formula(frm), data = mmdat)
    colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
    colnames(design)[seq((max(coeff) + 1), ncol(design))] <- colnames(pData(bs))[adjustCovariate]
    coeff.adj <- which(colnames(design) == colnames(pData(bs))[adjustCovariate])
  }
  else {
    design <- model.matrix(~testCov)
    colnames(design)[coeff] <- colnames(pData(bs))[testCovariate]
    coeff.adj <- NULL
  }
  if (length(coeff) > 1 && any(rowSums(design[, coeff]) > 1)) 
    stop("Interaction terms in testCovariate are not yet supported.")
  if (length(unique(testCov)) == 2) {
    message("Condition: ", unique(pData(bs)[, testCovariate][which(design[, 
                                                                          coeff] == 1)]), " vs ", unique(pData(bs)[, testCovariate][which(design[, 
                                                                                                                                                 coeff] == 0)]))
  }
  if (!is.null(adjustCovariate)) {
    message("Adjusting for covariate: ", paste(colnames(pData(bs))[adjustCovariate], 
                                               collapse = ", "))
  }
  if (!is.null(matchCovariate)) {
    if (length(matchCovariate) > 1) 
      stop("Covariate matching can only be carried out for one", 
           " covariate")
    if (length(unique(testCov)) > 2) 
      stop("Covariate matching can only be carried out for 2-group", 
           " comparisons")
    if (is.character(matchCovariate)) {
      if (sum(grepl(matchCovariate, colnames(pData(bs)))) == 
          0) {
        stop("Error: no column in pData() found that matches ", 
             "the matchCovariate")
      }
      else if (length(grep(matchCovariate, colnames(pData(bs)))) > 
               1) {
        stop("Error: matchCovariate matches more than one ", 
             "column in pData()")
      }
      mC <- grep(matchCovariate, colnames(pData(bs)))
    }
    else {
      stopifnot(matchCovariate <= ncol(pData(bs)))
    }
    message("Matching permutations on covariate: ", colnames(pData(bs))[mC])
  }
  if (fact) {
    lev <- unique(pData(bs)[[testCovariate]])
    filter <- NULL
    for (l in seq_along(lev)) {
      filter <- rbind(filter, 1 * (DelayedMatrixStats::rowSums2(getCoverage(bs)[, 
                                                                                pData(bs)[[testCovariate]] == lev[l]]) == 0))
    }
    filter <- which(apply(filter, 2, max) > 0)
    if (length(filter) > 0) {
      stop(length(filter), " loci have zero coverage in all samples ", 
           "of at least one condition. Please remove these loci ", 
           "before running dmrseq")
    }
  }
  else {
    filter <- DelayedMatrixStats::rowSums2(getCoverage(bs) == 
                                             0) >= ncol(bs) - 1
    if (sum(filter) > 0) 
      stop(sum(filter), " loci have zero coverage in at least ", 
           ncol(bs) - 1, " samples. Please remove these loci ", 
           "before running dmrseq")
  }
  BiocParallel::register(BPPARAM)
  backend <- paste0("BiocParallel:", class(bpparam())[1])
  if (bpparam()$workers == 1) {
    if (verbose) {
      mes <- "Using a single core (backend: %s)."
      message(sprintf(mes, backend))
    }
    parallel <- FALSE
  }
  else {
    if (verbose) {
      mes <- paste0("Parallelizing using %s workers/cores ", 
                    "(backend: %s).")
      message(sprintf(mes, bpparam()$workers, backend))
    }
    parallel <- TRUE
  }
  message("Computing on ", chrsPerChunk, " chromosome(s) at a time.\n")
  message("Detecting candidate regions with coefficient larger than ", 
          unique(abs(cutoff)), " in magnitude.")
  OBS <- bumphunt(bs = bs, design = design, coeff = coeff, 
                  coeff.adj = coeff.adj, minInSpan = minInSpan, minNumRegion = minNumRegion, 
                  cutoff = cutoff, maxGap = maxGap, maxGapSmooth = maxGapSmooth, 
                  smooth = smooth, bpSpan = bpSpan, verbose = verbose, 
                  parallel = parallel, block = block, blockSize = blockSize, 
                  chrsPerChunk = chrsPerChunk, fact = fact)
  if (length(OBS) > 0) {
    message("* ", nrow(OBS), " candidates detected")
    FLIP <- NULL
    if (length(unique(design[, coeff[1]])) == 2 && length(coeff) == 
        1 && choose(nrow(design), min(sampleSize)) < 5e+05) {
      if (verbose) {
        message("Performing balanced permutations of ", 
                "condition across samples ", "to generate a null distribution of region test statistics")
      }
      perms <- combn(seq(1, nrow(design)), min(sampleSize))
      if (length(unique(table(design[, coeff]))) == 1) {
        perms <- perms[, seq_len(ncol(perms)/2)]
      }
      rmv <- NULL
      for (p in seq_len(ncol(perms))) {
        if (length(unique(design[perms[, p], coeff])) == 
            1) {
          rmv <- c(rmv, p)
        }
      }
      if (length(rmv) > 0) 
        perms <- perms[, -rmv]
      if (maxPerms < ncol(perms)) {
        similarity <- apply(perms, 2, function(x) {
          max(table(design[x, coeff]))
        })
        perms.all <- perms
        perms <- NULL
        levs <- sort(unique(similarity))
        l <- 1
        num <- 0
        while (!(num == maxPerms) && l <= length(levs)) {
          keep <- sample(which(similarity == levs[l]), 
                         min(maxPerms - num, sum(similarity == levs[l])))
          perms <- cbind(perms, perms.all[, keep])
          l <- l + 1
          num <- ncol(perms)
        }
      }
    }
    else {
      if (verbose) {
        message("Performing unrestricted permutation of", 
                " covariate of interest across samples ", "to generate a null distribution of region test statistics")
      }
      perms <- as.matrix(seq_len(nrow(design)))
      for (p in seq_len(maxPerms)) {
        tries <- 0
        candidate <- sample(seq_len(nrow(design)), nrow(design))
        while ((sum(apply(perms, 2, function(x) all.equal(x, 
                                                          candidate)) == TRUE) > 0 || sum(apply(perms, 
                                                                                                2, function(x) all.equal(x, rev(candidate))) == 
                                                                                          TRUE) > 0) && tries <= 20) {
          candidate <- sample(seq(seq_len(nrow(design))), 
                              nrow(design))
          tries <- tries + 1
        }
        if (tries <= 20) {
          perms <- cbind(perms, candidate)
        }
      }
      perms <- perms[, -1]
    }
    pData.orig <- pData(bs)
    levs <- unique(pData.orig[[testCovariate]])
    for (j in seq_len(ncol(perms))) {
      if (verbose) {
        message("\nBeginning permutation ", j)
      }
      reorder <- perms[, j]
      designr <- design
      if (length(unique(design[, coeff[1]])) == 2 && length(coeff) == 
          1 && !nrow(perms) == nrow(designr)) {
        designr[, coeff] <- 0
        designr[reorder, coeff] <- 1
        pData(bs)[[testCovariate]] <- levs[1]
        pData(bs)[[testCovariate]][reorder] <- levs[2]
        if (!all(sort(pData.orig[[testCovariate]]) == 
                 sort(pData(bs)[[testCovariate]]))) {
          designr[, coeff] <- 1
          designr[reorder, coeff] <- 0
          pData(bs)[[testCovariate]] <- levs[2]
          pData(bs)[[testCovariate]][reorder] <- levs[1]
        }
        xr <- NULL
        for (rd in seq_len(nrow(pData.orig))) {
          match <- which(pData.orig[[testCovariate]] %in% 
                           pData(bs)[rd, ][[testCovariate]])
          taken <- which(match %in% xr)
          if (length(taken) > 0) 
            match <- match[-taken]
          if (length(match) > 0) 
            xr <- c(xr, match[1])
        }
        if (length(coeff.adj) > 0) {
          pData(bs)[, adjustCovariate] <- pData.orig[xr, 
                                                     adjustCovariate]
        }
      }
      else {
        designr[, coeff] <- designr[reorder, coeff]
        pData(bs) <- pData.orig[reorder, , drop = FALSE]
      }
      if (!is.null(matchCovariate)) {
        permLabel <- paste0(paste0(pData(bs)[designr[, 
                                                     coeff[1]] == 1, mC], collapse = "_"), "vs", 
                            paste0(pData(bs)[(1 - designr[, coeff[1]]) == 
                                               1, mC], collapse = "_"))
        c1 <- unlist(strsplit(permLabel, "vs"))[1]
        c2 <- unlist(strsplit(permLabel, "vs"))[2]
        c1 <- unlist(strsplit(c1, "_"))
        c2 <- unlist(strsplit(c2, "_"))
        keepPerm <- 1 * (sum(c1 %in% c2) > 0 && sum(c2 %in% 
                                                      c1) > 0)
        if (keepPerm == 0) {
          if (verbose) {
            message(paste0("Skipping permutation ", gsub("vs", 
                                                         " vs ", permLabel)))
          }
          next
        }
      }
      else {
        permLabel <- j
      }
      res.flip.p <- bumphunt(bs = bs, design = designr, 
                             coeff = coeff, coeff.adj = coeff.adj, minInSpan = minInSpan, 
                             minNumRegion = minNumRegion, cutoff = cutoff, 
                             maxGap = maxGap, maxGapSmooth = maxGapSmooth, 
                             smooth = smooth, bpSpan = bpSpan, verbose = verbose, 
                             parallel = parallel, block = block, blockSize = blockSize, 
                             chrsPerChunk = chrsPerChunk, fact = fact)
      if (verbose) {
        message("* ", j, " out of ", ncol(perms), " permutations completed (", 
                nrow(res.flip.p), " null candidates)")
      }
      if (!is.null(res.flip.p)) {
        res.flip.p$permNum <- permLabel
        FLIP <- rbind(FLIP, res.flip.p)
      }
    }
    pData(bs) <- pData.orig
    if (is.null(FLIP)) {
      warning("No candidate regions found in permutation, so inference ", 
              "can't be carried out. ", "Try decreasing the cutoff, or running on a larger ", 
              "dataset if you are currently using a subset.")
      OBS$pval <- NA
      OBS$qval <- NA
    }
    else if (nrow(FLIP) < 0.05 * nrow(OBS)) {
      message("Note: Very few null candidate regions were found.", 
              "For more accurate and sensitive inference, ", 
              "try decreasing the cutoff, or running on a larger ", 
              "dataset if you are currently using a subset.")
    }
    if (!is.null(FLIP)) {
      if (nrow(FLIP) > 1e+06) {
        rs <- sample(seq_len(nrow(FLIP)), 1e+06, replace = FALSE)
        FLIP <- FLIP[rs, ]
      }
      if (!(stat %in% c(colnames(OBS), "avg"))) {
        stop("Specified '", stat, "' as the test statistic which is not ", 
             "in the results. Please specify a valid name from one of ", 
             "L, area, beta, or stat")
      }
      else if (stat == "avg") {
        OBS$avg <- OBS$area/OBS$L
        FLIP$avg <- FLIP$area/FLIP$L
      }
      whichStatO <- which(colnames(OBS) == stat)
      whichStatF <- which(colnames(FLIP) == stat)
      perm.ordered <- c(sort(abs(FLIP[, whichStatF]), method = "quick"), 
                        Inf)
      pval <- rep(NA, nrow(OBS))
      pval[!is.na(OBS[, whichStatO])] <- (1 + vapply(abs(OBS[!is.na(OBS[, 
                                                                        whichStatO]), whichStatO]), function(x) length(perm.ordered) - 
                                                       min(which(x <= perm.ordered)), numeric(1)))/(1 + 
                                                                                                      sum(!is.na(FLIP[, whichStatF])))
      pval[abs(pval) == Inf] <- NA
      pval <- data.frame(x = pval, y = p.adjust(pval, method = "BH"))
      OBS$pval <- pval$x
      OBS$qval <- pval$y
    }
    indexIR <- IRanges(OBS$indexStart, OBS$indexEnd)
    OBS.gr <- makeGRangesFromDataFrame(OBS[, -c(4:5)], keep.extra.columns = TRUE)
    OBS.gr$index <- indexIR
    names(OBS.gr) <- NULL
    OBS.gr <- OBS.gr[order(OBS.gr$pval, -abs(OBS.gr$stat)), 
                     ]
    return(OBS.gr)
  }
  else {
    message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    return(NULL)
  }
}