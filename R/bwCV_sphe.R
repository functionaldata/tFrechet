# using CV to choose bw for local Fréchet regression on a unit hypersphere.
bwCV_sphe <- function(xin, yin, xout, optns) {
  yin <- yin[order(xin),]
  xin <- sort(xin)
  compareRange <- (xin > min(xin) + diff(range(xin))/5) & (xin < max(xin) - diff(range(xin))/5)
  
  # k-fold
  objFctn <- function(bw) {
    optns1 <- optns
    optns1$bw <- bw
    folds <- numeric(length(xin))
    n <- sum(compareRange)
    numFolds <- ifelse(n > 30, 10, sum(compareRange))
    
    tmp <- c(sapply(1:ceiling(n/numFolds), function(i)
      sample(x = seq_len(numFolds), size = numFolds, replace = FALSE)))
    tmp <- tmp[1:n]
    repIdx <- which(diff(tmp) == 0)
    for (i in which(diff(tmp) == 0)) {
      s <- tmp[i]
      tmp[i] <- tmp[i-1]
      tmp[i-1] <- s
    }
    #tmp <- cut(1:n,breaks = seq(0,n,length.out = numFolds+1), labels=FALSE)
    #tmp <- tmp[sample(seq_len(n), n)]
    
    folds[compareRange] <- tmp
    
    yout <- lapply(seq_len(numFolds), function(foldidx) {
      testidx <- which(folds == foldidx)
      res <- LocSpheGeoReg(xin = xin[-testidx], yin = yin[-testidx,], xout = xin[testidx], optns = optns1)
      res # each row is a spherical vector
    })
    yout <- do.call(rbind, yout)
    yinMatch <- yin[which(compareRange)[order(tmp)],]
    mean(sapply(1:nrow(yout), function(i) SpheGeoDist(yout[i,], yinMatch[i,])^2))
  }
  bwRange <- SetBwRange(xin = xin, xout = xout, kernel_type = optns$ker)
  #if (!is.null(optns$bwRange)) {
  #  if (min(optns$bwRange) < bwRange$min) {
  #    message("Minimum bandwidth is too small and has been reset.")
  #  } else bwRange$min <- min(optns$bwRange)
  #  if (max(optns$bwRange) >  bwRange$min) {
  #    bwRange$max <- max(optns$bwRange)
  #  } else {
  #    message("Maximum bandwidth is too small and has been reset.")
  #  }
  #}
  res <- optimize(f = objFctn, interval = c(bwRange$min, bwRange$max))
  res$minimum
}
