# L2 norm
l2norm <- function(x){
  #sqrt(sum(x^2))
  as.numeric(sqrt(crossprod(x)))
}

# gradient w.r.t. y
SpheGeoGrad <- function(x,y) { #, tol = 1e-10){
  tmp <- 1 - sum(x * y) ^ 2
  return(- (tmp) ^ (-0.5) * x)
  # if (tmp < tol) {
  #   return(- Inf * x)
  # } else {
  #   return(- (tmp) ^ (-0.5) * x)
  # }
}
# Hessian dy dy'
SpheGeoHess <- function(x,y) { #,tol = 1e-10){
  return(- sum(x * y) * (1 - sum(x * y) ^ 2) ^ (-1.5) * x %*% t(x))
}

# using CV to choose bw
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

# exp map
expSphere <- function(base,tg, tol = 1e-10) {
  tgNorm <- l2norm(tg)
  if (!is.na(tgNorm) & tgNorm <= tol) {
    base
  } else {
    sin(tgNorm) * tg / tgNorm + cos(tgNorm) * base
  }
}

# log map
logSphere <- function(base, x, tol = 1e-10) {
  tg <- (x - sum(x * base) * base)
  tgNorm <- l2norm(tg)
  if (!is.na(tgNorm) & tgNorm <= tol) {
    rep(0,length(base))
  } else {
    tg / l2norm(tg) * SpheGeoDist(base,x)
  }
}

# polar to Cartesian
pol2car <- function(phi, theta) {
  data <- data.frame(phi,theta)
  res <- with(data, cbind(
    sin(phi) * cos(theta),
    sin(phi) * sin(theta),
    cos(phi)
  ))
  if (nrow(res) == 1) {
    res <- as.vector(res)
  }
  res
}

# orthonormal bases of tangent plane at (phi, theta)
bss_tgsp <- function(phi, theta) {
  data <- data.frame(phi,theta)
  b1 <- with(data, cbind(
    cos(phi) * cos(theta),
    cos(phi) * sin(theta),
    - sin(phi)
  ))
  b2 <- with(data, cbind(
    sin(theta),
    -cos(theta),
    0
  ))
  if (nrow(b1) == 1) {
    b1 <- as.vector(b1)
    b2 <- as.vector(b2)
  }
  list(b1 = b1, b2 = b2)
}
