transf_thresholds_flexible <- function(gamma){
  if (length(gamma) > 1) {
    cumsum(c(gamma[1], exp(gamma[-1])))
  } else {
    gamma[1]
  }
}


transf_sigma <- function(tpar, ndim) {
  nu <- tpar
  angles <- pi * exp(nu)/(1 + exp(nu))
  cosmat <- diag(ndim)
  cosmat[lower.tri(cosmat)] <- cos(angles)
  S1 <- matrix(0, nrow = ndim, ncol = ndim)
  S1[lower.tri(S1, diag = TRUE)] <- c(rep(1, ndim), sin(angles))
  tLmat <- sapply(1:ndim,
                  function(j) cosmat[j, ] * cumprod(S1[j, ]))
  sigma <- crossprod(tLmat)
  sigma[lower.tri(sigma)]
}


rectbiv_norm_prob <- function(U, L, r) {
  # computes the rectangle probabilities for biv.normal-distribution
  p1 <- pbivnorm(U[,1], U[,2], r)
  p2 <- pbivnorm(L[,1], U[,2], r)
  p3 <- pbivnorm(U[,1], L[,2], r)
  p4 <- pbivnorm(L[,1], L[,2], r)
  ## replace NaN
  p1[is.nan(p1)] <- 0
  p2[is.nan(p2)] <- 0
  p3[is.nan(p3)] <- 0
  p4[is.nan(p4)] <- 0
  pr <- p1 - p2 - p3 + p4
  return(pr)
}







make_start_values <- function(y, X, family_type) {
  start_theta <- lapply(which(family_type == "ordinal"),
                        function(j) coef(ordinal::clm(factor(y[, j]) ~ 1)))
  ndim <- ncol(y)
  p <- ncol(X)
  pars <- c(unlist(lapply(start_theta, function(x) c(x[1], log(diff(x))))),
            # beta0 for normals
            colMeans(y[, family_type != "ordinal", drop = FALSE], na.rm = TRUE),
            # betas
            rep(0, ndim * p),
            # sigmas for normals
            rep(0, sum(family_type != "ordinal")),
            # correlation params
            rep(0, ndim * (ndim - 1)/2))
  return(pars)
}

jac_dttheta_dtheta_flexible <- function(theta, ndimo, ntheta) {
  lapply(seq_len(ndimo), function(j){
    emat <- diag(ntheta[j])
    theta_j <- theta[cumsum(c(0, ntheta[seq_len(j-1)]))[j] + seq_len(ntheta[j])]
    if (dim(emat)[1] > 1) {
      emat[cbind(2:ntheta[j], 1:(ntheta[j]-1))] <- -1
      emat <- sweep(emat, 1, c(1, diff(theta_j)), "/")
    }
    emat
  })
}

tr_r_function <- function(rvec, ndim, i) {
  R <- diag(ndim)
  R[lower.tri(R)] <- rvec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  l <- t(chol(R))
  angmat <- diag(ndim)

  angmat[-1,1] <- acos(l[-1,1])
  if (ndim > 2){
    for (j in 2:(ndim-1)){
      sinprod <- apply(sin(angmat[, seq_len(j-1), drop=FALSE]), 1, prod) ## denominator in division
      angmat[-(1:j),j]<-acos((l/sinprod)[-(1:j),j])
    }
  }
  angdivpi <- angmat[lower.tri(angmat)]/pi

  log(angdivpi/(1-angdivpi))[i]
}

backtransf_sigmas <- function(R){
  J <- nrow(R)
  l <- t(chol(R))
  angmat <- matrix(1,ncol=J,nrow=J)
  angmat[-1,1] <- acos(l[-1,1])
  for (j in 2:(J-1)){
    sinprod <- apply(sin(angmat[, seq_len(j-1), drop=F]), 1, prod) ## denominator in division
    angmat[-(1:j),j]<-acos((l/sinprod)[-(1:j),j])
  }
  angdivpi <- angmat[lower.tri(angmat)]/pi
  log(angdivpi/(1-angdivpi))
}



jac_dtr_dr <- function(rvec, ndim){
 t(sapply(seq_along(rvec), function(i)
    numDeriv::grad(function(x) tr_r_function(x, ndim, i),
                                          x=rvec)))
}

get_labels_theta <- function(yj, nthetaj) {
  lev <- 1:length(unique(yj)) # TODO
  sapply(seq_len(nthetaj), function(i){
    paste(lev[i], lev[i + 1], sep = "|")
  })
}
