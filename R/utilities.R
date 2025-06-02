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

make_start_values <- function(y, X, response_types) {
  start_theta <- lapply(which(response_types == "ordinal"),
                        function(j) (ordinal::clm(factor(y[, j]) ~ 1)$alpha))
  ndim <- ncol(y)
  p <- ncol(X)
  ido <- which(response_types == "ordinal")
  pars <- c(unlist(lapply(start_theta, function(x) c(x[1], log(diff(x))))),
            # beta0 for normals
            colMeans(y[, -ido, drop = FALSE], na.rm = TRUE),
            # betas
            rep(0, ndim * p),
            # sigmas for normals
            rep(0, sum(response_types != "ordinal")),
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

jac_dtheta_dttheta_flexible <- function(ttheta, ndimo, ntheta) {
  lapply(seq_len(ndimo), function(j){
        emat <- diag(ntheta[j])
        if (ncol(emat) >= 2) {
              emat[,1] <- 1
              for (k in 2:ncol(emat))
                 emat[(k:nrow(emat)), k] <-
                  exp(ttheta[cumsum(c(0, ntheta[seq_len(j-1)]))[j] + seq_len(ntheta[j])])[k]
        }
        emat
  })

}

tr_r_function <- function(rvec, ndim) {
  R <- diag(ndim)
  R[lower.tri(R)] <- rvec
  R[upper.tri(R)] <- t(R)[upper.tri(R)]
  chR <- chol(R) # tryCatch(chol(R), error = function(e) {
  #  chol(Matrix::nearPD(R, corr = TRUE, keepDiag = TRUE)$mat)
  #})
  l <- t(chR)
  angmat <- diag(ndim)

  angmat[-1,1] <- acos(l[-1,1])
  if (ndim > 2){
    for (j in 2:(ndim-1)){
      sinprod <- apply(sin(angmat[, seq_len(j-1), drop=FALSE]), 1, prod) ## denominator in division
      angmat[-(1:j),j]<-acos((l/sinprod)[-(1:j),j])
    }
  }
  angdivpi <- angmat[lower.tri(angmat)]/pi

  log(angdivpi/(1-angdivpi))
}

r_tr_function <- function(tpar, ndim, i) {
  nu <- tpar
  angles <- pi * exp(nu)/(1 + exp(nu))
  cosmat <- diag(ndim)
  cosmat[lower.tri(cosmat)] <- cos(angles)
  S1 <- matrix(0, nrow = ndim, ncol = ndim)
  S1[lower.tri(S1, diag = TRUE)] <- c(rep(1, ndim), sin(angles))
  tLmat <- sapply(1:ndim,
                  function(j) cosmat[j, ] * cumprod(S1[j, ]))
  sigma <- crossprod(tLmat)
  sigma[lower.tri(sigma)][i]
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


jac_dr_dtr <- function(tpar, ndim){
  t(sapply(seq_along(tpar), function(i)
    numDeriv::grad(function(x) r_tr_function(x, ndim, i),
                   x=tpar)))
}

#
# jac_dtr_dr <- function(rvec, ndim){
#  t(sapply(seq_along(rvec), function(i)
#     numDeriv::grad(function(x) tr_r_function(x, ndim),
#                                           x=rvec)))
# }

get_labels_theta <- function(yj, nthetaj) {
  lev <- 1:length(unique(yj)) # TODO
  sapply(seq_len(nthetaj), function(i){
    paste(lev[i], lev[i + 1], sep = "|")
  })
}

f_transf <- function(pars, response_types,ntheta,ndim, ndimn,p) {
  thetas <- lapply(1:sum(response_types == "ordinal"), function(j) {
    transf_thresholds_flexible(pars[cumsum(c(0, ntheta))[j] + seq_len(ntheta[j])])
  })
  ## regression coefs: intercepts for normals
  beta0n <- pars[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
  ## common regression coefs
  beta <- pars[sum(ntheta) + ndimn + seq_len(ndim * p)]
  ## sd parameters for normals
  sigman <- exp(pars[sum(ntheta) + ndimn + ndim * p + seq_len(ndimn)])
  ## correlations error structure
  tparerror <- pars[(sum(ntheta) + ndimn + ndim * p + ndimn + seq_len(ndim * (ndim - 1)/2))]
  rvec <- transf_sigma(tparerror, ndim)
  ##
  return(c(unlist(thetas), beta0n, beta, sigman, rvec))
}


group_by_missingness <- function(mat) {
  # Create a matrix indicating missingness (TRUE for missing, FALSE for non-missing)
  missing_pattern <- is.na(mat)
  # Convert each row's missingness pattern to a unique character string
  pattern_labels <- apply(missing_pattern, 1, paste, collapse = "-")
  # Find unique missingness patterns
  unique_patterns <- unique(pattern_labels)
  n_unique_patterns <- length(unique_patterns)
  # Loop through each unique missingness pattern
  grouped_rows <- vector("list", n_unique_patterns)
  for (i in seq_len(n_unique_patterns)) {
    # Identify the row id for the pattern
    ind_i  <- which(pattern_labels == unique_patterns[i])
    # Identify the column id for the non missing pattern
    combis <- which(strsplit(unique_patterns[i], "-")[[1]] == "FALSE")
    # If all are missing, skip to next
    if (length(combis) == 0) next;
    # Save as a list
    grouped_rows[[i]] <- list(ind_i = ind_i,
                              combis = combis)
  }
  # Remove empty elements
  grouped_rows <- grouped_rows[!sapply(grouped_rows, is.null)]
  # Return
  return(grouped_rows)
}

#
# neg_log_likelihood_multivariate <- function(par, Y, X) {
#   # Y: n x m matrix of response variables (n observations, m dependent variables)
#   # X: n x p matrix of independent variables (design matrix, n observations, p predictors)
#   # B: p x m matrix of regression coefficients
#   # Sigma: m x m covariance matrix of the errors
#
#   # Number of observations (n) and number of response variables (m)
#   n <- nrow(Y)
#   m <- ncol(Y)
#   p <- ncol(X)
#   beta <- par[seq_len(m * p)]
#   dim(beta) <- c(p, m)
#   sigma <- (par[(m * p + seq_len(m))])
#   rvec <- par[(m * p + m + 1)]
#
#   smat <- diag(nrow = m)
#   smat[lower.tri(smat)] <- rvec
#   R <- smat + t(smat) - diag(m)
#   D <- diag(sigma, nrow = ncol(Y))
#   Sigma <- D %*% R %*% D
#
#   # Compute the residual matrix (Y - XB)
#   B <- beta
#   residuals <- as.matrix(Y - X %*% B)
#
#   # Inverse of the covariance matrix Sigma
#   Sigma_inv <- solve(Sigma)
#
#   # Log determinant of Sigma
#   log_det_Sigma <- log(det(Sigma))
#
#   # Log-likelihood computation
#   log_likelihood <- 0
#
#   for (i in 1:n) {
#     # Extract the i-th residual (row vector)
#     resid_i <- residuals[i, , drop = FALSE]
#
#     # Compute the quadratic form: (Y_i - X_i B)^T Sigma_inv (Y_i - X_i B)
#     quad_form <- resid_i %*% Sigma_inv %*% t(resid_i)
#
#     # Update log-likelihood
#     log_likelihood <- log_likelihood - 0.5 * (quad_form + log_det_Sigma + m * log(2 * pi))
#   }
#
#   return(-log_likelihood)
# }
#
# f <- function(pars) {
#   #thetas <- lapply(seq_len(ndimo), function(j) {
#   #  transf_thresholds_flexible(pars[cumsum(c(0, ntheta))[j] + seq_len(ntheta[j])])
#   #})
#   ## regression coefs: intercepts for normals
#   #beta0n <- pars[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
#   ## common regression coefs
#   #beta <- pars[sum(ntheta) + ndimn + seq_len(ndim * p)]
#   ## sd parameters for normals
#   # sigman <- exp(pars)#[sum(ntheta) + ndimn + ndim * p + seq_len(ndimn)])
#   ## correlations error structure
#   #pos_corrs <- (sum(ntheta) + 2 * ndimn + ndim * p + seq_len(ndim * (ndim - 1)/2))
#   #tparerror <- pars[pos_corrs]
#   #rvec <- transf_sigma(tparerror, ndim)
#   #smat <- diag(nrow = ndim)
#   #smat[lower.tri(smat)] <- rvec
#   #R <- smat + t(smat) - diag(ndim)
#
#   beta_mat <- beta
#   dim(beta_mat) <- c(ndim, p)
#   #print(beta_mat)
#   Xbeta <- tcrossprod(X, beta_mat)
#   XU <- lapply(1:ndimo, function(j) {
#     B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
#     BU <- B2[,-(ntheta[j] + 1), drop = FALSE]
#     cbind(BU, - X)
#   })
#   XL <- lapply(1:ndimo, function(j) {
#     B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
#     BL <- B2[,-1, drop = FALSE]
#     cbind(BL, - X)
#   })
#
#   ## X matrix for normal variables to include intercept
#   Xn <- cbind(1, X)
#
#   eta_u_o <- sapply(1:ndimo, function(j)
#     c(thetas[[j]], 1e06)[y[, ido[j]]] - Xbeta[, ido[j]])
#
#   eta_l_o <- sapply(1:ndimo, function(j)
#     drop(c(-1e06, thetas[[j]])[y[, ido[j]]] - Xbeta[, ido[j]]))
#
#   eta_n <- sapply(1:ndimn, function(j) {
#     beta0n[j] + Xbeta[, idn[j]]
#   })
#
#   if (is.null(dim(eta_l_o))) {
#     dim(eta_l_o) <- dim(eta_u_o) <- c(1, ndimo)
#     dim(eta_n) <- c(1, ndimn)
#   }
#   log_pl_n_cond <- sapply(combis_fast$normal[2], function(x) {
#     # CASE 2: all normals ----
#     idx <- idn[x$combis]   # active cols in y in na pattern
#     id <- x$ind_i          # row indices for the na pattern
#     sx <- sigman[x$combis] # sigmas for these cols
#     Rx <- R[idx, idx]      # corresponding correlation matrix
#     chRu <- chol(Rx)       # Cholesky upper
#     D <- diag(x = sx, nrow = length(sx))
#     chS <- tcrossprod(D, chRu)
#     eta_n_x <- eta_n[id, x$combis, drop = FALSE]
#     y_x     <- y[id, , drop = FALSE]
#     Xbeta_x <- Xbeta[id, , drop = FALSE]
#     C <- ltMatrices(chS[lower.tri(chS, diag = TRUE)], diag = TRUE)
#     log_pl_n <- sum(ldmvnorm(t(y_x[, idx]), t(eta_n_x),
#                              logLik = TRUE, chol = C))
#
#     # CASE 3: each ordinal conditional on all normals ----
#     eps <- y_x[, idx] - eta_n_x
#     eps <- sweep(eps, 2, sx, "/")
#     Rxinv <- chol2inv(chRu)
#     log_pl_cond <- sapply(1:2, function(j) {
#       if (j == 1) rho <- pars[1:2] else rho <- pars[3:4]
#       #rho <- pars[idx-2] # R[j, c(idx)]
#       eta_cond <- Xbeta_x[,j] + tcrossprod(as.matrix(eps), rho %*% Rxinv)
#       sd_cond  <- drop(sqrt(1 - rho %*% Rxinv %*% rho))
#       eta_u_cond <- (c(thetas[[j]],  1e06)[y_x[, j]] - eta_cond)/sd_cond
#       eta_l_cond <- (c(-1e06, thetas[[j]])[y_x[, j]] - eta_cond)/sd_cond
#       p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)
#       p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps
#       sum(log(p_cond), na.rm = TRUE)
#     })
#     - sum(log_pl_cond)
#   })
#   sum(log_pl_n_cond)
# }
# pars2<-pars[19:20]
# par3 <- 0.8961249
#pars2<-c(0.7824650, 0.6576646,0.9199078,0.8149329)
# #pars2[pos_s] <- exp(pars2[pos_s])
#dnum2<-numDeriv::grad(function(p) f(p),
#            x = pars2)
#dnum2
# dnum2[pos_rho]
# 0.7831046 0.6588395
# par2 <- c(pars[c( 5,  9, 13, 17,  6, 10, 14, 18)], exp(pars[19:20]), 0.8961249)
# numDeriv::grad(function(p) neg_log_likelihood_multivariate(p, y[,3:4], Xn),
#                x = par2)



