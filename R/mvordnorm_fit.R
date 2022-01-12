mvordnorm_fit <- function(y, X, # w,  offset,
                         response_types,  control) {

  ndimo <- sum(response_types == "ordinal")
  ndim <- ncol(y)
  ndimn <- ndim  - ndimo
  n <- nrow(y)
  p <- ncol(X)
  idn <- which(response_types != "ordinal")
  ido <- which(response_types == "ordinal")
  ## number of thresholds
  ntheta <- apply(y[,response_types=="ordinal", drop = FALSE], 2,
                  function(x) nlevels(as.factor(x)) - 1)
  combis <- combn(ndim, 2)
  obj <- list()

  start_values <- make_start_values(y, X, family_type = response_types)
  ## Optimize negative log likelihood
  obj$res <- optimx(start_values, function(par)
    neg_log_lik_joint(par, response_types,
                      y, X, ntheta, p, ndimo, ndimn, ndim, combis,
                      idn, ido),
                    method = control$solver)
  ## TODO
  obj$parOpt <- unlist(obj$res[seq_along(start_values)])

  ## Finalize
  ynames <- names(y)

  tpar   <- unlist(obj$parOpt)
  thetas <- unlist(lapply(1:sum(response_types == "ordinal"), function(j) {
    transf_thresholds_flexible(tpar[cumsum(c(0, ntheta))[j] + 1:ntheta[j]])
  }))
  names(thetas) <- unlist(lapply(1:ndimo, function(j)
    paste(ynames[ido[j]], get_labels_theta(y[, ido[j]], ntheta[j]))))
  ## regression coefs: intercepts for normals
  beta0n <- tpar[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
  names(beta0n)   <- paste0("beta0.", ynames[idn])
  ## common regression coefs
  beta <- tpar[sum(ntheta) + ndimn + seq_len(ndim * p)]

  names(beta)  <-  c(outer(ynames, colnames(X), paste0))
  ## sd parameters for normals
  tsigman <- tpar[sum(ntheta) + ndimn + ndim * p + seq_len(ndimn)]
  sigman <- exp(tsigman)
  names(sigman)   <- paste0("sigma.", ynames[idn])
  ## correlations error structure
  tparerror <- tpar[(sum(ntheta) + ndimn + ndim * p + ndimn + seq_len(ndim * (ndim - 1)/2))]
  rvec <- transf_sigma(tparerror, ndim)
  names(rvec) <-
    apply(combn(ndim,2), 2,function(x)
    sprintf("corr_%s_%s", ynames[x[1]], ynames[x[2]]))


  obj$parameters <- list(thetas, beta0n, beta,  sigman, rvec)

  ## Standard errors
  if (control$se) {
    cat("Computing variability and hessian matrix numerically.\n")
    ## Compute Hessian numerically
    tparHess <- numDeriv::hessian(function(par)
      neg_log_lik_joint(par, response_types, y, X, ntheta,
                        p, ndimo, ndimn, ndim, combis, idn, ido),
      obj$parOpt)

    jac_list <- jac_dttheta_dtheta_flexible(thetas,
                                        ndimo, ntheta) ## d ttheta/d theta
    jac_list[ndimo + seq_len(ndimn)] <- 1 ## d tbeta0/dbeta0
    jac_list[ndim + seq_len(ndim * p)] <- 1 ## d tbeta/dbeta
    jac_list[[ndim + ndim * p + 1]] <- diag(1/sigman) ## d tsigma/dsigma
    jac_list[[ndim + ndim * p + 2]]<- jac_dtr_dr(rvec, ndim) # d tr/dr

    J <- as.matrix(Matrix::bdiag(jac_list))
    J.inv <- solve(J)
    Vi_num <- matrix(0, ncol = length(tpar), nrow = n)
    H <- crossprod(J, tparHess) %*% J
    for (i in seq_len(n)) {
         if (i %% 100 == 0)  cat('Computed gradient for', i, 'out of', n,'subjects\n')
         Vi_num[i, ] <- numDeriv::grad(function(par)
           neg_log_lik_joint(par, response_types,
                               y[i, ],
                               X[i, ],
                               ntheta, p, ndimo, ndimn, ndim,
                               combis, idn, ido), obj$parOpt,
           method = "Richardson")
    }
    Vi_num <-  Vi_num %*% J.inv
    V <- n/(n - length(tpar)) * crossprod(Vi_num)
    H.inv <- tryCatch(chol2inv(chol(H)),
                         error=function(e) {
                           warning("\nCondition number close to zero! Hessian is approximated by nearest positive semidefinite matrix.\n")
                           chol2inv(chol(Matrix::nearPD(H)$mat))
                         }
    )
    obj$vcov <- H.inv %*% V %*% H.inv
  }

  obj
}

neg_log_lik_joint <- function(pars, response_types, y, X,
                              ntheta, p, ndimo, ndimn, ndim,
                              combis, idn, ido) {

  ## thresholds
  thetas <- lapply(1:sum(response_types == "ordinal"), function(j) {
    transf_thresholds_flexible(pars[cumsum(c(0, ntheta))[j] + 1:ntheta[j]])
  })


  ## regression coefs: intercepts for normals
  beta0n <- pars[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
  ## coomon regression coefs
  beta <- pars[sum(ntheta) + ndimn + seq_len(ndim * p)]
  ## sd parameters for normals
  sigman <- exp(pars[sum(ntheta) + ndimn + ndim * p + seq_len(ndimn)])
  ## correlations error structure
  tparerror <- pars[(sum(ntheta) + ndimn + ndim * p + ndimn + seq_len(ndim * (ndim - 1)/2))]
  rvec <- transf_sigma(tparerror, ndim)

  beta_mat <- beta
  dim(beta_mat) <- c(ndim, p)

  Xbeta <- tcrossprod(X, beta_mat)

  eta_u_o <- lapply(1:ndimo, function(j)
      c(thetas[[j]], 1e06)[y[, ido[j]]] - Xbeta[, ido[j]])

  eta_l_o <- lapply(1:ndimo, function(j)
      drop(c(-1e06, thetas[[j]])[y[, ido[j]]] - Xbeta[, ido[j]]))

  eta_n <- lapply(1:ndimn, function(j) {
    beta0n[j] + Xbeta[, idn[j]]
  } )

  log_pl_vec <- NULL

  for (s in 1:ncol(combis)) {
    k <- combis[1, s]
    l <- combis[2, s]
    r <- rvec[s]
    ## CASE 1: 2 ordinals
    if (all(response_types[c(k,l)] == "ordinal")) {
      kido <- which(ido == k)
      lido <- which(ido == l)

      prs <- rectbiv_norm_prob(U = eta_u_o[c(kido, lido)],
                               L = eta_l_o[c(kido, lido)], r)
      prs[prs < .Machine$double.eps] <- .Machine$double.eps
      log_pl_vec[s] <- sum(log(prs))
    }

    ## CASE 2: 2 normals
    if (all(response_types[c(k,l)] != "ordinal")) {
      kidn <- which(idn == k)
      lidn <- which(idn == l)
      ykstd <-  (y[,k] - eta_n[[kidn]])
      ylstd  <- (y[,l] - eta_n[[lidn]])
      smat <- diag(2)
      smat[1,2] <- smat[2, 1] <- r
      smat <- tcrossprod(sigman[c(kidn, lidn)]) * smat
      log_pl_vec[s] <-  sum(mvtnorm::dmvnorm(x = cbind(ykstd, ylstd),
                                             mean = rep(0, 2),
                                             sigma = smat,
                                             log = TRUE))
    }

    ## CASE 3: ordinal + normal
    if (response_types[k] != response_types[l]) {
      ko <- c(k,l)[which(response_types[c(k,l)] == "ordinal")]
      kn <- c(k,l)[which(response_types[c(k,l)] != "ordinal")]
      knidn <- which(idn == kn)
      koido <- which(ido == ko)
      ld_marg <- dnorm(y[,kn],
                       mean = eta_n[[knidn]],
                       sd = sigman[knidn],
                       log = TRUE)

      sd_yn_cond <- sqrt(1 - r^2)

      eta_cond <-  Xbeta[,kn] + r * (y[,kn] - eta_n[[knidn]]) / sigman[knidn]

      eta_u_cond <- (c(thetas[[koido]],  1e06)[y[,ko]] - eta_cond)/sd_yn_cond

      eta_l_cond <- (c(-1e06, thetas[[koido]])[y[,ko]] - eta_cond)/sd_yn_cond

      p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)

      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps

      log_pl_vec[s] <-  sum(log(p_cond) + ld_marg)
    }
  }

  - sum(log_pl_vec)
}

