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



  ## Missings in y
  ind_univ <- which(!is.na(y) & rowSums(!is.na(y)) == 1, arr.ind = TRUE)
  n_univ <- NROW(ind_univ)
  ## index for subjects containing pair c(k,l)
  ind_kl <- lapply(1:ncol(combis), function(j) {
    kl <- combis[,j]
    rowSums(!is.na(y[, kl])) == 2
  })

  ind_i <- lapply(1:ncol(combis), function(h)  which(ind_kl[[h]]))

  combis_fast <- lapply(1:ncol(combis), function(h){
    list("combis" = combis[,h],
         "ind_i" =  ind_i[[h]],
         "rpos" = h)
  })

  ## Optimize negative log likelihood
  obj$res <- optimx(start_values, function(par)
    neg_log_lik_joint(par, response_types,
                      y, X, ntheta, p, ndimo, ndimn, ndim,
                      idn, ido, ind_univ, combis_fast),
    method = control$solver)
  ## TODO
  obj$parOpt <- unlist(obj$res[seq_along(start_values)])
  obj$combis_fast <- combis_fast

  ## Finalize
  ynames <- names(y)

  tpar   <- unlist(obj$parOpt)
  tparTheta <- tpar[1:sum(ntheta)]
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
    res_deriv_ana <- derivs_ana(obj$parOpt, y, X, response_types,
                                ind_univ, combis_fast)
    V <- n/(n - NCOL(res_deriv_ana$V)) * res_deriv_ana$V  ## correct for degrees of freedom
    H.inv <- tryCatch(chol2inv(chol(res_deriv_ana$H)), error=function(e) NA)
    if (length(H.inv) == 1 && is.na(H.inv)) {
      warning("Hessian is not positive semi-definite. Approximating to nearest PD matrix.")
      Happrox <- Matrix::nearPD(res_deriv_ana$H)$mat
      H.inv <- chol2inv(chol(Happrox))
    }
    obj$vcov <- H.inv %*% V %*% H.inv
    # cat("Computing variability and hessian matrix numerically.\n")
    # ## Compute Hessian numerically
    # tparHess <- numDeriv::hessian(function(par)
    #   neg_log_lik_joint(par, response_types, y, X, ntheta,
    #                     p, ndimo, ndimn, ndim, idn, ido,
    #                     ind_univ, combis_fast),
    #   obj$parOpt)
    #
    # # diag(solve(tparHess))
    #
    # cat("Hessian: Done.\n")
    # jac_list <- jac_dttheta_dtheta_flexible(thetas,
    #                                          ndimo, ntheta) ## d ttheta/d theta
    # jac_list[ndimo + seq_len(ndimn)] <- 1 ## d tbeta0/dbeta0
    # jac_list[ndim + seq_len(ndim * p)] <- 1 ## d tbeta/dbeta
    # jac_list[[ndim + ndim * p + 1]] <- diag(ndimn, x = 1/sigman) ## d tsigma/dsigma
    # jac_list[[ndim + ndim * p + 2]]<- jac_dtr_dr(rvec, ndim) # d tr/dr
    #
    # jac_list_inv <- jac_dtheta_dttheta_flexible(tparTheta,
    #                                         ndimo, ntheta) ## d theta/d ttheta
    # jac_list_inv[ndimo + seq_len(ndimn)] <- 1 ## d beta0/dtbeta0
    # jac_list_inv[ndim + seq_len(ndim * p)] <- 1 ## d ttbeta/dtbeta
    # jac_list_inv[[ndim + ndim * p + 1]] <- diag(ndimn, x = sigman) ## d sigma/dtsigma
    # jac_list_inv[[ndim + ndim * p + 2]]<- jac_dr_dtr(tparerror, ndim) # d tr/dr
    #
    #
    # J <- as.matrix(Matrix::bdiag(jac_list))
    # J.inv <- as.matrix(Matrix::bdiag(jac_list_inv))# solve(J)
    #
    # H <- crossprod(J, tparHess) %*% J
    # Vi_num <- matrix(0, ncol = length(tpar), nrow = n)
    #
    # # u <- numDeriv::grad(function(par)
    # #   neg_log_lik_joint(par, response_types,
    # #                     y, X,
    # #                     ntheta, p, ndimo, ndimn, ndim,
    # #                     idn, ido,
    # #                     ind_univ = ind_univ,
    # #                     combis_fast = combis_fast),
    # #   obj$parOpt,
    # #   method = "Richardson") # score
    # #
    # # V <- n/(n - length(tpar)) * crossprod(crossprod(u, J.inv))
    # for (i in seq_len(n)) {
    #   if (i %% 100 == 0)  cat('Computed gradient for', i, 'out of', n,'subjects\n')
    #   y_pos_nona <- which(!is.na(y[i,]))
    #   if (length(y_pos_nona) == 1) {
    #     ind_univ_i <- matrix(c(1, y_pos_nona), nrow = 1, ncol = 2)
    #     combis_fast_i <- NULL
    #   } else {
    #     ind_univ_i <- matrix(nrow = 0, ncol = 2)
    #     comb_i_tmp <- combn(y_pos_nona,  2, simplify = FALSE)
    #     combis_fast_i <- lapply(comb_i_tmp, function(x)
    #       list(combis = x,
    #            ind_i = 1,
    #            rpos = which(apply(combis, 2, function(k) all(x==k)))))
    #   }
    #   Vi_num[i, ] <- numDeriv::grad(function(par)
    #     neg_log_lik_joint(par, response_types,
    #                       y[i, ], X[i, ],
    #                       ntheta, p, ndimo, ndimn, ndim,
    #                       idn, ido,
    #                       ind_univ = ind_univ_i,
    #                       combis_fast = combis_fast_i),
    #     obj$parOpt,
    #     method = "simple")
    # }
    #
    # cat("Variability matrix: Done.\n")
    # Vi_num <-  Vi_num %*% J.inv
    # V <- n/(n - length(tpar)) * crossprod(Vi_num)
    # H.inv <- tryCatch(chol2inv(chol(H)),
    #                   error=function(e) {
    #                     warning("\nCondition number close to zero! Hessian is approximated by nearest positive semidefinite matrix.\n")
    #                     chol2inv(chol(Matrix::nearPD(H)$mat))
    #                   }
    # )
    #obj$vcov <- H.inv %*% V %*% H.inv
  }

  obj
}

neg_log_lik_joint <- function(pars, response_types, y, X,
                              ntheta, p, ndimo, ndimn, ndim,
                              idn, ido,
                              ind_univ,
                              combis_fast) {

  ## thresholds
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

  beta_mat <- beta
  dim(beta_mat) <- c(ndim, p)
  #print(beta_mat)
  Xbeta <- tcrossprod(X, beta_mat)

  eta_u_o <- sapply(1:ndimo, function(j)
    c(thetas[[j]], 1e06)[y[, ido[j]]] - Xbeta[, ido[j]])

  eta_l_o <- sapply(1:ndimo, function(j)
    drop(c(-1e06, thetas[[j]])[y[, ido[j]]] - Xbeta[, ido[j]]))

  eta_n <- sapply(1:ndimn, function(j) {
    beta0n[j] + Xbeta[, idn[j]]
  })

  if (is.null(dim(eta_l_o))) {
    dim(eta_l_o) <- dim(eta_u_o) <- c(1, ndimo)
    dim(eta_n) <- c(1, ndimn)
  }
  log_pl_vec_h <- unlist(sapply(combis_fast, function(x) {
    k <- x$combis[1]
    l <- x$combis[2]
    id <- x$ind_i

    r <- rvec[x$rpos]
    ## CASE 1: 2 ordinals
    if (all(response_types[c(k,l)] == "ordinal")) {
      kido <- which(ido == k)
      lido <- which(ido == l)
      prs <- rectbiv_norm_prob(U = eta_u_o[id, c(kido, lido), drop = FALSE],
                               L = eta_l_o[id, c(kido, lido), drop = FALSE], r)

      prs[prs < .Machine$double.eps] <- .Machine$double.eps

      log_pl_vec <- sum(log(prs))
    }

    ## CASE 2: 2 normals
    if (all(response_types[c(k,l)] != "ordinal")) {
      kidn <- which(idn == k)
      lidn <- which(idn == l)
      ykstd <-  (y[id,k] - eta_n[id, kidn, drop = FALSE])
      ylstd  <- (y[id,l] - eta_n[id, lidn, drop = FALSE])
      smat <- diag(2)
      smat[1,2] <- smat[2, 1] <- r
      smat <- tcrossprod(sigman[c(kidn, lidn)]) * smat
      log_pl_vec <-  sum(mvtnorm::dmvnorm(x = cbind(ykstd, ylstd),
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
      ld_marg <- dnorm(y[id,kn],
                       mean = eta_n[id, knidn,drop = FALSE],
                       sd = sigman[knidn],
                       log = TRUE)

      sd_yn_cond <- sqrt(1 - r^2)

      eta_cond <-  Xbeta[id, ko] +
        r * (y[id,kn] - eta_n[id, knidn, drop = FALSE]) / sigman[knidn]

      eta_u_cond <- (c(thetas[[koido]],  1e06)[y[id,ko]] - eta_cond)/sd_yn_cond
      eta_l_cond <- (c(-1e06, thetas[[koido]])[y[id,ko]] - eta_cond)/sd_yn_cond

      p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)
      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps
      log_pl_vec <-  sum(log(p_cond) + ld_marg)
    }
    log_pl_vec

  }))
  # This are the observations which have only one response (aka univariate)
  colu <- ind_univ[, 2] # col of univariate
  rowu <- ind_univ[, 1] # row of univariate
  idn_univ <- cbind(rowu, match(colu, idn))
  ido_univ <- cbind(rowu,  match(colu, ido))
  ndens <- dnorm(y[cbind(rowu, colu)],
                 eta_n[idn_univ],
                 sigman[match(colu, idn)], log =TRUE)
  lndens <-  ifelse(is.na(ndens), 0, ndens)
  odens  <- log(pnorm(eta_u_o[ido_univ]) - pnorm(eta_l_o[ido_univ]))
  lodens <-  ifelse(is.na(odens), 0, odens)
  univ_nll <- lndens + lodens
  ## Final neg log lik
  - sum(log_pl_vec_h) - sum(univ_nll)
}

#
# neg_log_lik_joint_i <- function(pars, response_types, yi, Xi,
#                               ntheta, p, ndimo, ndimn, ndim,
#                               idn, ido,
#                               ind_univ,
#                               combis_fast) {
#
#   ## thresholds
#   thetas <- lapply(1:sum(response_types == "ordinal"), function(j) {
#     transf_thresholds_flexible(pars[cumsum(c(0, ntheta))[j] + seq_len(ntheta[j])])
#   })
#
#
#   ## regression coefs: intercepts for normals
#   beta0n <- pars[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
#   ## common regression coefs
#   beta <- pars[sum(ntheta) + ndimn + seq_len(ndim * p)]
#   ## sd parameters for normals
#   sigman <- exp(pars[sum(ntheta) + ndimn + ndim * p + seq_len(ndimn)])
#   ## correlations error structure
#   tparerror <- pars[(sum(ntheta) + ndimn + ndim * p + ndimn + seq_len(ndim * (ndim - 1)/2))]
#   rvec <- transf_sigma(tparerror, ndim)
#
#   beta_mat <- beta
#   dim(beta_mat) <- c(ndim, p)
#   #print(beta_mat)
#   Xbeta <- tcrossprod(X, beta_mat)
#
#   eta_u_o <- sapply(1:ndimo, function(j)
#     c(thetas[[j]], 1e06)[yi[ido[j]]] - Xbeta[,ido[j]])
#
#   eta_l_o <- sapply(1:ndimo, function(j)
#     drop(c(-1e06, thetas[[j]])[yi[ido[j]]] - Xbeta[,ido[j]]))
#
#   eta_n <- sapply(1:ndimn, function(j) {
#     beta0n[j] + Xbeta[,idn[j]]
#   })
#
#   log_pl_vec_h <- unlist(sapply(combis_fast, function(x) {
#     k <- x$combis[1]
#     l <- x$combis[2]
#     id <- x$ind_i
#
#     r <- rvec[x$rpos]
#     ## CASE 1: 2 ordinals
#     if (all(response_types[c(k,l)] == "ordinal")) {
#       kido <- which(ido == k)
#       lido <- which(ido == l)
#       prs <- rectbiv_norm_prob(U = eta_u_o[ c(kido, lido), drop = FALSE],
#                                L = eta_l_o[ c(kido, lido), drop = FALSE], r)
#
#       prs[prs < .Machine$double.eps] <- .Machine$double.eps
#
#       log_pl_vec <- (log(prs))
#     }
#
#     ## CASE 2: 2 normals
#     if (all(response_types[c(k,l)] != "ordinal")) {
#       kidn <- which(idn == k)
#       lidn <- which(idn == l)
#       ykstd <-  (y[k] - eta_n[kidn, drop = FALSE])
#       ylstd  <- (y[l] - eta_n[lidn, drop = FALSE])
#       smat <- diag(2)
#       smat[1,2] <- smat[2, 1] <- r
#       smat <- tcrossprod(sigman[c(kidn, lidn)]) * smat
#       log_pl_vec <-  sum(mvtnorm::dmvnorm(x = c(ykstd, ylstd),
#                                           mean = rep(0, 2),
#                                           sigma = smat,
#                                           log = TRUE))
#     }
#
#     ## CASE 3: ordinal + normal
#     if (response_types[k] != response_types[l]) {
#       ko <- c(k,l)[which(response_types[c(k,l)] == "ordinal")]
#       kn <- c(k,l)[which(response_types[c(k,l)] != "ordinal")]
#       knidn <- which(idn == kn)
#       koido <- which(ido == ko)
#       ld_marg <- dnorm(y[kn],
#                        mean = eta_n[knidn],
#                        sd = sigman[knidn],
#                        log = TRUE)
#
#       sd_yn_cond <- sqrt(1 - r^2)
#
#       eta_cond <-  Xbeta[, ko] +
#         r * (y[kn] - eta_n[knidn, drop = FALSE]) / sigman[knidn]
#
#       eta_u_cond <- (c(thetas[[koido]],  1e06)[y[ko]] - eta_cond)/sd_yn_cond
#       eta_l_cond <- (c(-1e06, thetas[[koido]])[y[ko]] - eta_cond)/sd_yn_cond
#
#       p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)
#       p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps
#       log_pl_vec <-  (log(p_cond) + ld_marg)
#     }
#     log_pl_vec
#
#   }))
#   # This are the observations which have only one response (aka univariate)
#   colu <- ind_univ[, 2] # col of univariate
#   rowu <- ind_univ[, 1] # row of univariate
#   idn_univ <- cbind(rowu, match(colu, idn))
#   ido_univ <- cbind(rowu,  match(colu, ido))
#   ndens <- dnorm(y[cbind(rowu, colu)],
#                  eta_n[idn_univ],
#                  sigman[match(colu, idn)], log =TRUE)
#   lndens <-  ifelse(is.na(ndens), 0, ndens)
#   odens  <- log(pnorm(eta_u_o[ido_univ]) - pnorm(eta_l_o[ido_univ]))
#   lodens <-  ifelse(is.na(odens), 0, odens)
#   univ_nll <- lndens + lodens
#   ## Final neg log lik
#   - sum(log_pl_vec_h) - sum(univ_nll)
# }
