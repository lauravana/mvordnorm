mvordnorm_fit <- function(y, X, # w,  offset,
                          response_types,  control) {
  # Setup ----
  ndimo <- sum(response_types == "ordinal")
  ndim <- ncol(y)
  ndimn <- ndim  - ndimo
  n <- nrow(y)
  p <- ncol(X)
  idn <- which(response_types != "ordinal")
  ido <- which(response_types == "ordinal")

  ntheta <- apply(y[,response_types == "ordinal", drop = FALSE], 2,
                  function(x) nlevels(as.factor(x)) - 1)   # number of thresholds

  start_values <- make_start_values(y, X, response_types = response_types)

  # Missings in y ----
  if (control$type_composite_log_lik == "type_1") {
    ## Univariate ones ----
    ind_univ <- which(!is.na(y) & rowSums(!is.na(y)) == 1, arr.ind = TRUE)
    n_univ <- NROW(ind_univ)

    ## Bivariate ones ----
    combis <- combn(ndim, 2) # bivariate pairs
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
  }
  if (control$type_composite_log_lik == "type_2") {
    ### ordinal univariate
    yo <- y[, ido]
    ind_univ <- which(!is.na(yo) & rowSums(!is.na(yo)) == 1, arr.ind = TRUE)
    n_univ <- NROW(ind_univ)

    ## ordinal pairs
    combis <- combn(which(response_types == "ordinal"), 2)
    ind_kl <- lapply(1:ncol(combis), function(j) {
      kl <- combis[,j]
      rowSums(!is.na(y[, kl])) == 2
    })
    ind_i <- lapply(1:ncol(combis), function(h)  which(ind_kl[[h]]))
    combis_fast_o <- lapply(1:ncol(combis), function(h){
      list("combis" = combis[,h],
           "ind_i" =  ind_i[[h]],
           "rpos" = h)})

    ## Normals - find na patterns
    na_patterns <- group_by_missingness(y[, idn])
    combis_fast <- c(ordinal = list(combis_fast_o),
                     normal = list(na_patterns))
  }

  # Optimize negative log likelihood ----
  neg_log_lik_joint <- switch(control$type_composite_log_lik,
                              "type_1" = neg_log_lik_joint_type_1,
                              "type_2" = neg_log_lik_joint_type_2,
                              stop("Provided type of composite likelihood not supported.")
  )
  grad_neg_log_lik_joint <- switch(control$type_composite_log_lik,
                                   "type_1" = grad_neg_log_lik_joint_type_1,
                                   "type_2" = grad_neg_log_lik_joint_type_2,
                                   stop("Provided type of composite likelihood not supported.")
  )
  if (control$usegrfun) {
    grfun <- function(par) grad_neg_log_lik_joint(par, response_types,
                                                  y, X, ntheta, p, ndimo, ndimn, ndim,
                                                  idn, ido, ind_univ, combis_fast)
  } else {
    grfun <- NULL
  }


  obj <- list()
  obj$res <- optimx(start_values, function(par) neg_log_lik_joint(par, response_types,
                      y, X, ntheta, p, ndimo, ndimn, ndim,
                      idn, ido, ind_univ, combis_fast),
    gr = grfun,
    method = control$solver)

  # grfun(par = fit2$parOpt)
  # dnum<-numDeriv::grad(function(par) neg_log_lik_joint(par, response_types,
  #                                               y, X, ntheta, p, ndimo, ndimn, ndim,
  #                                               idn, ido, ind_univ, combis_fast),
  #                x = fit1$parOpt)
  # cbind(grfun(par = fit1$parOpt), dnum)
  # microbenchmark::microbenchmark( obj$res <- optimx(start_values, function(par)
  #   neg_log_lik_joint(par, response_types,
  #                     y, X, ntheta, p, ndimo, ndimn, ndim,
  #                     idn, ido, ind_univ, combis_fast),
  #   method = control$solver),
  #
  #   obj$res_2 <- optimx(start_values, function(par)
  #     neg_log_lik_joint(par, response_types,
  #                       y, X, ntheta, p, ndimo, ndimn, ndim,
  #                       idn, ido, ind_univ, combis_fast),
  #     gr = function(par)
  #       grad_neg_log_lik_joint(par, response_types,
  #                              y, X, Xn,ntheta, p, ndimo, ndimn, ndim,
  #                              idn, ido, ind_univ, combis_fast),
  #     method = control$solver), times = 1L)

  obj$objective <- obj$res[["value"]]
  obj$parOpt <- unlist(obj$res[seq_along(start_values)])
  obj$combis_fast <- combis_fast

  # Finalize ----
  ynames <- names(y)
  tpar   <- unlist(obj$parOpt)
  tparTheta <- tpar[1:sum(ntheta)]
  thetas <- unlist(lapply(seq_len(ndimo), function(j) {
    transf_thresholds_flexible(tpar[cumsum(c(0, ntheta))[j] + 1:ntheta[j]])
  }))
  names(thetas) <- unlist(lapply(seq_len(ndimo), function(j)
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

  names_all <- c(names(thetas), names(beta0n), names(beta), names(sigman),
                 names(rvec))
  obj$parameters <- list(thetas, beta0n, beta,  sigman, rvec)


  ## Standard errors ----
  if (control$se) {
    rowwise_pairwise_grad_neg_log_lik_joint <- switch(
      control$type_composite_log_lik,
      "type_1" = rowwise_pairwise_grad_neg_log_lik_joint_type_1,
      "type_2" = rowwise_pairwise_grad_neg_log_lik_joint_type_2,
      stop("Provided type of composite likelihood not supported.")

    )
    gradients_row_pairs_orig <-
      rowwise_pairwise_grad_neg_log_lik_joint(obj$parOpt,
                                              response_types,
                                              y, X,
                                              ntheta, p, ndimo,
                                              ndimn, ndim,
                                              idn, ido, ind_univ,
                                              combis_fast)
    ## Gradients of original parameters to enter optimizer
    ## Correct for transformations in theta, sd and r:
    jthetas <-
      jac_dttheta_dtheta_flexible(thetas,
                                  ndimo = ndimo, ntheta = ntheta)
    Jsd <- 1/sigman
    Jcor <- jac_dtr_dr(rvec, ndim)

    jacobian_neg_log_lik <- as.matrix(Matrix::bdiag(
      c(jthetas,
        rep(list(1), ndimn + ndim * p),
        Jsd,
        list(Jcor))
    ))

    gradients_row_pairs <- lapply(gradients_row_pairs_orig,
                                  function(x) x %*% jacobian_neg_log_lik)

    ## Variability matrix
    Vi <- Reduce("+", gradients_row_pairs)
    V <- crossprod(Vi)
    ## Hessian matrix
    H <- Reduce("+", lapply(gradients_row_pairs, crossprod))
    V <- n/(n - NCOL(V)) * V  ## correct for degrees of freedom
    H.inv <- tryCatch(chol2inv(chol(H)), error=function(e) NA)
    if (length(H.inv) == 1 && is.na(H.inv)) {
      warning("Hessian is not positive semi-definite. Approximating to nearest PD matrix.")
      Happrox <- Matrix::nearPD(H)$mat
      H.inv <- chol2inv(chol(Happrox))
    }
    obj$H.inv <-  H.inv
    obj$V <- V
    rownames(obj$V) <- colnames(obj$V) <-
      rownames(obj$H.inv) <- colnames(obj$H.inv) <-
      names_all
    obj$claic <- 2 * obj$res[["value"]] + 2 * sum(diag(V %*% H.inv))
    obj$clbic <- 2 * obj$res[["value"]] + log(n) * sum(diag(V %*% H.inv))
  }
  obj
}

neg_log_lik_joint_type_1 <- function(pars, response_types, y, X,
                                     ntheta, p, ndimo, ndimn, ndim,
                                     idn, ido,
                                     ind_univ,
                                     combis_fast) {

  ## thresholds
  thetas <- lapply(seq_len(ndimo), function(j) {
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
      prs[prs > 1] <- 1 - .Machine$double.eps

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

  # These are the observations which have only one response (aka univariate)
  colu <- ind_univ[, 2] # col of univariate
  rowu <- ind_univ[, 1] # row of univariate
  idn_univ <- cbind(rowu, match(colu, idn))
  ido_univ <- cbind(rowu,  match(colu, ido))
  ndens <- dnorm(y[ind_univ],
                 eta_n[idn_univ],
                 exp(sigman[match(colu, idn)]), log =TRUE)
  lndens <-  ifelse(is.na(ndens), 0, ndens)
  odens  <- log(pnorm(eta_u_o[ido_univ]) - pnorm(eta_l_o[ido_univ]))
  lodens <-  ifelse(is.na(odens), 0, odens)
  univ_nll <- lndens + lodens
  ## Final neg log lik
  - sum(log_pl_vec_h) - sum(univ_nll)
}

#pars <- fit2$parOpt
neg_log_lik_joint_type_2 <- function(pars, response_types, y, X,
                                     ntheta, p, ndimo, ndimn, ndim,
                                     idn, ido,
                                     ind_univ,
                                     combis_fast) {
  ## thresholds
  thetas <- lapply(seq_len(ndimo), function(j) {
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
  smat <- diag(nrow = ndim)
  smat[lower.tri(smat)] <- rvec
  R <- smat + t(smat) - diag(ndim)

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

  # CASE 1: ordinals ----
  # Case 1.a. univariate ordinal
  odens  <- log(pnorm(eta_u_o[ind_univ]) - pnorm(eta_l_o[ind_univ]))
  # Case 1.b. univariate ordinal
  log_pl_o <- sapply(combis_fast$ordinal, function(x) {
    k <- x$combis[1]
    l <- x$combis[2]
    r <- rvec[x$rpos]
    id <- x$ind_i
    kido <- which(ido == k)
    lido <- which(ido == l)
    prs <- rectbiv_norm_prob(U = eta_u_o[id, c(kido, lido), drop = FALSE],
                             L = eta_l_o[id, c(kido, lido), drop = FALSE],
                             r)

    prs[prs < .Machine$double.eps] <- .Machine$double.eps

    sum(log(prs))
  })

  # Iterate over all na patterns in normals
  log_pl_n_cond <- sapply(combis_fast$normal, function(x) {
    # CASE 2: all normals ----
    idx <- idn[x$combis]   # active cols in y in na pattern
    id <- x$ind_i          # row indices for the na pattern
    sx <- sigman[x$combis] # sigmas for these cols
    Rx <- R[idx, idx]      # corresponding correlation matrix
    chRu <- chol(Rx)       # Cholesky upper
    D <- diag(x = sx, nrow = length(sx))
    chS <- tcrossprod(D, chRu)
    eta_n_x <- eta_n[id, x$combis, drop = FALSE]
    y_x     <- y[id, , drop = FALSE]
    Xbeta_x <- Xbeta[id, , drop = FALSE]
    C <- ltMatrices(chS[lower.tri(chS, diag = TRUE)], diag = TRUE)
    log_pl_n <- sum(ldmvnorm(t(y_x[, idx]), t(eta_n_x),
                             logLik = TRUE, chol = C))

    # CASE 3: each ordinal conditional on all normals ----
    eps <- y_x[, idx] - eta_n_x
    eps <- sweep(eps, 2, sx, "/")
    Rxinv <- chol2inv(chRu)
    log_pl_cond <- sapply(ido, function(j) {
      rho <- R[j, c(idx)]
      eta_cond <- Xbeta_x[,j] + tcrossprod(as.matrix(eps), rho %*% Rxinv)
      sd_cond  <- drop(sqrt(1 - rho %*% Rxinv %*% rho))
      eta_u_cond <- (c(thetas[[j]],  1e06)[y_x[, j]] - eta_cond)/sd_cond
      eta_l_cond <- (c(-1e06, thetas[[j]])[y_x[, j]] - eta_cond)/sd_cond
      p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)
      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps
      sum(log(p_cond), na.rm = TRUE)
    })
    log_pl_n + sum(log_pl_cond)
  })

  ## Final neg log lik
  - sum(log_pl_n_cond) - sum(log_pl_o) - sum(odens)
}

# pars <- fit2$parOpt
# dnum<-numDeriv::grad(function(par)
#   neg_log_lik_joint_type_2(par, response_types, y,
#                            X, ntheta, p, ndimo, ndimn, ndim, idn,ido,
#                            ind_univ,combis_fast), x = fit2$parOpt)
# dnum
