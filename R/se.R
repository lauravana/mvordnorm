d_rect <- function(U1, U2, L1, L2, r,
                   dUmat, dLmat, d_biv_fun) {
  UU <- d_biv_fun(U1, U2, r)
  UL <- d_biv_fun(U1, L2, r)
  LU <- d_biv_fun(L1, U2, r)
  LL <- d_biv_fun(L1, L2, r)
  ((UU - UL) * dUmat  - (LU - LL) * dLmat)
}

d_corr_rect <- function(Uk, Ul, Lk, Ll, r, fun) {
  fun(Uk, Ul, r) - fun(Uk, Ll, r) - fun(Lk, Ul, r) + fun(Lk, Ll, r)
}

dF2dx <- function(x1, x2, r) dnorm(x1) * pnorm((x2 - r * x1)/sqrt(1 - r^2))

dF2dr <- function(x, y, r){
  1/(2 * pi * sqrt(1 - r^2)) *
    exp(-(x^2 - 2 * r * x * y + y^2)/(2 * (1 - r^2)))
}

dphi2dx <- function(x, y, rho, sigma_x, sigma_y){
  # derivative of centered bivariate normal pdf wrt to x
  1/(1 - rho^2) * (1/sigma_x^2 * x - rho/(sigma_x * sigma_y) * y)
}
dphi2dsigmax <- function(x, y, rho, sigma_x, sigma_y){
  # derivative of centered bivariate normal pdf wrt to sd parameter of x
  1/sigma_x + 1/(1 - rho^2) * (- 1/sigma_x^3 * x^2 + rho/(sigma_x^2 * sigma_y) * x * y)
}

dphi2drho <- function(x, y, rho, sigma_x, sigma_y){
  # derivative of centered bivariate normal pdf wrt to correlation parameter
  a <- (x/sigma_x)^2 - 2 * rho * (x/sigma_x) *  (y/sigma_y) +  (y/sigma_y)^2
  - rho/(1 - rho^2) + rho/(1 - rho^2)^2 * a - 1/(1 - rho^2) * (x/sigma_x) *  (y/sigma_y)
}


# pars <- fit2$parOpt
# y <- fit2$y
# X <- fit2$X
# response_types <- fit$response_types
rowwise_pairwise_grad_neg_log_lik_joint_type_2 <- function(pars, response_types, y, X, Xn,
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
  pos_corrs <- (sum(ntheta) + 2 * ndimn + ndim * p + seq_len(ndim * (ndim - 1)/2))
  tparerror <- pars[pos_corrs]
  rvec <- transf_sigma(tparerror, ndim)
  smat <- diag(nrow = ndim)
  smat[lower.tri(smat)] <- rvec
  R <- smat + t(smat) - diag(ndim)

  beta_mat <- beta
  dim(beta_mat) <- c(ndim, p)
  #print(beta_mat)
  Xbeta <- tcrossprod(X, beta_mat)
  XU <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BU <- B2[,-(ntheta[j] + 1), drop = FALSE]
    cbind(BU, - X)
  })
  XL <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BL <- B2[,-1, drop = FALSE]
    cbind(BL, - X)
  })

  ## X matrix for normal variables to include intercept
  Xn <- cbind(1, X)

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

  jtthetas <-
    jac_dtheta_dttheta_flexible(pars[seq_len(sum(ntheta))],
                                ndimo = ndimo, ntheta = ntheta)
  Jpsi <- as.matrix(Matrix::bdiag(jtthetas))
  Jcor <- (jac_dr_dtr(tparerror, ndim))

  # Univariate ordinal ----
  d_ana_univariate_o <- if (NROW(ind_univ) == 0) NULL else {
    ## univariate observations:
    lapply(unique(ind_univ[, 2]), function(j) {
      gradmat <- matrix(0, nrow = nrow(y), ncol = length(pars))
      idj <- ind_univ[ind_univ[, 2] == j, , drop = F]
      subj <- idj[, 1] # subject ids for j-th response
      #if (response_types[j] == "ordinal") {
      jido <- which(ido == j)
      U <- eta_u_o[subj, jido]
      L <- eta_l_o[subj, jido]
      xUj <- XU[[jido]][subj, ]
      xLj <- XL[[jido]][subj, ]
      pr <- pnorm(U) - pnorm(L)
      pr[pr < .Machine$double.eps] <- .Machine$double.eps

      pos_theta_jo <- cumsum(c(0, ntheta))[jido] + seq_len(ntheta[jido])
      pos_beta_jo <- sum(ntheta) + ndimn + j + ndim * (1:p - 1)
      gradmat[subj, c(pos_theta_jo, pos_beta_jo)] <-
        1/pr * (dnorm(L) * xLj -  dnorm(U) * xUj)

      gradmat[subj,  pos_theta_jo]  <- gradmat[subj,  pos_theta_jo] %*%
        Jpsi[pos_theta_jo, pos_theta_jo, drop=FALSE]
      gradmat
    })
  }
  # Bivariate ordinal ----
  d_ana_o <- lapply(combis_fast$ordinal, function(x) {
    comb <- x$combis
    ## for each pair make an indicator for each subject where the pair applies
    indkl <- x$ind_i
    ## correlation
    rkl <- rvec[x$rpos]
    sigma_c <- sqrt(1 - rkl^2)
    k <- comb[1]
    l <- comb[2]

    gradmat <- matrix(0,
                      ncol = length(pars),
                      nrow =  nrow(y))

    kido <- which(ido == k)
    lido <- which(ido == l)

    Uk <- eta_u_o[indkl, kido] ## upper predictor k-th response
    Ul <- eta_u_o[indkl, lido] ## upper predictor l-th response
    Lk <- eta_l_o[indkl, kido] ## lower predictor k-th response
    Ll <- eta_l_o[indkl, lido] ## lower predictor l-th response
    XUk <- XU[[kido]][indkl, ]
    XUl <- XU[[lido]][indkl, ]
    XLk <- XL[[kido]][indkl, ]
    XLl <- XL[[lido]][indkl, ]
    ## pr_{kl}
    pr <- rectbiv_norm_prob(U = cbind(Uk, Ul),
                            L = cbind(Lk, Ll), rkl)

    pr[pr < .Machine$double.eps] <- .Machine$double.eps
    ###############################
    ## dtheta and dbeta for pair kl
    ###############################
    pos_theta_k <- cumsum(c(0, ntheta))[kido] + seq_len(ntheta[kido])
    pos_theta_l <- cumsum(c(0, ntheta))[lido] + seq_len(ntheta[lido])
    pos_beta_k <- sum(ntheta) + ndimn + k + ndim * (seq_len(p) - 1)
    pos_beta_l <- sum(ntheta) + ndimn + l + ndim * (seq_len(p) - 1)

    gradmat[indkl, c(pos_theta_k, pos_beta_k)] <- - 1/pr * d_rect(
      U1 = Uk, L1 = Lk,
      U2 = Ul, L2= Ll,
      r = rkl,
      dUmat = XUk,
      dLmat = XLk,
      d_biv_fun = dF2dx)
    gradmat[indkl, c(pos_theta_l, pos_beta_l)] <- - 1/pr * d_rect(
      U1 = Ul, L1 = Ll, U2 = Uk, L2 = Lk, r = rkl,
      dUmat = XUl,
      dLmat = XLl,
      d_biv_fun = dF2dx)

    id_theta_kl <- c(pos_theta_k, pos_theta_l)

    grad_theta <- gradmat[indkl,  id_theta_kl] %*% Jpsi[id_theta_kl, id_theta_kl]
    gradmat[indkl,  id_theta_kl] <- grad_theta

    ##################
    ## dcorr
    ##################
    pos_corr_kl <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
    gradmat[indkl, pos_corr_kl] <- - 1/pr * d_corr_rect(Uk, Ul, Lk, Ll, r = rkl, dF2dr)
    gradmat[indkl,  pos_corrs]  <- gradmat[indkl,  pos_corrs]  %*% Jcor

    gradmat
  })

  gradmat_n <- matrix(0,
                      ncol = length(pars),
                      nrow =  nrow(y))
  gradmat_cond <- rep(list(matrix(0, ncol = length(pars),
                                  nrow =  nrow(y))), ndimo)
  for (i in seq_along(combis_fast$normal)) {
    x <- combis_fast$normal[[i]]
    # CASE 2 DONE: normals ----
    idx <- idn[x$combis]   # active cols in y in na pattern
    nidx <- length(idx)
    id <- x$ind_i          # row indices for the na pattern
    sx <- sigman[x$combis] # sigmas for these cols
    Rx <- R[idx, idx]      # corresponding correlation matrix
    chRu <- chol(Rx)       # Cholesky upper
    D <- diag(x = sx, nrow = length(sx))
    Dinv <- diag(x = 1/sx, nrow = length(sx))
    chS <- tcrossprod(D, chRu)
    Rxinv <- chol2inv(chRu)
    Sinv <- chol2inv(t(chS))
    eta_n_x <- eta_n[id, x$combis, drop = FALSE]
    y_x     <- y[id, , drop = FALSE]
    z <- as.matrix(y_x[, idx] - eta_n_x)

    ## dbeta ----
    s_mu <- z %*% Sinv
    s_beta <- - s_mu[, rep(seq_len(ncol(s_mu)), each = ncol(Xn))] *
      Xn[id, rep(seq_len(ncol(Xn)), times = ncol(s_mu))]

    pos_beta_n <- unlist(lapply(seq_along(idx), function(j)
      c(sum(ntheta) + x$combis[j], sum(ntheta) + ndimn + idx[j] +
          ndim * (seq_len(p) - 1))))

    gradmat_n[id, pos_beta_n] <- s_beta

    ## dR ----
    eps <- as.matrix(y_x[, idx] - eta_n_x)
    eps <- sweep(eps, 2, sx, "/")
    epsRinv <- (eps %*% Rxinv)
    ncorn <- nidx * (nidx - 1)/2
    if (ncorn > 0) {
      s_R <- matrix(nrow = length(id), ncol = nidx * (nidx - 1)/2)
      h <- 0
      for (j in 1:(nidx - 1)) {
        for (k in (j + 1):nidx) {
          h <- h + 1
          s_R[, h] <- (Rxinv[k, j] - epsRinv[, k] * epsRinv[, j])
        }
      }
      id_pos_n_r <- which(apply(combn(ndim, 2), 2, function(x) all(x %in% idx)))
      pos_r   <- sum(ntheta) + 2 * ndimn + ndim * p + id_pos_n_r
      gradmat_n[id,  pos_r]  <- s_R
    }
    # correct back to orig param
    gradmat_n[id,  pos_corrs]  <- gradmat_n[id,  pos_corrs]  %*% Jcor

    ## dsigma ----
    epsRinveps <- eps *  (eps %*% Rxinv)
    epsRinveps <- sweep(epsRinveps, 2, - sx,  "/")
    s_sigma <- sweep(epsRinveps, 2, 1/sx,  "+")
    # correct back to orig param
    s_sigma <- sweep(s_sigma, 2, sx, "*")
    pos_s <- sum(ntheta) + ndimn + ndim * p + seq_len(nidx)
    gradmat_n[id, pos_s] <-  s_sigma

    # CASE 3: each ordinal conditional on all normals ----
    for (jo in ido) {
      rho <- R[jo, c(idx)]
      rhoRinv <- (rho %*% Rxinv)
      epsrhoRinv <- tcrossprod(eps, rhoRinv)
      eta_c <- Xbeta[id, jo, drop = FALSE] + epsrhoRinv
      sd_c  <- drop(sqrt(1 - rho %*% Rxinv %*% rho))
      eta_u_c <- (c(thetas[[jo]],  1e06)[y_x[, jo]] - eta_c)/sd_c
      eta_l_c <- (c(-1e06, thetas[[jo]])[y_x[, jo]] - eta_c)/sd_c
      p_cond <- pnorm(eta_u_c) - pnorm(eta_l_c)
      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps

      XUko <- XU[[jo]][id, ]
      XLko <- XL[[jo]][id, ]

      pos_theta_j <- cumsum(c(0, ntheta))[jo] + seq_len(ntheta[jo])
      pos_beta_j  <- sum(ntheta) + ndimn + jo + ndim * (seq_len(p) - 1)

      ## d beta j and d theta j ----
      gradmat_cond[[jo]][id, c(pos_theta_j, pos_beta_j)] <-
        c(1/p_cond) * (c(dnorm(eta_l_c)/sd_c) * XLko -  c(dnorm(eta_u_c)/sd_c) *
                         XUko)

      grad_theta <- gradmat_cond[[jo]][id, pos_theta_j] %*%
        Jpsi[pos_theta_j, pos_theta_j, drop=FALSE]

      gradmat_cond[[jo]][id, pos_theta_j] <- grad_theta

      ## d betas of the normals ----
      s_beta <- - s_mu[, rep(seq_len(ncol(s_mu)), each = ncol(Xn))] *
        Xn[id, rep(seq_len(ncol(Xn)), times = ncol(s_mu))]
      id_rep_response <- rep(seq_len(ncol(rhoRinv)), each = ncol(Xn))
      id_rep_p        <- rep(seq_len(ncol(Xn)), times = ncol(rhoRinv))
      vecrhoRinv <- ((rho %*% Rxinv)/sx)[id_rep_response]
      s_beta_x_n <- sweep(Xn[id, id_rep_p], 2, vecrhoRinv, "*")

      s_beta_n <-  -  c(1/p_cond) * c(dnorm(eta_u_c) - dnorm(eta_l_c)) / sd_c  *
        s_beta_x_n
      gradmat_cond[[jo]][id, pos_beta_n] <- s_beta_n

      ## drho ----
      id_pos_j_rho <- which(apply(combn(ndim, 2), 2, function(x) {
        x[1] == jo & x[2]  %in% idx
      }))
      pos_rho   <- sum(ntheta) + 2 * ndimn + ndim * p + id_pos_j_rho
      deta_u_c_drho <- (eta_u_c %*% rhoRinv / sd_c^2 -  epsRinv / sd_c)
      deta_l_c_drho <- (eta_l_c %*% rhoRinv / sd_c^2  - epsRinv / sd_c)
      s_rho <- - c(1/p_cond) *
        (drop(dnorm(eta_u_c)) * deta_u_c_drho  -
           drop(dnorm(eta_l_c)) * deta_l_c_drho)

      gradmat_cond[[jo]][id, pos_rho] <- s_rho

      ## dR ----
      if (ncorn > 0) {
        s_R <- sapply(seq_len(nidx * (nidx - 1)/2), function(k) {
          id_rhos_n <- combn(nidx, 2)
          R_k_deriv <- matrix(0, nrow = nrow(Rx), ncol = ncol(Rx))
          R_k_deriv[id_rhos_n[2, k], id_rhos_n[1, k]] <- 1
          R_k_deriv[id_rhos_n[1, k], id_rhos_n[2, k]] <- 1
          dRinvdr <- Rxinv %*% R_k_deriv %*% Rxinv
          parta <- tcrossprod(eps, rho %*% dRinvdr)/sd_c
          partb <- drop(rho %*% dRinvdr %*% rho) / (sd_c^2)

          deta_u_c_dR <- - parta + 0.5 * eta_u_c * partb
          deta_l_c_dR <- - parta + 0.5 * eta_l_c * partb

          c(1/p_cond) * (dnorm(eta_u_c) * deta_u_c_dR - dnorm(eta_l_c) * deta_l_c_dR)
        })
        gradmat_cond[[jo]][id, pos_r] <- s_R
      }

      gradmat_cond[[jo]][id, pos_corrs]  <-  gradmat_cond[[jo]][id, pos_corrs]  %*% Jcor

      ## d sigma ----
      z <- as.matrix(y_x[, idx, drop = FALSE] - eta_n_x)
      deta_c_dsigma <- z %*% diag(drop(rho  %*%  Rxinv %*% Dinv^2),
                                  nrow = nrow(Rxinv))/sd_c
      s_sigma <- - drop(1/p_cond) *
        drop(dnorm(eta_u_c) - dnorm(eta_l_c)) * deta_c_dsigma
      ## correct for transformation exp(lsigma) = sigma
      s_sigma <- sweep(s_sigma, 2, sx, "*")
      gradmat_cond[[jo]][id, pos_s] <- s_sigma

      ## Finalize ----
      gradmat_cond[[jo]][is.na(gradmat_cond[[jo]])] <- 0
    }
  }
  # cbind(colSums(Reduce("+", res)),dnum)
  res <- c(d_ana_univariate_o, d_ana_o,
           list(gradmat_n), gradmat_cond)

}


# pars<-fit2$parOpt
rowwise_pairwise_grad_neg_log_lik_joint_type_1 <- function(pars, response_types, y, X, Xn,
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
  pos_corrs <- (sum(ntheta) + 2 * ndimn + ndim * p + seq_len(ndim * (ndim - 1)/2))
  tparerror <- pars[pos_corrs]
  rvec <- transf_sigma(tparerror, ndim)

  beta_mat <- beta
  dim(beta_mat) <- c(ndim, p)
  #print(beta_mat)
  Xbeta <- tcrossprod(X, beta_mat)
  XU <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BU <- B2[,-(ntheta[j] + 1), drop = FALSE]
    cbind(BU, - X)
  })
  XL <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BL <- B2[,-1, drop = FALSE]
    cbind(BL, - X)
  })

  ## X matrix for normal variables to include intercept
  Xn <- cbind(1, X)

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

  jtthetas <-
    jac_dtheta_dttheta_flexible(pars[seq_len(sum(ntheta))],
                                ndimo = ndimo, ntheta = ntheta)
  Jpsi <- as.matrix(Matrix::bdiag(jtthetas))
  Jcor <- (jac_dr_dtr(tparerror, ndim))

  d_ana_univariate <- if (NROW(ind_univ) == 0) NULL else {
    ## univariate observations:
    lapply(unique(ind_univ[, 2]), function(j) {
      gradmat <- matrix(0, nrow = nrow(y), ncol = length(pars))
      idj <- ind_univ[ind_univ[, 2] == j, , drop = F]
      subj <- idj[, 1] # subject ids for j-th response
      if (response_types[j] == "ordinal") {
        jido <- which(ido == j)
        U <- eta_u_o[subj, jido]
        L <- eta_l_o[subj, jido]
        xUj <- XU[[jido]][subj, ]
        xLj <- XL[[jido]][subj, ]
        pr <- pnorm(U) - pnorm(L)
        pr[pr < .Machine$double.eps] <- .Machine$double.eps

        pos_theta_jo <- cumsum(c(0, ntheta))[jido] + seq_len(ntheta[jido])
        pos_beta_jo <- sum(ntheta) + ndimn + j + ndim * (1:p - 1)
        gradmat[subj, c(pos_theta_jo, pos_beta_jo)] <-
          1/pr * (dnorm(L) * xLj -  dnorm(U) * xUj)

        gradmat[subj,  pos_theta_jo]  <- gradmat[subj,  pos_theta_jo] %*% Jpsi[pos_theta_jo, pos_theta_jo, drop=FALSE]

      }

      if (response_types[j] != "ordinal") {
        jidn <- which(idn == j)
        dsigmaj <- 1/sigman[jidn] - sigman[jidn]^(-3) *
          (y[subj,j] - eta_n[subj,jidn])^2
        dbetaj <- - (y[subj,j] - eta_n[subj,jidn])/(sigman[jidn]^2) * Xn[subj, ]
        pos_sigma_j   <- sum(ntheta) + ndimn + ndim * p + jidn
        pos_beta_jn <- c(sum(ntheta) + jidn,
                         sum(ntheta) + ndimn + j + ndim * (1:p - 1))

        gradmat[subj, pos_sigma_j] <- dsigmaj * sigman[jidn]
        gradmat[subj, pos_beta_jn] <- dbetaj
      }
      gradmat
    })
  }

  # iterate over pairs ----
  d_ana <- lapply(seq_along(combis_fast), function(i) {
    x <- combis_fast[[i]]
    comb <- x$combis
    ## for each pair make an indicator for each subject where the pair applies
    indkl <- x$ind_i
    ## correlation
    rkl <- rvec[x$rpos]
    sigma_c <- sqrt(1 - rkl^2)
    k <- comb[1]
    l <- comb[2]

    gradmat <- matrix(0,
                      ncol = length(pars),
                      nrow =  nrow(y))

    ## CASE 1: 2 ordinals ---
    if (all(response_types[c(k,l)] == "ordinal")) {
      kido <- which(ido == k)
      lido <- which(ido == l)

      Uk <- eta_u_o[indkl, kido] ## upper predictor k-th response
      Ul <- eta_u_o[indkl, lido] ## upper predictor l-th response
      Lk <- eta_l_o[indkl, kido] ## lower predictor k-th response
      Ll <- eta_l_o[indkl, lido] ## lower predictor l-th response
      XUk <- XU[[kido]][indkl, ]
      XUl <- XU[[lido]][indkl, ]
      XLk <- XL[[kido]][indkl, ]
      XLl <- XL[[lido]][indkl, ]
      ## pr_{kl}
      pr <- rectbiv_norm_prob(U = cbind(Uk, Ul),
                              L = cbind(Lk, Ll), rkl)

      pr[pr < .Machine$double.eps] <- .Machine$double.eps
      ###############################
      ## dtheta and dbeta for pair kl
      ###############################
      pos_theta_k <- cumsum(c(0, ntheta))[kido] + seq_len(ntheta[kido])
      pos_theta_l <- cumsum(c(0, ntheta))[lido] + seq_len(ntheta[lido])
      pos_beta_k <- sum(ntheta) + ndimn + k + ndim * (seq_len(p) - 1)
      pos_beta_l <- sum(ntheta) + ndimn + l + ndim * (seq_len(p) - 1)

      gradmat[indkl, c(pos_theta_k, pos_beta_k)] <- - 1/pr * d_rect(
        U1 = Uk, L1 = Lk,
        U2 = Ul, L2= Ll,
        r = rkl,
        dUmat = XUk,
        dLmat = XLk,
        d_biv_fun = dF2dx)
      gradmat[indkl, c(pos_theta_l, pos_beta_l)] <- - 1/pr * d_rect(
        U1 = Ul, L1 = Ll, U2 = Uk, L2 = Lk, r = rkl,
        dUmat = XUl,
        dLmat = XLl,
        d_biv_fun = dF2dx)

      id_theta_kl <- c(pos_theta_k, pos_theta_l)

      grad_theta <- gradmat[indkl,  id_theta_kl] %*% Jpsi[id_theta_kl, id_theta_kl]
      gradmat[indkl,  id_theta_kl] <- grad_theta

      ##################
      ## dcorr
      ##################
      pos_corr_kl <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl, pos_corr_kl] <- - 1/pr * d_corr_rect(Uk, Ul, Lk, Ll, r = rkl, dF2dr)
      gradmat[indkl,  pos_corrs]  <- gradmat[indkl,  pos_corrs]  %*% Jcor
    }

    ## CASE 2: 2 normals ----
    if (all(response_types[c(k,l)] != "ordinal")) {
      kidn <- which(idn == k)
      lidn <- which(idn == l)
      ###############################
      ## dbeta0 and dbeta for pair kl
      ###############################
      epsk <- (y[indkl, k] - eta_n[indkl, kidn])
      epsl <- (y[indkl, l] - eta_n[indkl, lidn])

      dbetak <-  dphi2dx(epsk, epsl, rkl, sigman[kidn], sigman[lidn]) * (-Xn[indkl, ])
      dbetal <-  dphi2dx(epsl, epsk, rkl, sigman[lidn], sigman[kidn]) * (-Xn[indkl, ])

      pos_beta_k <- c(sum(ntheta) + kidn, sum(ntheta) + ndimn + k + ndim * (1:p - 1))
      pos_beta_l <- c(sum(ntheta) + lidn, sum(ntheta) + ndimn + l + ndim * (1:p - 1))
      gradmat[indkl, pos_beta_k] <- dbetak
      gradmat[indkl, pos_beta_l] <- dbetal

      ##################
      ## dcorr and dsigma
      ##################
      # see https://stats.stackexchange.com/questions/27436/how-to-take-derivative-of-multivariate-normal-density
      dr <- dphi2drho(epsk, epsl, rkl, sigman[kidn], sigman[lidn])
      pos_r_kl   <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl,  pos_r_kl]  <- dr
      # correct back to orig param
      gradmat[indkl,  pos_corrs]  <- gradmat[indkl,  pos_corrs]  %*% Jcor

      dsigmak <- dphi2dsigmax(epsk, epsl, rkl, sigman[kidn], sigman[lidn])
      dsigmal <- dphi2dsigmax(epsl, epsk, rkl, sigman[lidn], sigman[kidn])

      # correct back to orig param
      dsigmak <- dsigmak * sigman[kidn]
      dsigmal <- dsigmal * sigman[lidn]
      pos_sdn_kl <- sum(ntheta) + ndimn + ndim * p + c(kidn, lidn)
      gradmat[indkl, pos_sdn_kl] <- cbind(dsigmak, dsigmal)
    }
    ## CASE 3: 1 normals + 1 ordinal ----
    if (response_types[k] != response_types[l]) {
      ko <- c(k,l)[which(response_types[c(k,l)] == "ordinal")]
      kn <- c(k,l)[which(response_types[c(k,l)] != "ordinal")]
      knidn <- which(idn == kn)
      koido <- which(ido == ko)
      ld_marg <- dnorm(y[indkl,kn],
                       mean = eta_n[indkl, knidn,drop = FALSE],
                       sd = sigman[knidn],
                       log = TRUE)

      eps_kn <- y[indkl,kn] - eta_n[indkl, knidn, drop = FALSE]

      mu_cond <-  Xbeta[indkl, ko] + rkl * eps_kn / sigman[knidn]

      eta_u_c <- (c(thetas[[koido]],  1e06)[y[indkl,ko]] - mu_cond)/sigma_c
      eta_l_c <- (c(-1e06, thetas[[koido]])[y[indkl,ko]] - mu_cond)/sigma_c

      p_cond <- pnorm(eta_u_c) - pnorm(eta_l_c)
      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps

      #dbeta ko
      XUko <- XU[[koido]][indkl, ]
      XLko <- XL[[koido]][indkl, ]
      pos_theta_ko <- cumsum(c(0, ntheta))[koido] + seq_len(ntheta[koido])
      pos_beta_ko <- sum(ntheta) + ndimn + ko + ndim * (1:p - 1)

      dpsiko <- c(1/p_cond) * (c(dnorm(eta_l_c)/sigma_c) * XLko -
                                 c(dnorm(eta_u_c)/sigma_c) * XUko)

      id_psi_ko <- c(pos_theta_ko, pos_beta_ko)

      gradmat[indkl, id_psi_ko] <- dpsiko

      grad_theta <- gradmat[indkl, pos_theta_ko] %*% Jpsi[pos_theta_ko,pos_theta_ko,drop=FALSE]
      gradmat[indkl, pos_theta_ko] <- grad_theta


      # dbeta0 and dbeta kn
      pos_beta_kn <- c(sum(ntheta) + knidn, sum(ntheta) + ndimn + kn +
                         ndim * (seq_len(p) - 1))
      parta <- c(dnorm(eta_u_c) - dnorm(eta_l_c)) * rkl/(sigma_c * sigman[knidn])
      partb <-  - drop(eps_kn)/(sigman[knidn]^2) * Xn[indkl, ]

      dbetakn <-  - c(1/p_cond) * parta * Xn[indkl, ] + partb

      gradmat[indkl, pos_beta_kn] <- dbetakn

      ### sigmal
      # deriv of marginal normal
      parta <- 1/sigman[knidn] - sigman[knidn]^(-3) *  (drop(eps_kn))^2

      partb <- -c(1/p_cond) * (dnorm(eta_u_c) - dnorm(eta_l_c)) *
        rkl * sigman[knidn]^(-2) * drop(eps_kn)/sigma_c

      pos_sigma_n   <- sum(ntheta) + ndimn + ndim * p + knidn

      gradmat[indkl, pos_sigma_n] <- (parta + partb) * sigman[knidn]

      ### r

      pos_r_kl   <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      zn <- - (y[indkl,kn] - eta_n[indkl, knidn])/sigman[knidn]
      deta_u_c_dbeta <- zn / sigma_c + rkl * eta_u_c/sigma_c^2
      deta_l_c_dbeta <- zn / sigma_c + rkl * eta_l_c/sigma_c^2

      gradmat[indkl, pos_r_kl] <- -c(1/p_cond) *
        (dnorm(eta_u_c) * deta_u_c_dbeta  - dnorm(eta_l_c) * deta_l_c_dbeta)
      gradmat[indkl,  pos_corrs]  <- gradmat[indkl,  pos_corrs] %*% Jcor

    }
    gradmat
  })
  c(d_ana, d_ana_univariate)
}


grad_neg_log_lik_joint_type_1 <- function(pars, response_types, y, X, Xn,
                                          ntheta, p, ndimo, ndimn, ndim,
                                          idn, ido,
                                          ind_univ,
                                          combis_fast) {
  colSums(Reduce("+", rowwise_pairwise_grad_neg_log_lik_joint_type_1(pars, response_types, y, X, Xn,
                                                                     ntheta, p, ndimo, ndimn, ndim,
                                                                     idn, ido,
                                                                     ind_univ,
                                                                     combis_fast)))
}
grad_neg_log_lik_joint_type_2 <- function(pars, response_types,
                                          y, X,  Xn, ntheta, p, ndimo, ndimn, ndim,
                                          idn, ido, ind_univ, combis_fast){
  colSums(Reduce("+", rowwise_pairwise_grad_neg_log_lik_joint_type_2(pars, response_types, y, X, Xn,
                                                                     ntheta, p, ndimo, ndimn, ndim,
                                                                     idn, ido,
                                                                     ind_univ,
                                                                     combis_fast)))
}
