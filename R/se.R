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


# pars <- fit$parOpt
# y <- fit$y
# X <- fit$X
# response_types <- fit$response_types
rowwise_pairwise_grad_neg_log_lik_joint <- function(pars, response_types, y, X, Xn,
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
  Jcor <- jac_dr_dtr(tparerror, ndim)

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

    ## CASE 2: 2 normals ---
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
    ## CASE 3: 1 normals + 1 ordinal ---
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

grad_neg_log_lik_joint <- function(pars, response_types, y, X, Xn,
                                   ntheta, p, ndimo, ndimn, ndim,
                                   idn, ido,
                                   ind_univ,
                                   combis_fast) {
  colSums(Reduce("+", rowwise_pairwise_grad_neg_log_lik_joint(pars, response_types, y, X, Xn,
                                                              ntheta, p, ndimo, ndimn, ndim,
                                                              idn, ido,
                                                              ind_univ,
                                                              combis_fast)))
}
