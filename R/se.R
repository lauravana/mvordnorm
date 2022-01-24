d_rect <- function(Uk, Ul, Lk, Ll, r,
                   dUmat, dLmat, d_biv_fun) {
  UU <- d_biv_fun(Uk, Ul, r)
  UL <- d_biv_fun(Uk, Ll, r)
  LU <- d_biv_fun(Lk, Ul, r)
  LL <- d_biv_fun(Lk, Ll, r)
  - ((UU - UL) * dUmat  - (LU - LL) * dLmat)
}

d_psi_rect_kl <- function(Uk, Ul, Lk, Ll, r,
                          dUkmat, dLkmat,
                          dUlmat, dLlmat,
                          d_biv_fun) {
  - d_rect(Uk = Uk, Ul = Ul, Lk = Lk, Ll = Ll, r = r,
           dUmat =  dUkmat, dLmat = dLkmat, d_biv_fun = d_biv_fun) +
    d_rect(Uk = Uk, Ul = Ul, Lk = Lk, Ll = Ll, r = r,
           dUmat =  dUlmat, dLmat = dLlmat, d_biv_fun = d_biv_fun)
}

d_corr_rect <- function(Uk, Ul, Lk, Ll, r, fun) {
  - fun(Uk, Ul, r) + fun(Uk, Ll, r) + fun(Lk, Ul, r) - fun(Lk, Ll, r)
}

dF2dx <- function(x, y, r) dnorm(x) * pnorm((y - r * x)/sqrt(1 - r^2))

dF2dr <- function(x, y, r){
  1/(2 * pi * sqrt(1 - r^2)) *
    exp(-(x^2 - 2 * r * x * y + y^2)/(2 * (1 - r^2)))}

# pars <- obj$parOpt
derivs_ana <- function(pars, y, X, response_types, ind_univ,
                       combis_fast){

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
  ############################################
  ## function for analytic gradient and hessian
  #############################################
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
  XU <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BU <- B2[,-(ntheta[j] + 1), drop = FALSE]
    cbind(BU, - X)
  })
  XL <- lapply(1:ndimo, function(j) {
    B2 <- (col(matrix(0, nrow(y), ntheta[j] + 1)) == c(unclass(y[, ido[j]])))
    BL <- B2[,-1, drop = FALSE]
    cbind(BL, -X)
  })
  Xn <- lapply(1:ndimn, function(j) {
    cbind(1, X)
  })

  eta_u_o <- sapply(1:ndimo, function(j)
    c(thetas[[j]], 1e06)[y[, ido[j]]] - Xbeta[, ido[j]])

  eta_l_o <- sapply(1:ndimo, function(j)
    drop(c(-1e06, thetas[[j]])[y[, ido[j]]] - Xbeta[, ido[j]]))

  eta_n <- sapply(1:ndimn, function(j) {
    beta0n[j] + Xbeta[, idn[j]]
  })

  y <- as.matrix(y)
  ######################################################
  ## First the univariate case (q_i = 1)
  ######################################################
  g_list <- NULL
  if (NROW(ind_univ) > 0){
    ## univariate observations:
    g_list <- lapply(unique(ind_univ[, 2]), function(j) {
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
      }

      if (response_types[j] != "ordinal") {
        jidn <- which(idn == j)
        dsigmaj <- -1/sigman[jidn]^2 - sigman[jidn]^(-3) *
          (y[subj,j] - eta_n[subj,jidn])^2
        dbetaj <- - (y[subj,j] - eta_n[subj,jidn])/sigman[jidn] * Xn[[jidn]][subj, ]
        pos_sigma_j   <- sum(ntheta) + ndimn + ndim * p + jidn
        pos_beta_jn <- c(sum(ntheta) + jidn,
                         sum(ntheta) + ndimn + j + ndim * (1:p - 1))

        gradmat[subj, pos_sigma_j] <- dsigmaj
        gradmat[subj, pos_beta_jn] <- dbetaj
      }
      gradmat
    })
  }
  #####################################
  ## take each possible pair (k, l)
  ######################################
  it0 <- length(g_list)
  #it <- 1
  for (i in 1:length(combis_fast)) {
    x <- combis_fast[[i]]
    comb <- x$combis
    ## for each pair make an indicator for each subject where the pair applies
    indkl <- x$ind_i
    ## correlation
    rkl <- rvec[x$rpos]
    k <- comb[1]
    l <- comb[2]
    gradmat <- matrix(0,
                   ncol = length(pars),
                   nrow = n)
    ## CASE 1: 2 ordinals
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

      ## vector h_kl will contain the gradient for all d -log p_{kl}/d pars

      ###############################
      ## dtheta and dbeta for pair kl
      ###############################
      pos_theta_k <- cumsum(c(0, ntheta))[kido] + seq_len(ntheta[kido])
      pos_theta_l <- cumsum(c(0, ntheta))[lido] + seq_len(ntheta[lido])
      pos_beta_k <- sum(ntheta) + ndimn + k + ndim * (1:p - 1)
      pos_beta_l <- sum(ntheta) + ndimn + l + ndim * (1:p - 1)

      gradmat[indkl, c(pos_theta_k, pos_beta_k)] <- 1/pr * d_rect(Uk = Uk, Lk = Lk,
                                                           Ul = Ul, Ll= Ll,
                                                           r = rkl,
                                                           dUmat = XUk,
                                                           dLmat = XLk,
                                                           d_biv_fun = dF2dx)
      gradmat[indkl, c(pos_theta_l, pos_beta_l)] <- 1/pr * d_rect(Uk = Uk, Lk = Lk,
                                                           Ul = Ul, Ll= Ll,
                                                           r = rkl,
                                                           dUmat = XUl,
                                                           dLmat = XLl,
                                                           d_biv_fun = dF2dx)

      ##################
      ## dcorr
      ##################
      pos_corr_kl <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl, pos_corr_kl] <- 1/pr * d_corr_rect(Uk, Ul, Lk, Ll, r = rkl, dF2dr)
    }
    ## CASE 2: 2 normals
    if (all(response_types[c(k,l)] != "ordinal")) {
      kidn <- which(idn == k)
      lidn <- which(idn == l)
      Xnkl <- Xn[[kidn]][indkl, ]
      smat <- diag(2)
      smat[1,2] <- smat[2, 1] <- rkl
      smat <- tcrossprod(sigman[c(kidn, lidn)]) * smat
      smatinv <- chol2inv(chol(smat))
      ###############################
      ## dbeta0 and dbeta for pair kl
      ###############################
      dbeta <- lapply(indkl, function(i){
        (-2 * tcrossprod(Xn[[kidn]][i, ], y[i, c(k,l)])  +
          2 * tcrossprod(Xn[[kidn]][i, ]) %*% t(cbind(beta0n, beta_mat[c(k,l),]))) %*% smatinv
      })
      dbetak <- t(sapply(dbeta, function(x) x[, 1]))
      dbetal <- t(sapply(dbeta, function(x) x[, 2]))

      ###############################
      ## dtheta and dbeta for pair kl
      ###############################
      pos_beta_k <- c(sum(ntheta) + kidn, sum(ntheta) + ndimn + k + ndim * (1:p - 1))
      pos_beta_l <- c(sum(ntheta) + lidn, sum(ntheta) + ndimn + l + ndim * (1:p - 1))
      gradmat[indkl, pos_beta_k] <- dbetak
      gradmat[indkl, pos_beta_l] <- dbetal

      ##################
      ## dcorr and dsigma
      ##################
      # see https://stats.stackexchange.com/questions/27436/how-to-take-derivative-of-multivariate-normal-density
      dLdSigma <- do.call("rbind", lapply(indkl, function(j){
        A <- -1/2 * (smatinv - smatinv %*% tcrossprod(y[j, c(k,l)] - eta_n[j, c(kidn, lidn)])
                       %*% smatinv)
        B <- 2 * A - diag(diag(A))
        c(B)
      }))

      f <- function(x) {
        xmat <- matrix(x, nrow=2L)
        r <- cov2cor(xmat)[1,2]
        sigma <- sqrt(diag(xmat))
        c(sigma, r)
      }
      dSigmadR <- numDeriv::jacobian(f, x=c(smat))
      dLdR <- tcrossprod(dLdSigma,  dSigmadR)

      pos_sdn_kl <- sum(ntheta) + ndimn + ndim * p + c(kidn, lidn)
      pos_r_kl   <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl, c(pos_sdn_kl, pos_r_kl)] <- dLdR
    }

    ## CASE 3: 1 normals + 1 ordinal
    if (response_types[k] != response_types[l]) {
      ko <- c(k,l)[which(response_types[c(k,l)] == "ordinal")]
      kn <- c(k,l)[which(response_types[c(k,l)] != "ordinal")]
      knidn <- which(idn == kn)
      koido <- which(ido == ko)
      ld_marg <- dnorm(y[indkl,kn],
                       mean = eta_n[indkl, knidn,drop = FALSE],
                       sd = sigman[knidn],
                       log = TRUE)

      sd_yn_cond <- sqrt(1 - rkl^2)

      eta_cond <-  Xbeta[indkl, ko] +
        rkl * (y[indkl,kn] - eta_n[indkl, knidn, drop = FALSE]) / sigman[knidn]

      eta_u_cond <- (c(thetas[[koido]],  1e06)[y[indkl,ko]] - eta_cond)/sd_yn_cond
      eta_l_cond <- (c(-1e06, thetas[[koido]])[y[indkl,ko]] - eta_cond)/sd_yn_cond

      p_cond <- pnorm(eta_u_cond) - pnorm(eta_l_cond)
      p_cond[p_cond < .Machine$double.eps] <- .Machine$double.eps

      #dbeta ko
      XUko <- XU[[koido]][indkl, ]
      XLko <- XL[[koido]][indkl, ]
      pos_theta_ko <- cumsum(c(0, ntheta))[koido] + seq_len(ntheta[koido])
      pos_beta_ko <- sum(ntheta) + ndimn + ko + ndim * (1:p - 1)

      dpsiko <- c(1/p_cond) * (c(dnorm(eta_l_cond)/sd_yn_cond) * XUko -
                               c(dnorm(eta_u_cond)/sd_yn_cond) * XLko)
      gradmat[indkl, c(pos_theta_ko, pos_beta_ko)] <- dpsiko

      # dbeta0 and dbeta kn
      pos_beta_kn <- c(sum(ntheta) + knidn, sum(ntheta) + ndimn + kn + ndim * (1:p - 1))

      dbetakn <- (c(1/p_cond) * c(dnorm(eta_u_cond)/(sd_yn_cond * sigman[knidn]) -
                     dnorm(eta_l_cond)/(sd_yn_cond * sigman[knidn]))  * Xn[[knidn]][indkl, ]) -
        (y[indkl,kn] - eta_n[indkl,knidn])/sigman[knidn] * Xn[[knidn]][indkl, ]

      gradmat[indkl, pos_beta_kn] <- dbetakn

      ### sigmal
      # deriv of marginal normal
      parta <- -1/sigman[knidn]^2 - sigman[knidn]^(-3) *  (y[indkl,kn] - eta_n[indkl,knidn])^2
      #
      partb <- c(1/p_cond) * (dnorm(eta_u_cond) - dnorm(eta_l_cond)) *
        rkl * sigman[knidn]^(-2) * (y[indkl,kn] - eta_n[indkl,knidn])/sd_yn_cond
      pos_sigma_n   <- sum(ntheta) + ndimn + ndim * p + knidn
      gradmat[indkl, pos_sigma_n] <- parta + partb
      ### r
      pos_r_kl   <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl, pos_r_kl] <- c(1/p_cond) * (dnorm(eta_u_cond) - dnorm(eta_l_cond)) *
        ((y[indkl,kn] - eta_n[indkl, knidn])/sigman[knidn] + rkl/sqrt(1-rkl^2))/(1-rkl^2)
    }
    g_list[[i + it0]] <- gradmat
  }
  ## matrix containing the gradients for each subject
  Vi <- Reduce("+", g_list)
  ## Variability matrix
  V <- crossprod(Vi)
  ## Hessian matrix
  H <- Reduce("+", lapply(g_list, crossprod))
  list(V = V, H = H)
}
