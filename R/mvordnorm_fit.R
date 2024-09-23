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

  ## Optimize negative log likelihood ----
  obj$res <- optimx(start_values, function(par)
    neg_log_lik_joint(par, response_types,
                      y, X, ntheta, p, ndimo, ndimn, ndim,
                      idn, ido, ind_univ, combis_fast),
    method = control$solver)

  obj$res <- optimx(start_values, function(par)
    neg_log_lik_joint(par, response_types,
                      y, X, ntheta, p, ndimo, ndimn, ndim,
                      idn, ido, ind_univ, combis_fast),
    method = control$solver)


  obj$objective <- obj$res[["value"]]
  ## TODO
  obj$parOpt <- unlist(obj$res[seq_along(start_values)])
  obj$combis_fast <- combis_fast

  ## Finalize ----
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
    obj$H.inv <-  H.inv
    obj$V <- V
    obj$claic <- 2 * obj$res[["value"]] + 2 * sum(diag(V %*% H.inv))
    obj$clbic <- 2 * obj$res[["value"]] + log(n) * sum(diag(V %*% H.inv))
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

  # These are the observations which have only one response (aka univariate)
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

grad_neg_log_lik_joint <- function(pars, response_types, y, X,
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
  Jpsi <- as.matrix(Matrix::bdiag(c(jtthetas, rep(list(1), ndimn + length(beta)))))
  Jcor <- jac_dr_dtr(tparerror, ndim)


  d_ana_univariate <- if (NROW(ind_univ) == 0) NULL else {
    ## univariate observations:
    d_ana_univariate <- colSums(Reduce("+", lapply(unique(ind_univ[, 2]), function(j) {
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

        gradmat[subj,  pos_theta_jo]  <- gradmat[subj,  pos_theta_jo] %*% Jpsi[pos_theta_jo, pos_theta_jo]

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
    })))
  }

  # iterate over pairs ----
  d_ana <- colSums(Reduce("+", lapply(1:length(combis_fast), function(i) {
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
                      nrow = n)

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
      # smat <- diag(2)
      # smat[1,2] <- smat[2, 1] <- rkl
      # smat <- tcrossprod(sigman[c(kidn, lidn)]) * smat
      # smatinv <- chol2inv(chol(smat))
      ###############################
      ## dbeta0 and dbeta for pair kl
      ###############################
      sigma_c2 <- sigma_c^2
      epsk <- (y[indkl, k] - eta_n[indkl, kidn])
      epsl <- (y[indkl, l] - eta_n[indkl, lidn])
      A <- epsk^2/sigman[kidn]^2 -
        2*rkl*epsk*epsl/(sigman[kidn]*sigman[lidn]) +
        epsl^2/sigman[lidn]^2

      dbetak <-  - 1/sigma_c2*  Xn[indkl, ] *
        (epsk/(sigman[kidn]^2) - rkl * epsl/(sigman[kidn] * sigman[lidn]))

      dbetal <-  - 1/sigma_c2 * Xn[indkl, ] *
        (epsl/(sigman[lidn]^2) - rkl * epsk/(sigman[lidn] * sigman[kidn]))


      pos_beta_k <- c(sum(ntheta) + kidn, sum(ntheta) + ndimn + k + ndim * (1:p - 1))
      pos_beta_l <- c(sum(ntheta) + lidn, sum(ntheta) + ndimn + l + ndim * (1:p - 1))
      gradmat[indkl, pos_beta_k] <- dbetak
      gradmat[indkl, pos_beta_l] <- dbetal

      ##################
      ## dcorr and dsigma
      ##################
      # see https://stats.stackexchange.com/questions/27436/how-to-take-derivative-of-multivariate-normal-density
      # dLdSigma <- do.call("rbind", lapply(indkl, function(j){
      #   SinvxxtSinv <- smatinv %*% tcrossprod(y[j, c(k,l)] - eta_n[j, c(kidn, lidn)]) %*% smatinv
      #   res <- (1/2 * (2 * smatinv - diag(diag(smatinv)) - 2 * SinvxxtSinv + diag(diag(SinvxxtSinv))))
      #   res[lower.tri(res, diag = TRUE)] # it is symmetric, we only need the lower triangle
      # }))
      # dSigmadR <- matrix(c(2 * sigman[kidn], 0, 0,
      #     rkl * sigman[lidn],rkl * sigman[kidn], sigman[kidn] * sigman[lidn],
      #     0, 2 * sigman[lidn], 0), byrow = TRUE, ncol = 3)
      # dLdR <- tcrossprod(dLdSigma,  dSigmadR)
      dr <- -rkl/sigma_c2 + rkl/(sigma_c2)^2 * A - 1/sigma_c2 * epsl*epsk/((sigman[kidn]*sigman[lidn]))
      pos_r_kl   <- sum(ntheta) + 2 * ndimn + ndim * p + x$rpos
      gradmat[indkl,  pos_r_kl]  <- dr
      gradmat[indkl,  pos_corrs]  <- gradmat[indkl,  pos_corrs]  %*% Jcor

      dsigmak <- 1/sigman[kidn] + (- epsk^2/(sigman[kidn]^3) + rkl * epsk * epsl/(sigman[kidn]^2 * sigman[lidn]))/sigma_c2
      dsigmal <- 1/sigman[lidn] + (- epsl^2/(sigman[lidn]^3) + rkl * epsk * epsl/(sigman[kidn] * sigman[lidn]^2))/sigma_c2
      dsigmak <- dsigmak * sigman[kidn]
      dsigmal <- dsigmal * sigman[lidn]
      pos_sdn_kl <- sum(ntheta) + ndimn + ndim * p + c(kidn, lidn)
      gradmat[indkl, pos_sdn_kl] <- cbind(dsigmak, dsigmal)
      # colSums(cbind(dbetak,dbetal,dsigmak, dsigmal,dr))
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

      grad_theta <- gradmat[indkl, pos_theta_ko] %*% Jpsi[pos_theta_ko,pos_theta_ko]
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
  })))
  d_ana + d_ana_univariate
}

# pars <- fit$parOpt
# d_num<-numDeriv::grad(function(par) neg_log_lik_joint(par, response_types, y, X,
#                                                       ntheta, p, ndimo, ndimn, ndim,
#                                                       idn, ido,
#                                                       ind_univ,
#                                                       combis_fast),x = pars, method = "Richardson")
#
# #cbind(d_ana, d_num)# TODO: grad_neg_log_lik_joint
# cbind(d_num, grad_neg_log_lik_joint(pars, response_types, y, X,
#                                     ntheta, p, ndimo, ndimn, ndim,
#                                     idn, ido,
#                                     ind_univ,
#                                     combis_fast)

)

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
