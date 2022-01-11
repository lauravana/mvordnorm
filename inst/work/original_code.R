## Libraries
library(pbivnorm)
library(mvtnorm)
library(optimx)

#############################
### 1. Simulation of data ###
#############################

set.seed(12345)
n <- 1000
p <- 3
qo <- 2 # no of ordinal
qn <- 2 # no of normals
X <- matrix(rnorm(p * n), ncol = p)
colnames(X) <- paste0("X", 1:ncol(X))
beta <- c(2, 0, -2)

beta0n <- c(-1, 1)
theta1 <-  c(-1, 1)
theta2 <- c(-2, 2)


# S <- rWishart(1, df =6,Sigma = diag(qo + qn))[,,1]

S <- matrix(c(1.0,  0.7,  0.8,  0.7,
  0.7,  1.0,  0.9,  0.8,
  0.8,  0.9,  1.0,  0.9,
  0.7,  0.8,  0.9,  1.0), nrow = 4)


err <- mvtnorm::rmvnorm(n, mean = rep(0, qo + qn), sigma = cov2cor(S))
y1tilde <- X %*% beta + err[,1]
y2tilde <- X %*% beta + err[,2]
y1 <- as.numeric(cut(y1tilde, c(-Inf, theta1, Inf)))
y2 <- as.numeric(cut(y2tilde, c(-Inf, theta2, Inf)))

sigma1 <- 1
sigma2 <- 2
z1 <- beta0n[1] + X %*% beta + sigma1 * err[,3]
z2 <- beta0n[2] + X %*% beta + sigma2 * err[,4]

df <- cbind.data.frame(y1, y2, z1, z2, X)
head(df)


################


##
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
  p1 <- pbivnorm(U[[1]], U[[2]], r)
  p2 <- pbivnorm(L[[1]], U[[2]], r)
  p3 <- pbivnorm(U[[1]], L[[2]], r)
  p4 <- pbivnorm(L[[1]], L[[2]], r)
  ## replace NaN
  p1[is.nan(p1)] <- 0
  p2[is.nan(p2)] <- 0
  p3[is.nan(p3)] <- 0
  p4[is.nan(p4)] <- 0
  pr <- p1 - p2 - p3 + p4
  return(pr)
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





make_start_values <- function(y, X, family_type) {
  start_theta <- lapply(which(family_type == "ordinal"),
                       function(j) coef(ordinal::clm(factor(y[, j])~1)))

  pars <- c(unlist(lapply(start_theta, function(x) c(x[1], log(diff(x))))),
            colMeans(y[, family_type=="normal"]),
            # beta0 for normals
            rep(0, sum(family_type == "normal")),
            # betas
            rep(0, ncol(X)),
            # sigmas for normales
            rep(0, ncol(X)),
            # correlation params
            rep(0,length(S[lower.tri(S)])))
  return(pars)
}



log_lik_joint <- function(pars, family_type, y, X) {
  ndimo <- sum(family_type=="ordinal")
  ndim <- ncol(y)
  ndimn <- ndim  - ndimo

  p <- ncol(X)
  ## number of thresholds
  ntheta <- apply(y[,family_type=="ordinal", drop = FALSE], 2,
                  function(x) nlevels(as.factor(x)) - 1)

  ## thresholds
  thetas <- lapply(1:sum(family_type=="ordinal"), function(j) {
    transf_thresholds_flexible(pars[cumsum(c(0,ntheta))[j] + 1:ntheta[j]])
  })
  ## regression coefs: intercepts for normals
  beta0n <- pars[sum(ntheta) + seq_len(ndimn)] # we need intercepts for the normal variables
  ## coomon regression coefs
  beta <- pars[sum(ntheta) + ndimn + 1:p]
  ## sd parameters for normals
  sigman <- exp(pars[sum(ntheta) + ndimn + p + seq_len(ndimn)])
  ## correlations error structure
  tparerror <- pars[(sum(ntheta) + ndimn + p + ndimn + seq_len(ndim * (ndim - 1)/2))]
  rvec <- transf_sigma(tparerror, ndim)

  Xbeta <- drop(X %*% beta)
  eta_u <- lapply(1:ndim, function(j) {
    if (family_type[j] == "ordinal") {
      c(thetas[[j]], 1e06)[y[,j]] - Xbeta
    } else {
      beta0n[j-ndimo] + Xbeta
    }
  })
  eta_l <- lapply(1:ndim, function(j) {
    if (family_type[j] == "ordinal") {
      drop(c(-1e06, thetas[[j]])[y[,j]] - Xbeta)
    }
  })

  log_pl_vec <- NULL
  combis <- combn(1:ncol(y), 2)
  for (s in 1:ncol(combis)) {
    k <- combis[1, s]
    l <- combis[2, s]
    r <- rvec[s]

    ## CASE 1: 2 ordinals
    if (all(family_type[c(k,l)] == "ordinal")) {
      prs <- rectbiv_norm_prob(U = eta_u[c(k, l)], L = eta_l[c(k,l)], r)
      prs[prs < .Machine$double.eps] <- .Machine$double.eps
      log_pl_vec[s] <- sum(log(prs))
    }
    ## CASE 2: 2 normals
    if (all(family_type[c(k,l)] == "normal")) {
      ykstd <- (y[,k] - eta_u[[k]])
      ylstd  <- (y[,l] -  eta_u[[l]])
      smat <- diag(2)
      smat[1,2] <- smat[2, 1] <- r
      smat <- tcrossprod(sigman[c(k, l) - ndimo]) * smat
      log_pl_vec[s] <-  sum(mvtnorm::dmvnorm(x = cbind(ykstd, ylstd),
                                             mean =rep(0,2),
                                             sigma = smat,
                                             log = TRUE))
    }
    ## CASE 3: ordinal + normal
    if (family_type[k] != family_type[l]) {
      ko <- c(k,l)[which(family_type[c(k,l)]=="ordinal")]
      kn <- c(k,l)[which(family_type[c(k,l)]!="ordinal")]

      ld_marg <- dnorm(y[,kn], mean = eta_u[[kn]],
                       sd = sigman[kn-ndimo],
                       log = TRUE)

      sd_yn_cond <- sqrt(1 - r^2)
      y_tilde_cond <- Xbeta + r * (y[,kn] - eta_u[[kn]])/(sigman[kn-ndimo])
      y_cond_eta_u <- (c(thetas[[ko]],  1e06)[y[,ko]] - y_tilde_cond)/sd_yn_cond
      y_cond_eta_l <- (c(-1e06, thetas[[ko]])[y[,ko]] - y_tilde_cond)/sd_yn_cond
      p_marg <- pnorm(y_cond_eta_u) - pnorm(y_cond_eta_l)
      p_marg[p_marg < .Machine$double.eps] <- .Machine$double.eps
      log_pl_vec[s] <-  sum(log(p_marg) + ld_marg)
    }
  }
  - sum(log_pl_vec)
}


###### FUNCTION CALL

X <- df[, grepl("^X", colnames(df))]
y <- df[, c("y1",
            "y2",
            "z1", "z2")]
family_type <- c("ordinal",
                 "ordinal",
                 "normal",
                 "normal")
pars0 <- make_start_values(y, X, family_type)
pars <- pars0

optRes <- optimx(pars0,  function(x)
  log_lik_joint(x, family_type, y = y, X = as.matrix(X)),
  method = "CG")


pars_true <- c(theta1[1], log(diff(theta1)),
  theta2[1], log(diff(theta2)),
               beta0n,
               beta,
               log(sigma1), log(sigma2),
               backtransf_sigmas(cov2cor(S[1:4, 1:4])))

cbind(estimate = unlist(optRes[seq_along(pars_true)]), true = pars_true)
pars <- unlist(optRes[seq_along(pars_true)])


#####################
## Standard errors ##
#####################

tparHess <- numDeriv::hessian(function(x)
  log_lik_joint(x, family_type, y = y, X = as.matrix(X)), x = pars)

J <- as.matrix(Matrix::bdiag(list(
  ## dttheta/dtheta
  diag(obj$dim$Pstar), ## d tbeta0/d beta0
  diag(obj$dim$Pstar), ## d tbeta0/d beta0
  dttau2.tau(tau), ## d ttau2/d tau
  dtcor.cor(obj$cor_structure, gamma),
        dtomega.omega(omega),# diag(dtomega.omega(omega))),## d tomega/d omega
                                  diag(obj$dims$nlambda) # d tlambda/d lambda
)))
H <- crossprod(J, tparHess) %*% J
obj$vcov <- tryCatch(chol2inv(chol(H)),
                     error=function(e) {
                       warning("\nCondition number close to zero! Hessian is approximated by nearest positive semidefinite matrix.\n")
                       chol2inv(chol(Matrix::nearPD(H)$mat))
                     }
)
