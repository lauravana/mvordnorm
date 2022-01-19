## DATASETS
#' @title Data set toy example
#' @description A simulated data set containing 8 variables (5 responses variables and 3 covariates).
#' \itemize{
#'   \item \code{z1} numeric, gaussian response variable
#'   \item \code{z2} numeric, gaussian response variable
#'   \item \code{y1} integer, ordinal response variable
#'   \item \code{y2} integer, ordinal response variable
#'   \item \code{X1} continuous covariate
#'   \item \code{X2} continuous covariate
#'   \item \code{X3} continuous covariate
#'}
#' @name data_toy
#' @docType data
#' @usage data("data_toy", package = "mvordnorm")
#' @format A data frame with 1000 rows and 7 variables
NULL

## IMPORTS
#' @importFrom stats family gaussian as.formula dnorm dcauchy model.offset pnorm printCoefmat model.matrix model.weights nlminb setNames
#' @importFrom utils write.table combn
#' @importFrom ordinal clm
#' @importFrom numDeriv hessian grad
#' @importFrom pbivnorm pbivnorm
#' @importFrom mvtnorm dmvnorm
#' @importFrom Matrix bdiag nearPD


#############################################################################################
#' @title Fitting a multivariate model of mixed normal and ordinal responses
#' @description mvordnorm is used to fit a multivariate model of mixed normal and ordinal random variables. Probit link is used for the ordinal models
#' @param formula A formula as in the Formula package?
#' @param data data
#' @param response_types a (named) vector of characters. Each element of the vector is
#'  either "gaussian" or "ordinal".
#' @param na.action eliminate NAs or keep them
#' @param weights  weights which need to be constant across multiple measurements. Negative weights are not allowed.
#' @param offset this can be used to specify an a priori known component to be included in the linear predictor during fitting.
#' @param contrasts an optional list. See the \code{contrasts.arg} of \code{\link{model.matrix.default}}.
#' @param control  list of parameters for controlling the fitting process. See \code{\link{mixoglmm.control}} for details.
#' @param ...  additional arguments.
#' @return an object of class mvordnorm
#' @examples
#' sum(1:10)
#'
#' \dontrun{
#' sum("a")
#' }
#' @export
#'
mvordnorm <- function(formula, data,
                      response_types = NULL,
                      na.action, weights,
                      offset = NULL,
                      contrasts = NULL,
                      control = mvordnorm.control(),
                      ...) {
  ## call
  call <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  oformula <- as.formula(formula)
  formula <- Formula::as.Formula(formula)
  if(length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
    simple_formula <- TRUE
  } else {
    if(length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
    simple_formula <- FALSE
  }
  mf$formula <- formula

  ## evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## extract terms, model matrix, model response
  mt <- attr(mf, "terms")
  # mtX <- terms(formula, data = data, rhs = 1L)
  y <- Formula::model.part(formula, data = mf, lhs = 1L)
  X <- model.matrix(formula, data = mf, rhs = 1L, contrasts = contrasts)
  ## if there are subjects without any response, eliminate them
  row_no_response <- rowSums(is.na(y)) == NCOL(y)
  y <- y[!row_no_response, , drop = FALSE]
  X <- X[!row_no_response, , drop = FALSE]

  ## TODO OFFSETS
  # offset_form <- model.offset(mf)
  # offset_arg  <- offset
  # if (!is.null(offset_form) & !is.null(offset_arg)) stop("Offsets should be specified either in the formula or in the argument.")
  # if (!is.null(offset_form)) {
  #   offset <- rep(offset_form, NCOL(Y))
  # }
  # if (is.null(offset)) offset <- 0
  # if (length(offset) == 1) offset <- rep(rep.int(offset, NROW(y)), NCOL(y))
  # if (!is.null(offset_arg)) {
  #   if (length(offset_arg) != NCOL(y)) stop(sprintf("Offsets should be a list of length equal to the number of responses. Each element of the list must be a vector of length %i", NROW(Y)))
  #   if (any(sapply(offset, function(x) length(x) > 0 & length(x) != NROW(y)))) stop(sprintf("Offsets should be a list of length equal to the number of responses. Each element of the list must be a vector of length %i", NROW(Y)))
  #   offset <- c(sapply(offset_arg, function(x) if (length(X) == 0) rep.int(0, NROW(Y)) else x))
  # }

  ## weights
  # weights <- model.weights(mf)
  # if(is.null(weights)) weights <- 1
  # if(length(weights) == 1) weights <- rep.int(weights, NROW(y))
  # weights <- as.vector(weights)
  # names(weights) <- rownames(mf)

  ## response_types
  if (is.null(response_types)) {
    response_types <- rep("gaussian", NCOL(y))
  } else {
    response_types <- check_families(response_types, X, y)
  }

  ## Fit model
  fit <- mvordnorm_fit(y = y, X = X,
                      # w =  weights, # offset = offset,
                      control = control,
                      response_types = response_types)

  fit$response_types <- response_types
  fit$call <- call
  fit$mf <- mf
  fit$y <- y
  fit$X <- X
  # fit$weights <- weights
  # fit$offset <- offset
  fit$formula <- formula
  fit$control <- control
  class(fit) <- "mvordnorm"
  fit
}

#' @title Print Method for class mvordnorm
#' @description Prints regression coefficients, parameters of the error structure of class mvordnorm
#' @param x object of class \code{'mvordnorm'}
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method print mvordnorm
#' @export
print.mvordnorm <- function(x, call = FALSE, ...) {
   if(call){
     cat("\nCall:\n",
         paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
   }
   cat("Formula:\n")
   print(x$formula)
   cat("\n")
   cat("Log-likelihood at convergence:", - x$res$value, "\n")
   pars <- x$parameters
   ## ------------------------------
   cat("\nThresholds for ordinal variables:\n")
   beta <- pars[[1]]
   print(beta, ...)
   ## ------------------------------
   cat("\n Intercepts  for normal variables:\n")
   beta0n <- pars[[2]]
   print(beta0n, ...)
   ## ------------------------------
   cat("\nRegression coefficients common to all responses:\n")
   print(pars[[3]], ...)
   ## ------------------------------
   cat("\n Standard deviation parameters  for normal variables:\n")
   print(pars[[4]], ...)
   ## ------------------------------
   cat("\nCorrelation parameters:\n")
   print(pars[[5]])

}


#' @title Summary method for mvordnorm
#' @description Summary of thresholds, regression coefficients
#' and parameters of the error structure of class \code{'mvordnorm'}.
#' @param object object of class \code{'mvordnorm'}
#' @param call displays function call if \code{TRUE}
#' @param ... further arguments passed to or from other methods.
#' @method summary mvordnorm
#' @export
summary.mvordnorm <- function(object, call = FALSE, ...)
  {
   mat <- cbind.data.frame(c("nunits", nrow(object$y)),
                             c("ndim", ncol(object$y)),
                           c("logLik", round(-object$res$value,2)))
   pars <- object$parameters
   cf <- unlist(object$parameters)
   if (is.null(object$vcov)) {
     se <- NA
   } else {
     se <- sqrt(diag(object$vcov))
   }
   cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
   colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
   summary.output <- list()
   summary.output$formula <- object$formula
#   ## ------------------------
   cat("Formula: ")
   print(summary.output$formula, ...)
   cat("\n")
   summary.output$info <- format(mat, justify="right")
   write.table(summary.output$info, row.names = FALSE, col.names = FALSE, quote = FALSE)
   if(call){
     summary.output$call <- object$call
     cat("\nCall: ",
         paste(deparse(object$call, width.cutoff = getOption("deparse.cutoff")),
               sep = "\n", collapse = "\n"),
         "\n\n", sep = "")
   }
   ## ------------------------
   cat("\nThresholds:\n")
   cftheta <- cf[seq_along(pars[[1]]), ,drop = FALSE]
   summary.output$thresholds <- printCoefmat(cftheta, signif.legend = FALSE)
   ## ------------------------
   cat("\nIntercept for normals \n")
   cfbeta0 <- cf[length(pars[[1]]) + seq_along(pars[[2]]), ,drop = FALSE]
   summary.output$intercepts <- printCoefmat(cfbeta0, signif.legend = FALSE)
   ## ------------------------
   cat("\nCoefficients:\n")
   cfbeta <-  cf[length(pars[[1]]) + length(pars[[2]]) +
                   seq_along(pars[[3]]), ,drop = FALSE]
   summary.output$coefficients <- printCoefmat(cfbeta, signif.legend = FALSE)
   #---------------------------------------
   cat("\nStandard deviation of the Gaussian response variables:\n")
   cfs <-  cf[length(pars[[1]]) + length(pars[[2]]) +length(pars[[3]]) +
                seq_along(pars[[4]]), ,drop = FALSE]
   summary.output$sigmas <- printCoefmat(cfs, signif.legend = FALSE)

#   ## ------------------------------
   cat("\nCorrelation params:\n")
   cfr <-  cf[length(pars[[1]]) + length(pars[[2]]) +length(pars[[3]]) +length(pars[[4]]) +
                seq_along(pars[[5]]), ,drop = FALSE]
   summary.output$corrpars <- printCoefmat(cfr, signif.legend = FALSE)
#   ## ------------------------------

   class(summary.output) <- "summary.mvordnorm"
   return(invisible(summary.output))

}


print.summary.mvordnorm <- function(summary.output, ...){
   cat("Formula: ")
   print(summary.output$formula, ...)
   cat("\n")
   write.table(summary.output$info, row.names = FALSE, col.names = FALSE, quote = FALSE)

   cat("\nThresholds:\n")
   print(summary.output$thresholds)
   ## ------------------------
   cat("\nIntercept for normals \n")
   print(summary.output$intercepts)
   ## ------------------------
   cat("\nCoefficients:\n")
   print(summary.output$coefficients)

   #---------------------------------------
   cat("\nStandard deviation of the Gaussian response variables:\n")
   print(summary.output$sigmas)

   #   ## ------------------------------
   cat("\nCorrelation params:\n")
   print(summary.output$corrpars)
}

#' @title Control function for mvordnorm()
#' @description Control arguments are set.
# #' @param start.values.beta vector of (optional) starting values for the regression coefficients.
#' @param solver name of solver.
# #' @param scale If \code{scale = TRUE}, then for each response the corresponding covariates of \code{\link{class}} \code{"numeric"} are standardized before fitting,
# #'  i.e., by substracting the mean and dividing by the standard deviation.
#' @param se logical, should standard errors be calculated?
#' @export
mvordnorm.control <- function(se = TRUE, solver = "CG") {
#   if (is.null(solver.nlminb.control$eval.max)) solver.nlminb.control$eval.max <- 10000
#   if (is.null(solver.nlminb.control$iter.max)) solver.nlminb.control$iter.max <- 5000
#   if (is.null(solver.nlminb.control$trace)) solver.nlminb.control$trace <- 0
  list(se = se,
       solver = solver)

}
#' @title vcov of Multivariate Models with Ordinal and Gaussian Responses.
#' @description
#' \code{vcov} is a generic function which extracts the Godambe information matrix from objects of class \cr
#' \code{'mvordnorm'}.
#' @param object an object of class \code{'mvordnorm'}.
#' @param ... further arguments passed to or from other methods.
#' @method vcov mvordnorm
#' @export
vcov.mvordnorm <- function(object, ...) object$vcov

#' #' #' @title Number of units in the GLMM with mixed responses.
#' #' #' @description \code{nobs} is a generic function which extracts the number of units in the GLMM with mixed responses.
#' #' #' @param object an object of class \code{'mixoglmm'}.
#' #' #' @param ... further arguments passed to or from other methods.
#' #' #' @method nobs mixoglmm
#' #' #' @export
#' #' nobs.mixoglmm <- function(object, ...) object$dims$N
#' #'
