check <- function(...){
  stopifnot(...)
}

check_families <- function(families, X, y) {
  if (length(families) != NCOL(y)) stop("Number of specified families does not coincide with the number of responses.")
  if (anyNA(match(names(families), colnames(y)))) stop("The names of the specified families do not fit the names of the response variables in formula.")
  if (!all(families %in% c("gaussian", "ordinal")))  stop("One or more of the specified families are not supported. Only ordinal and gaussian responses are supported by the mvordnorm package.")
  return(families)
}

## TODO implement checks for the columns
