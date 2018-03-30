#' Function used in the \code{proLik} function to produce a profile likelihood
#' for a variance component
#' 
#' Given a model object from \code{asreml} and a range of estimates of the
#' parameter, the function will supply the likelihood ratio test statistic for
#' the comparison of the full model to one where the parameter of interest is
#' constrained.
#' 
#' Used internally in the \code{proLik} function to call \code{constrainFun}
#' 
#' @param x section of all parameter values to analyze
#' @param parameters a value for which the log Likelihood of a model is to be
#'   calculated
#' @param full the full model \code{asreml} object
#' @param fm2 starting values for the full model
#' @param comp which variance component to constrain
#' @param G logical indicating if the component is part of the G structure
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso See Also \code{\link{proLik}}, \code{\link{constrainFun}}
#' @export
parConstrainFun <- function(x, parameters, full, fm2, comp, G)
  {
   vapply(parameters[min(x):max(x)], FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full, fm2 = fm2, comp = comp, G = G)
 }

