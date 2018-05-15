###########################################################

# Code for the pin function is taken and modified from R code
# with the same name that was created and posted online by 
# Ian White at the University of Edinburgh

##########################################################


#' Approximate standard errors for linear functions of variance components
#' 
#' This function is similar to the pin calculations performed by the standalone
#' ASReml.  This function, written by Ian White, applies the delta method for
#' the estimation of approximate standard errors on linear functions of
#' variance components from a REML mixed model
#' 
#' Object is intended to be an asreml-R model output.
#' 
#' The formula can use \code{V1,..., Vn} to specify any one of the \code{n}
#' variance components.  These should be in the same order as they are in the
#' object (e.g., see the row order of \code{summary(object)$varcomp} for
#' asreml-R models.
#' 
#' @param object A list with at least the following elements: \code{gammas},
#'   \code{gammas.type}, and \code{ai} from a REML mixed model
#' @param transform A formula specifying the linear transformation of variance
#'   components to conduct
#'
#' @return A \code{data.frame} with row names corresponding to the operator on 
#'   the left hand side of the \code{transform} formula and the entries 
#'   corresponding to the \code{Estimate} and approximate \code{SE} of the 
#'   linear transformation.
#' @author Ian White
#' @seealso See Also as \code{\link{aiCI}}, \code{\link{aiFun}}
#' @examples
#' 
#' # Below is the heritability calculation for tait1 of the warcolak dataset
#' # Re-create the output from a basic, univariate animal model in asreml-R
#'    asrMod <- list(gammas = c(0.6383285, 1.00),
#' 	gammas.type = c(2, 1),
#' 	ai = c(0.0044461106, -0.0011754947, 0.0004424668))
#'    namevec <- c("ped(ID)!ped", "R!variance")
#'    names(asrMod[[1]]) <- names(asrMod[[2]]) <- namevec
#' 
#'    nadiv:::pin(asrMod, h2 ~ V1 / (V1 + V2))
#' @export
pin <- function (object, transform){
  pframe <- as.list(object$gammas)
  names(pframe) <- paste("V", seq(1, length(pframe)), sep = "")
  tvalue <- eval(deriv(transform[[length(transform)]], names(pframe)), pframe)
  X <- as.vector(attr(tvalue, "gradient"))
  X[object$gammas.type == 1] <- 0
  tname <- if (length(transform) == 3) transform[[2]] else ""
  n <- length(pframe)
  i <- rep(1:n, 1:n)
  j <- sequence(1:n)
  k <- 1 + (i > j)
  Vmat <- object$ai
  se <- sqrt(sum(Vmat * X[i] * X[j] * k))
 return(data.frame(row.names = tname, Estimate = tvalue, SE = se))
}




#' REML convergence checks
#' 
#' Mainly checks to ensure the variance components in a REML mixed model do not
#' change between the last two iterations more than what is allowed by the
#' tolerance value.  See details for extra check on asreml-R models.
#' 
#' Object is intended to be an asreml-R model output. NOTE, The first 3 rows
#' are ignored and thus should not be variance components from the model (e.g.,
#' they should be the loglikelihood or degrees of freedom, etc.).  Also, the
#' last column is ignored and should not be an iteration of the model (e.g., it
#' indicates the constraint).
#' 
#' The function also checks \code{object} to ensure that the output from the
#' asreml-R model does not contain a log-likelihood value of exactly 0.00.  An
#' ASReml model can sometimes fail while still returning a \code{monitor}
#' object and \code{TRUE} value in the \code{converge} element of the output.
#' This function will return \code{FALSE} if this is the case.
#' 
#' @param object A list with at least one element named: \code{monitor} (see
#'   Details)
#' @param traces Optionally, a matrix to substitute instead of the monitor
#'   element to \code{object}.  Each row corresponds to a different variance
#'   component in the model and each column is a different iteration of the
#'   likelihood calculation (column 1 is the first iterate).
#' @param tol The tolerance level for which to check against all of the changes
#'   in variance component parameter estimates
#' @param silent Optional argument to silence the output of helpful (indicating
#'   default underlying behavior) messages
#'
#' @return Returns \code{TRUE} if all variance parameters change less than the
#'   value specified by \code{tol}, otherwise returns \code{FALSE}. Also see the
#'   \code{details} section for other circumstances when \code{FALSE} might be
#'   returned.
#' @author \email{matthewwolak@@gmail.com}
#' @examples
#' 
#' # Below is the last 3 iterations from the trace from an animal model of 
#' # tait1 of the warcolak dataset.
#' # Re-create the output from a basic, univariate animal model in asreml-R
#'    tracein <- matrix(c(0.6387006, 1, 0.6383099, 1, 0.6383294, 1, 0.6383285, 1),
#' 	nrow = 2, ncol = 4, byrow = FALSE)
#'    dimnames(tracein) <- list(c("ped(ID)!ped", "R!variance"), c(6, 7, 8, 9))
#' 
#'    pcc(object = NULL, trace = tracein)
#' 
#' 
#' @export
pcc <- function(object, traces = NULL, tol = 0.01, silent = FALSE){
    if(is.null(object) & is.null(traces)){
       stop("one of object or traces must be non-NULL")
    }

    if(!is.null(object)){
       if(!silent) message("trimming the asreml monitor matrix")
       rc <- dim(object$monitor)
       traces <- object$monitor[4:rc[1], 1:(rc[2]-1)]
       if(object$loglik == 0.00 | object$converge == FALSE) return(FALSE)
    } 
    rc <- dim(traces)
    penultimate <- rc[2] - 1
    pchange <- abs(traces[, rc[2]] - traces[, penultimate]) / c(apply(traces[, penultimate:rc[2]], MARGIN = 1, FUN = max))
     return(all(pchange < tol))
}

   
