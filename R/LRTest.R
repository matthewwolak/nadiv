# Likelihood Ratio Test:
# Tests the hypothesis that the reduced model offers a better fit
# Helpful reading: section 6.4.1.1 of Bolker 2008. pp. 189-194
#NOTE: sometimes ASReml-R returns positive log-likelihoods, other times negative log-likelihoods. My experience so far is that positive values are returned when there are often boundary parameters or the model is having trouble fitting the data. 
#	 ***!!!ASSUMPTION:*** 
##  	when both are positive, the log-likelihood is being maximized (i.e., the greater value is a better fit).


#' log-Likelihood Ratio Test
#' 
#' Test the null hypothesis that the two models fit the data equally well.
#' 
#' Boundary correction should be applied if the parameter that is dropped from
#' the full model was on the boundary of its parameter space. In this instance,
#' the distribution of the log-likelihood ratio test statistic is approximated
#' by a mix of chi-square distributions (Self and Liang 1987). A \code{TRUE}
#' value will implement the boundary correction for a one degree of freedom
#' test. This is equivalent to halving the p-value from a test using a
#' chi-square distribution with one degree of freedom (Dominicus et al. 2006).
#' 
#' Currently, the test assumes that both log-likelihoods are negative or both
#' are positive and will stop if they are of opposite sign. The interpretation
#' is that the model with a greater negative log-likelihood (closer to zero) or
#' greater positive log-likelihood provides a better fit to the data.
#' 
#' @usage LRTest(full, reduced, df = 1, boundaryCorrection = FALSE)
#' @param full A numeric variable indicating the log-likelihood of the full
#' model
#' @param reduced A numeric variable indicating the log-likelihood of the
#' reduced model
#' @param df The number of degrees of freedom to use, representing the
#' difference between the full and reduced model in the number of parameters
#' estimated
#' @param boundaryCorrection A logical argument indicating whether a boundary
#' correction under one degree of freedom should be included. If the parameter
#' that is dropped from the reduced model is estimated at the boundary of its
#' parameter space in the full model, the boundary correction is often
#' required. See Details for more.
#' @return \item{lambda }{a numeric log-likelihood ratio test statistic}
#' \item{Pval }{a numeric p-value given the \code{lambda} tested against a
#' chi-squared distribution with the number of degrees of freedom as specified.
#' May have had a boundary correction applied.} \item{corrected.Pval }{a
#' logical indicating if the p-value was derived using a boundary correction.
#' See \code{Details}}
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{constrainFun}}
#' @references Self, S. G., and K. Y. Liang. 1987. Asymptotic properties of
#' maximum likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. Journal of the American Statistical Association 82:605-610.
#' 
#' Dominicus, A., A. Skrondal, H. K. Gjessing, N. L. Pedersen, and J. Palmgren.
#' 2006. Likelihood ratio tests in behavioral genetics: problems and solutions.
#' Behavior Genetics 36:331-340.
#' @examples
#' 
#' # No boundary correction
#' (noBC <- LRTest(full = -2254.148, reduced = -2258.210,
#' 	df = 1, boundaryCorrection = FALSE))
#' # No boundary correction
#' (withBC <- LRTest(full = -2254.148, reduced = -2258.210,
#' 	df = 1, boundaryCorrection = TRUE))
#' stopifnot(noBC$Pval == 2*withBC$Pval)
#' 
#' @export LRTest
LRTest <- function(full, reduced, df = 1, boundaryCorrection = FALSE){
    if(sign(full) != sign(reduced)){
       stop("Signs of the log-likelihoods are opposite - or 1 log-likelihood is zero...don't know what to do")
    }
    # positive log-likelihoods: better fit has higher log-likelihood (more positive)
    if(sign(full) > 0){
      lambda <- 2*(full - reduced)
      warning("Positive log-likelihoods:\nASSUMING full model has greater log-likelihood if it fits the data better than the reduced model")
    }
    # negative log-likelihoods: better fit has greater log-likelihood (less negative)
    if(sign(full) < 0){
       lambda <- 2*(full - reduced)
    }

    if(boundaryCorrection & df == 1){
       lrtP <- 0.5*(pchisq(lambda, df = df, lower.tail = FALSE))
    } else{
         lrtP <- pchisq(lambda, df = df, lower.tail = FALSE)
      }
         
 return(list(lambda = lambda, Pval = lrtP, corrected.Pval = boundaryCorrection))
}

