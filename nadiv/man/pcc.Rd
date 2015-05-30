\name{pcc}
\alias{pcc}
\title{
REML convergence checks
}
\description{
Mainly checks to ensure the variance components in a REML mixed model do not change between the last two iterations more than what is allowed by the tolerance value.  See details for extra check on asreml-R models.}
\usage{
pcc(object, traces = NULL, tol = 0.01, silent = FALSE)
}
\arguments{
  \item{object }{
A list with at least one element named: \code{monitor} (see Details)
}
  \item{traces }{
Optionally, a matrix to substitute instead of the monitor element to \code{object}.  Each row corresponds to a different variance component in the model and each column is a different iteration of the likelihood calculation (column 1 is the first iterate).
}
  \item{tol }{
The tolerance level for which to check against all of the changes in variance component parameter estimates
}
  \item{silent }{
Optional argument to silence the output of helpful (indicating default underlying behavior) messages
}
}

\details{
Object is intended to be an asreml-R model output. NOTE, The first 3 rows are ignored and thus should not be variance components from the model (e.g., they should be the loglikelihood or degrees of freedom, etc.).  Also, the last column is ignored and should not be an iteration of the model (e.g., it indicates the constraint). 

The function also checks \code{object} to ensure that the output from the asreml-R model does not contain a log-likelihood value of exactly 0.00.  An ASReml model can sometimes fail while still returning a \code{monitor} object and \code{TRUE} value in the \code{converge} element of the output. This function will return \code{FALSE} if this is the case.
}
\value{
Returns \code{TRUE} if all variance parameters change less than the value specified by \code{tol}, otherwise returns \code{FALSE}. Also see the \code{details} section for other circumstances when \code{FALSE} might be returned.
}
\author{
\email{matthewwolak@gmail.com}
}

\examples{
# Below is the last 3 iterations from the trace from an animal model of 
# tait1 of the warcolak dataset.
# Re-create the output from a basic, univariate animal model in asreml-R
   tracein <- matrix(c(0.6387006, 1, 0.6383099, 1, 0.6383294, 1, 0.6383285, 1),
	nrow = 2, ncol = 4, byrow = FALSE)
   dimnames(tracein) <- list(c("ped(ID)!ped", "R!variance"), c(6, 7, 8, 9))

   pcc(object = NULL, trace = tracein)

}

