#' Simulated design random effects
#' 
#' This function simulates effects for random terms in a linear mixed model
#' based on design matrices. The intended purpose is for simulating
#' environmental effects from a pedigree.
#' 
#' If G = x, where 'x' is a single number, then 'x' should still be specified
#' as a 1-by-1 matrix (e.g., \code{matrix(x)}).  Note, the G-matrix should
#' never have a structure which produces a correlation exactly equal to 1 or
#' -1.  Instead, covariances should be specified so as to create a correlation
#' of slightly less than (greater than) 1 (-1).  For example: 0.9999 or
#' -0.9999.
#' 
#' @param G The variance-covariance matrix to model the effects after
#' @param fac A character indicating the factor in \code{dataf} with which to
#'   construct the design matrix
#' @param dataf A dataframe with \code{fac} in it
#' @param ...  Arguments to be passed to the internal use of \code{\link{grfx}}
#'
#' @return \item{fx }{A matrix with 'd' columns of random effects} \item{Z }{A
#' design matrix (of the format 'Matrix') from which the random effects in
#'   \code{fx} were assigned }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{grfx}}
#' @examples
#' 
#' # Create maternal common environment effects for 2 traits
#' # with perfectly correlated effects
#'   Gmat <- matrix(c(10, 7.071, 7.071, 5), 2, 2)
#'   cfx <- drfx(G = Gmat, fac = "Dam", dataf = warcolak[1:200, ])
#' 
#' 
#' @export
drfx <- function(G, fac, dataf, ...){
   dataf[, fac] <- as.factor(dataf[, fac])
   d <- nrow(G)
   if(all(G == G[1,1]) & d > 1){
      warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
   }
   Z <- sparse.model.matrix(as.formula(paste0("~", fac, " - 1")), dataf)
   M <- grfx(n = ncol(Z), G = G, incidence = Diagonal(ncol(Z)), ...)
   fx <- sapply(seq.int(d), FUN = function(c){ (Z %*% M[, c])@x})
 return(list(fx = fx, Z = Z))
}

