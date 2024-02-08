#' Simulated genetic random effects
#' 
#' This function simulates effects for random terms in a linear mixed model
#' based on relatedness matrices. The intended purpose is for simulating
#' genetic and environmental effects from a pedigree.
#' 
#' The total number of effects simulated will be n*d, where d is the number of
#' columns in the 'G' matrix. The standard normal deviates can be supplied
#' instead of generated within the function when \code{stdnorms != NULL}. The
#' length of this vector must be \code{n*nrow(G)}.
#' 
#' Supplied incidence matrices should be n-by-n symmetric matrices or cholesky
#' factorizations that resulted from a call to \code{Matrix::Cholesky()}.  For
#' simulated random effects using design matrices, see \code{\link{drfx}}.  If
#' no incidence matrix is supplied, \code{incidence = NULL}, the Identity matrix
#' is used, which assumes that all 'n' random effects are independently and
#' identically distributed (default to Identity matrix).
#' 
#' See examples for how to make and use a Cholesky factorized incidence matrix,
#' for instance in a Monte Carlo simulation. Whether such an approach results
#' in performance of speed improvements within the Monte Carlo simulation, by
#' avoiding a Cholesky decomposition of a large matrix at each iteration, has
#' not been tested. Setting \code{warn = FALSE} will suppress the warnings that
#' the function is assuming a Cholesky factorization is contained in the object
#' supplied to the \code{incidence} argument. Currently, Cholesky factorizations
#' must inheriti from the class \dQuote{CHMfactor}.
#' 
#' If G = x, where 'x' is a single number, then 'x' should still be specified
#' as a 1-by-1 matrix (e.g., \code{matrix(x)}).  Note, the G-matrix should
#' never have a structure which produces a correlation exactly equal to 1 or
#' -1.  Instead, covariances should be specified so as to create a correlation
#' of slightly less than (greater than) 1 (-1).  For example: 0.9999 or
#' -0.9999.
#' 
#' @param n The number of individuals for which to simulate effects
#' @param G The variance-covariance matrix to model the effects after
#' @param incidence A matrix of the covariance structure of the 'n' individuals
#'   or the Cholesky factorization of class \code{CHMfactor} for this structure.
#' @param output Format for the output
#' @param stdnorms Standard normal deviates to use
#' @param warn Should a warning message be produced when the function interprets
#'   what to do based on the object class supplied to \code{incidence}
#'
#' @return The random effects coerced to be in the format specified by output.
#' The default is a "matrix".
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[MCMCglmm]{MCMCglmm}}, \code{\link{drfx}},
#' \code{\link{makeA}}, \code{\link{makeAA}}, \code{\link{makeD}},
#' \code{\link{makeDomEpi}}, \code{\link{makeDsim}}, \code{\link{makeS}}
#' @examples
#' 
#' # Create additive genetic breeding values for 2 uncorrelated traits
#' # with different additive genetic variances
#'   A <- makeA(warcolak[1:200, 1:3])
#'   Gmat <- matrix(c(20, 0, 0, 10), 2, 2)
#'   breedingValues <- grfx(n = 200, G = Gmat, incidence = A)
#' 
#'  # Now with a user supplied set of standard normal deviates
#'   snorms <- rnorm(nrow(warcolak[1:200,]) * ncol(Gmat))
#'   breedingValues2a <- grfx(n = 200, G = Gmat, incidence = A, stdnorms = snorms)
#'   breedingValues2b <- grfx(n = 200, G = Gmat, incidence = A, stdnorms = snorms)
#'   identical(breedingValues2a, breedingValues2b)  #<-- TRUE
#'   var(breedingValues2a)
#'   var(breedingValues2b)
#' 
#'  # User supplied Cholesky factorization of the incidence matrix from above
#'   cA <- Cholesky(A, LDL = FALSE, super = FALSE) 
#'     inherits(cA, "CHMfactor")  #<-- TRUE
#'   breedingValues3 <- grfx(n = 200, G = Gmat, incidence = cA, stdnorms = snorms)
#'   all.equal(breedingValues2a, breedingValues3)  #<-- TRUE
#' @export
grfx <- function(n, G, incidence = NULL, output = "matrix", stdnorms = NULL,
		  warn = TRUE){
  d <- nrow(G)
  if(d > 1 && all(G == G[1,1])){
    warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
  }
  
  Mg <- as(as(chol(G), "triangularMatrix"), "CsparseMatrix")

  if(is.null(incidence)){
    chol_incidence <- Diagonal(n, 1)
    if(warn) warning("Incidence matrix used = Identity matrix")
  } else{
      if(inherits(incidence, "CHMfactor")){
        if(warn) warning("Object given to incidence inherits from class 'CHMfactor'. Object in incidence being used as a Cholesky factor of an incidence matrix")
        chol_incidence <- Reduce("%*%",
               expand2(incidence, LDL = FALSE)[c("P1.", "L.", "P1")])
      } else chol_incidence <- chol(incidence)
    }  

  M <- kronecker(chol_incidence, Mg)
  if(is.null(stdnorms)){
     Z <- Matrix(rnorm(n*d), nrow = 1)
  } else{
      if(length(stdnorms) != n*d){
        stop("length(stdnorms) must be equal to 'n' times the order of 'G'")
      }
      Z <- Matrix(stdnorms, nrow = 1)
    }
  X <- Matrix((Z %*% M)@x, ncol = d, byrow = TRUE)

 return(as(X, output))
}

