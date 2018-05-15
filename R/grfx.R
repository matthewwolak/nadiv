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
#' Supplied incidence matrices should be n-by-n symmetric matrices.  For
#' simulated random effects using design matrices, see \code{\link{drfx}}.  If
#' no incidence matrix is supplied, \code{incidence = NULL}, the function first
#' checks the environment to see if anything with the name
#' 'nadiv_prev_Mincidence' exists when checking \code{ls()}.  If so, this saved
#' version is used with a warning.  Otherwise the Identity matrix is used,
#' which assumes that all 'n' random effects are independently and identically
#' distributed (default to Identity matrix).
#' 
#' BE CAREFUL with \code{saveIncidence = TRUE} as this will save the incidence
#' matrix outside of the function environment so as to be accessed within the
#' function at a later call.  This can be useful for Monte Carlo simulation, to
#' avoid performing the cholesky decomposition on a large matrix at each
#' iteration.  Setting \code{warn = FALSE} will suppress the warnings that this
#' is occurring. DO NOT turn this warning off unless you are sure which
#' incidence matrix will be used by \code{grfx}.
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
#' @param incidence The covariance structure of the 'n' individuals
#' @param saveIncidence A logical or NULL, indicating if the cholesky
#'   decomposition of the incidence matrix should be saved/retrieved to the
#'   global environment.  SEE details below!!
#' @param output Format for the output
#' @param stdnorms Standard normal deviates to use
#' @param warn Should a warning be produced when \code{incidence = NULL}
#'   indicating that a previous incidence matrix is being used if one exists,
#'   otherwise an Identity matrix
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
#'   breedingValues <- grfx(n = 200, G = Gmat, incidence = A, saveIncidence = FALSE)
#' 
#'  # Now with a user supplied set of standard normal deviates
#'   snorms <- rnorm(nrow(warcolak[1:200,]) * ncol(Gmat))
#'   breedingValues2a <- grfx(n = 200, G = Gmat, incidence = A, stdnorms = snorms)
#'   breedingValues2b <- grfx(n = 200, G = Gmat, incidence = A, stdnorms = snorms)
#'   identical(breedingValues2a, breedingValues2b)  # TRUE
#'   var(breedingValues2a)
#'   var(breedingValues2b)
#' 
#' 
#' @export
grfx <- function(n, G, incidence = NULL, saveIncidence = FALSE, output = "matrix", stdnorms = NULL, warn = TRUE){
  d <- nrow(G)
  if(d > 1 && all(G == G[1,1])) warning("variance-covariance matrix 'G' may have caused 'chol.default(G)' error.  If so, consider subtracting 0.0001 from the covariances to make correlations < 1 or >-1")
  Mg <- as(chol(G), "dtCMatrix")
  if(is.null(incidence)){
     if(any(ls(envir = globalenv() ) == "nadiv_prev_Mincidence")){
       if(warn) warning("using previous incidence matrix")
     } else{
          if(saveIncidence){
	     nadiv_prev_Mincidence <<- Diagonal(n, 1)
             } else{
                  nadiv_prev_Mincidence <- Diagonal(n, 1)
               }
          if(warn) warning("Incidence matrix used = Identity matrix")
       }
  } else{
       if(saveIncidence){
          nadiv_prev_Mincidence <<- chol(incidence)
       } else{
            nadiv_prev_Mincidence <- chol(incidence)
         }
    }

  M <- suppressMessages(kronecker(nadiv_prev_Mincidence, Mg))
  if(is.null(stdnorms)){
     Z <- Matrix(rnorm(n*d), nrow = 1)
  } else{
       if(length(stdnorms) != n*d) stop("length(stdnorms) must be equal to 'n' times the order of 'G'")
       Z <- Matrix(stdnorms, nrow = 1)
    }
  X <- Matrix((Z %*% M)@x, ncol = d, byrow = TRUE)

 return(as(X, output))
}

