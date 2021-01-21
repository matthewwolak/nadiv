#' Creates the additive genetic relationship matrix
#' 
#' This returns the additive relationship matrix in sparse matrix format.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' Used as a support function to \code{\link{makeD}}.
#' 
#' See function \code{\link{makeAinv}} for directly obtaining the inverse of
#' the additive genetic relationship matrix.
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#'
#' @return Returns A, or the numerator relationship matrix, in sparse 
#'   matrix form.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeD}}, \code{\link{makeS}}
#' @examples
#' 
#'  makeA(Mrode2)
#' 
#' @export
makeA <- function(pedigree)
{
  nPed <- numPed(pedigree)
  N <- dim(nPed)[1]
  nPed[nPed == -998] <- N + 1
  #TODO: implement functionality/argument to supply f-coeffs (see `makeDiiF()`)
  f <- c(rep(0, N), -1)
  Cout <- .C("fcoeff", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 			#dam
	    as.integer(nPed[, 3] - 1),  		#sire
	    as.double(f),				#f
            as.double(rep(0, N)),  			#dii
            as.integer(N),   				#n
	    as.integer(1))	#FIXME placeholder	#f missing or supplied
  nPed[nPed == (N+1)] <- -998
  A <- as(chol2inv(t(crossprod(makeTinv(nPed), Diagonal(x = sqrt(1 / Cout[[4]]), n = N)))), "symmetricMatrix")
    A@Dimnames <- list(as.character(pedigree[, 1]),
			as.character(pedigree[, 1]))
 A
}

