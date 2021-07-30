#' Creates the (additive) mutational effects relationship matrix
#' 
#' This returns the (additive) mutational effects relationship matrix in sparse
#'   matrix format.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' See function \code{\link{makeMinv}} for directly obtaining the inverse of
#' the (additive) mutational effects genetic relationship matrix.
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#'
#' @return Returns M, or the mutational effects relationship matrix, in sparse 
#'   matrix form.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeA}}, \code{\link{makeS}}
#' @examples
#' 
#'  makeM(Mrode2)
#' 
#' @export
makeM <- function(pedigree){

listMinv <- makeMinv(pedigree)

  nPed <- numPed(pedigree)
  N <- dim(nPed)[1]
#  nPed[nPed == -998] <- N + 1
#  f <- c(rep(0, N), -1)
#  Cout <- .C("fcoeff", PACKAGE = "nadiv",
#	    as.integer(nPed[, 2] - 1), 			#dam
#	    as.integer(nPed[, 3] - 1),  		#sire
#	    as.double(f),				#f
#           as.double(rep(0, N)),  			#dii
#          as.integer(N),   				#n
#	    as.integer(1))	#FIXME placeholder	#f missing or supplied
#  nPed[nPed == (N+1)] <- -998
  M <- as(chol2inv(t(crossprod(makeTinv(nPed),
#      Diagonal(x = sqrt(1 / Cout[[4]]), n = N)))), "symmetricMatrix")
      Diagonal(x = sqrt(1 / listMinv$dii), n = N)))), "symmetricMatrix")

    M@Dimnames <- list(as.character(pedigree[, 1]),
			as.character(pedigree[, 1]))

  
#sM <- drop0(zapsmall(solve(listMinv$Minv), digits = 5))

M
}

