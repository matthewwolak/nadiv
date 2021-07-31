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

  nPed <- numPed(pedigree)
  N <- nrow(nPed)
  Tinv <- makeTinv(nPed)
  nPed[nPed == -998] <- N + 1
  Cout <- .C("mdiif", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 			#dam
	    as.integer(nPed[, 3] - 1),  		#sire
	    as.double(rep(0, N)),			#h
            as.double(rep(0, N)),  			#dii
            as.integer(N))   				#n

  M <- as(chol2inv(crossprod(Diagonal(x = sqrt(1 / Cout[[4]]), n = N),
        Tinv)),
      "symmetricMatrix")

    M@Dimnames <- list(as.character(pedigree[, 1]),
			as.character(pedigree[, 1]))
 M
}

