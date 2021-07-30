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
  Dinv <- makeDiiF(nPed)$D
    Dinv@x <- sqrt(1 / Dinv@x)
  A <- as(chol2inv(t(crossprod(makeTinv(nPed), Dinv))), "symmetricMatrix")
    A@Dimnames <- list(as.character(pedigree[, 1]),
			as.character(pedigree[, 1]))
 A
}

