# Generic

#' Creates components of the additive genetic relationship matrix and its inverse
#' 
#' This returns the Cholesky decomposition of the numerator relationship matrix
#' and its inverse. It can also be used to obtain coefficients of inbreeding for
#' the pedigreed population.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' The function implements an adaptation of the Meuwissen and Luo (1992)
#' algorithm (particularly, following the description of the algorithm in
#' Mrode 2005) with some code borrowed from the \code{inverseA} function by
#' Jarrod Hadfield in the \code{MCMCglmm} package. 
#' 
#' At the moment, providing the inbreeding level of individuals or the base
#' population has not been implemented. However, this argument is a placeholder
#' for now.
#' 
#' @aliases makeTinv makeTinv.default makeTinv.numPed makeDiiF makeDiiF2
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param f A numeric indicating the level of inbreeding. See Details
#' @param \dots Arguments to be passed to methods
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{Tinv }{the inverse of the Cholesky decomposition of the additive
#'       genetic relationship matrix (Ainv=Tinv' Dinv Tinv) in sparse matrix form}
#'     \item{D }{the diagonal D matrix of the A=TDT' Cholesky decomposition.
#'       Contains the variance of Mendelian sampling. Matches
#'       the order of the first/ID column of the pedigree.} 
#'     \item{f }{the individual coefficients of inbreeding for each individual 
#'       in the pedigree (matches the order of the first/ID column of the
#'       pedigree).}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeAinv}}, \code{\link{makeA}}
#' @references Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding 
#' coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.
#' 
#' Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
#' Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' 
#' @examples
#' 
#'  Tinv <- makeTinv(Mrode2)
#'  # Method for a numeric pedigree (of `nadiv` class "numPed")
#'  nPed <- numPed(Mrode2)
#'  Tinv2 <- makeTinv(nPed)
#'
#'  ########
#'  DF <- makeDiiF(Mrode2)
#'  # manually construct the inverse of the relatedness matrix `Ainv`
#'  Dinv <- DF$D  #<-- not the inverse yet, just copying the object
#'  Dinv@x <- 1 / DF$D@x  #<-- inverse of a diagonal matrix
#'  handAinv <- crossprod(Tinv, Dinv) %*% Tinv
#'    # make the A-inverse directly
#'    Ainv <- makeAinv(Mrode2)$Ainv
#'    # Compare
#'    handAinv
#'    Ainv
#'    stopifnot(all(abs((Ainv - handAinv)@x) < 1e-6))
#' 
#'
#' @export
makeTinv <- function(pedigree, ...){
  UseMethod("makeTinv", pedigree)
}

################################
# Methods:
#' @rdname makeTinv
#' @method makeTinv default
#' @export
makeTinv.default <- function(pedigree, ...){
  nPed <- numPed(pedigree)
  N <- nrow(nPed)
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv <- as(new("dtTMatrix", i = as.integer(c(nPed[, 1][dnmiss], nPed[, 1][snmiss])-1),
    j = as.integer(c(nPed[, 2][dnmiss], nPed[, 3][snmiss])-1),
    x = as.double(rep(-0.5, length(dnmiss) + length(snmiss))),
    Dim = c(N, N), Dimnames = list(as.character(pedigree[, 1]), NULL),
    uplo = "L", diag = "U"), "dtCMatrix")

 return(Tinv)
}

################################
#' @rdname makeTinv
#' @method makeTinv numPed
#' @export
makeTinv.numPed <- function(pedigree, ...){
  N <- nrow(pedigree)
  dnmiss <- which(pedigree[, 2] != -998)
  snmiss <- which(pedigree[, 3] != -998)
  Tinv <- as(new("dtTMatrix", i = as.integer(c(pedigree[, 1][dnmiss], pedigree[, 1][snmiss])-1),
    j = as.integer(c(pedigree[, 2][dnmiss], pedigree[, 3][snmiss])-1),
    x = as.double(rep(-0.5, length(dnmiss) + length(snmiss))),
    Dim = c(N, N), Dimnames = list(as.character(pedigree[, 1]), NULL),
    uplo = "L", diag = "U"), "dtCMatrix")

 return(Tinv)
}



###############################################################################
###############################################################################


#' @export
makeDiiF <- function(pedigree, ...){
  UseMethod("makeDiiF", pedigree)
}

################################
# Methods:
#' @rdname makeTinv
#' @method makeDiiF default
#' @export
makeDiiF.default <- function(pedigree, f = NULL, ...){
  renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
  nPed <- numPed(pedigree[renPed, ])
  N <- nrow(nPed)
  nPed[nPed == -998] <- N + 1
  f <- c(rep(0, N), -1)
  Cout <- .C("fcoeff", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N))   					#n
  fsOrd <- as(as.integer(renPed), "pMatrix")
  f <- Cout[[3]][t(fsOrd)@perm]
  dii <- Cout[[4]][t(fsOrd)@perm]


 return(list(D = Diagonal(x = dii, n = N),
	f = f))
}


################################
#' @rdname makeTinv
#' @method makeDiiF numPed
#' @export
makeDiiF.numPed <- function(pedigree, f = NULL, ...){
  renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
  nPed <- ronPed(pedigree, renPed)
  N <- nrow(nPed)
  nPed[nPed == -998] <- N + 1
  f <- c(rep(0, N), -1)
  Cout <- .C("fcoeff", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N))   					#n
  fsOrd <- as(as.integer(renPed), "pMatrix")
  f <- Cout[[3]][t(fsOrd)@perm]
  dii <- Cout[[4]][t(fsOrd)@perm]


 return(list(D = Diagonal(x = dii, n = N),
	f = f))
}

