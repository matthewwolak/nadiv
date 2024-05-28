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
#' The inbreeding level of individuals can be provided instead of calculated.
#' \code{f} must be a vector that is the same length as individuals in the
#' pedigree. Supplied coefficients of inbreeding are used instead of being 
#' calculated until a \code{NA} is encountered in the vector. From this position
#' on, then coefficients of inbreeding are calculated and replace entries in 
#' \code{f}. This can be used, for example, to calculate coefficients of
#' inbreeding for later generations when coefficients of inbreeding in the
#' previous generations have already been calculated. To specify an average
#' coefficient of inbreeding for the base population, modify the pedigree to
#' include a single phantom parent and specify this individual's non-zero
#' coefficient of inbreeding in \code{f} with the rest of the terms as NA.
#' 
#' @aliases makeT makeT.default makeT.numPed
#'   makeTinv makeTinv.default makeTinv.numPed makeDiiF
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param genCol An integer value indicating the generation up to which the
#'   \code{T} matrix is to be created (corresponding to columns of the lower
#'   triangle \code{T} matrix). The first generation is numbered 0, default is
#'   all generations.
#' @param f A numeric vector indicating the level of inbreeding. See Details
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
#'  # supply previous generation coefficients of inbreeding (f)
#'  ## to keep from re-calculating their f when analyzing subsequent generations
#'  DF <- makeDiiF(Mrode2[, 1:3])
#'  Mrode2$gen <- genAssign(Mrode2)
#'  Mrode2$f_full <- DF$f
#'  Mrode2$f_in <- with(Mrode2, c(f_full[gen <= 1], rep(NA, sum(gen > 1))))
#'  DF2 <- makeDiiF(Mrode2[, 1:3], f = Mrode2$f_in) 
#'  stopifnot(identical(DF, DF2))
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
    uplo = "L", diag = "U"), "CsparseMatrix")

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
    uplo = "L", diag = "U"), "CsparseMatrix")

 return(Tinv)
}



###############################################################################
###############################################################################
#' @export
makeT <- function(pedigree, genCol = NULL, ...){
  UseMethod("makeT", pedigree)
}

################################
# Methods:
#' @rdname makeTinv
#' @method makeT default
#' @export
makeT.default <- function(pedigree, genCol = NULL, ...){

  gen <- genAssign(pedigree)
    if(is.null(genCol)) genCol <- max(gen)
  pedOrd <- order(gen, pedigree[, 2], pedigree[, 3])
  nPed <- numPed(pedigree[pedOrd, 1:3])
  ogen <- gen[pedOrd]
  
  N <- dim(nPed)[1]
  n0 <- sum(ogen == 0)
  ncol <- sum(ogen <= genCol)
  # size of upper part of Tcol (lower-triangle, incl. diagonal, matrix)
  nlt <- ncol * (ncol + 1) / 2
  # size of rectangle that is rest of ncols below upper trianglular portion
  nrect <- (N - ncol) * ncol


  Cout <- .C("Trow", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
            as.double(c(rep(1.0, n0),
              rep(0.0, (nlt - n0) + nrect))),  		        #xTrow
	    as.integer(c(seq(n0),
	      rep(0, (nlt - n0) + nrect)) - 1),			#iTrow
	    as.integer(c((1:n0),
	      rep(n0 + 1, N - n0 + 1)) - 1), 			#pTrow
            as.integer(c(ncol, n0, N)))				


  nnz <- Cout[[5]][N+1]
  Trow <- sparseMatrix(i = as.integer(Cout[[4]][1:nnz]),
    p = as.integer(Cout[[5]]),
    x = as.double(Cout[[3]][1:nnz]),
    index1 = FALSE, dims = c(ncol, N), symmetric = FALSE,
    dimnames = list(NULL, NULL))

  fsOrd <- as(as.integer(pedOrd), "pMatrix")
  if(ncol == N){
    T <- crossprod(fsOrd, crossprod(Trow, fsOrd))
  } else{
      T <- crossprod(fsOrd,
        crossprod(rbind(Trow,
          sparseMatrix(i = as.integer(seq(N-ncol) - 1),
            j = as.integer(seq(ncol+1, N, 1) - 1),
            x = rep(1.0, N-ncol),
            index1 = FALSE, dims = c(N-ncol, N), symmetric = FALSE,
            dimnames = list(NULL, NULL))),
        fsOrd))[, gen <= genCol]
    }          
 return(T)
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
  if(fmiss <- missing(f)){
    f <- rep(0, N)
  } else{
      f <- f[renPed]
      nafInd <- which(is.na(f))
        f[nafInd] <- 0.0
      fmiss <- min(nafInd)     
    }
  f <- c(f, -1)
  Cout <- .C("diif", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 		#dam
	    as.integer(nPed[, 3] - 1),  	#sire
	    as.double(f),			#f
            as.double(rep(0, N)),  		#dii
            as.integer(N),   			#n
            as.integer(0),			#number genetic groups
	    as.integer(fmiss - 1))		#first f to calculate (not supplied)
  fsOrd <- as(as.integer(renPed), "pMatrix")
  f <- Cout[[3]][invPerm(fsOrd@perm)]
  dii <- Cout[[4]][invPerm(fsOrd@perm)]


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
  if(fmiss <- missing(f)){
    f <- rep(0, N)
  } else{
      f <- f[renPed]
      nafInd <- which(is.na(f))
        f[nafInd] <- 0.0
      fmiss <- min(nafInd)     
    }
  f <- c(f, -1)
  Cout <- .C("diif", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 		#dam
	    as.integer(nPed[, 3] - 1),  	#sire
	    as.double(f),			#f
            as.double(rep(0, N)),  		#dii
            as.integer(N),   			#n
            as.integer(0),			#number genetic groups
	    as.integer(fmiss))			#first f to calculate (not supplied)
  fsOrd <- as(as.integer(renPed), "pMatrix")
  f <- Cout[[3]][invPerm(fsOrd@perm)]
  dii <- Cout[[4]][invPerm(fsOrd@perm)]


 return(list(D = Diagonal(x = dii, n = N),
	f = f))
}

