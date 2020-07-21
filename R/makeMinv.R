#' Create the inverse (additive) mutational effects relationship matrix
#' 
#' Returns the inverse of the (additive) mutational effects relationship matrix.
#' It can also be used to obtain components needed for the calculations in the
#' underlying algorithm.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' Note the assumption under the infinitesimal model, that mutation has essentially
#' zero probability of affecting an inbred locus (hence removing inbred
#' identity-by-descent), however, mutations may themselves be subject to
#' inbreeding (Wray 1990).
#' 
#' @aliases makeMinv
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param \dots Arguments to be passed to methods
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{Minv }{the inverse of the (additive) mutational effects
#'       relationship matrix in sparse matrix form}
#'     \item{listMinv }{the three column list of the non-zero elements for the 
#'       inverse of the (additive) mutational effects relationship matrix.
#'       \code{attr(*, "rowNames")} links the integer for rows/columns to the ID
#'       column from the pedigree.}
#'     \item{h }{the amount by which segregtation variance is reduced by
#'       inbreeding. Similar to the individual coefficients of inbreeding (f)
#'       derived during the construction of the inverse numerator relatedness matrix.
#'       in the pedigree (matches the order of the first/ID column of the
#'       pedigree).}
#'     \item{logDet }{the log determinant of the M matrix}
#'     \item{dii }{the (non-zero) elements of the diagonal D matrix of the M=TDT'
#'       decomposition. Contains the variance of Mendelian sampling. Matches
#'       the order of the first/ID column of the pedigree. Note Wray (1990) and
#'       Casellas and Medrano (2008) algorithms use \code{v=sqrt(dii)}.} 
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @references Casellas, J. and J.F. Medrano. 2008. Within-generation mutation
#' variance for litter size in inbred mice. Genetics. 179:2147-2155. 
#'
#' Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding
#' coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.
#' 
#' Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
#' Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' 
#' Wray, N.A. 1990. Accounting for mutation effects in the additive genetic
#' variance-covariance matrix and its inverse. Biometrics. 46:177-186.
#' @examples
#' 
#'  ##  Example pedigree from Wray 1990
#'    Wray90 <- data.frame(id = seq(8),
#'	dam = c(NA, NA, 1, 1, 4, 4, 6, 7),
#'	sire = c(NA, NA, 2, 2, 3, 2, 5, 5),
#'	time = c(0, 0, 1, 1, 2, 2, 3, 4))
#'    Mout <- makeMinv(Wray90)
#' 
#' @export
makeMinv <- function(pedigree, ...){
#TODO / FIXME; turned off for now
  renPed <- seq(nrow(pedigree))
#  renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3],
#    na.last = FALSE)
  nPed <- numPed(pedigree[renPed, ])
  N <- nrow(nPed)
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
    p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
    index1 = FALSE, dims = c(N, N), symmetric = FALSE,
    dimnames = list(as.character(nPed[, 1]), NULL))

  Minv <- t(crossprod(sTinv)) # transpose gives lower triangle
  #TODO First checks to see if individual k has same dam and sire as k-1, if so then just assigns k-1's f 
  nPed[nPed == -998] <- N + 1
  h <- c(rep(0, N), -1)
  Cout <- .C("minv", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(h),					#h (f)
            as.double(rep(0, N)),  				#dii/v
            as.integer(N),   					#n
            as.double(rep(0, length(Minv@i))),  			#xMinv
	    as.integer(Minv@i), 				#iMinv
	    as.integer(Minv@p), 				#pMinv
	    as.double(rep(0, 1)))				#logDet of M
  Minv <- as(Minv, "dsCMatrix")
  Minv@x <- Cout[[6]]
  fsOrd <- as(as.integer(renPed), "pMatrix")
  Minv <- as(crossprod(fsOrd, Minv) %*% fsOrd, "dgCMatrix")
    Minv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
  h <- Cout[[3]][t(fsOrd)@perm][1:N]
  dii <- Cout[[4]][t(fsOrd)@perm][1:N]

 return(list(Minv = Minv,
	listMinv = sm2list(Minv, rownames = rownames(Minv),
	  colnames = c("row", "column", "Minv")),
	h = h,
	logDet = Cout[[9]],
	dii = dii))
}






