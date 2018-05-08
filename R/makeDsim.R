#' Create the dominance genetic relationship matrix through an iterative
#' (simulation) process
#' 
#' Alleles are explicitly traced through a pedigree to obtain coefficients of
#' fraternity between pairs of individuals (the probability of sharing both
#' alleles identical by descent) - for either autosomes or sex chromosomes.
#' This is accomplished in an iterative process to account for the various
#' routes by which an allele will progress through a pedigree due to Mendelian
#' sampling at either autosomes or sex chromosomes. The autosomal case is an
#' implementation of the simulation approach of Ovaskainen et al. (2008).
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' \code{parallel} = TRUE should only be used on Linux or Mac operating systems
#' (i.e., not Windows).
#' 
#' Ovaskainen et al. (2008) indicated that the method of calculating the D
#' matrix (see \code{\link{makeD}}) is only an approximation.  They proposed a
#' simulation method that is implemented here.  This should be more
#' appropriate, especially when inbreeding occurs in the pedigree.
#' 
#' The objects \code{listDsim} and \code{listSdsim} will list both the
#' approximate values (returned from \code{\link{makeD}} or
#' \code{\link{makeSd}}) as well as the simulated values.  If \code{calcSE} is
#' TRUE, these values will be listed in \code{listDsim} or \code{listSdsim}.
#' 
#' @aliases makeDsim makeSdsim
#' @param pedigree A pedigree with columns organized: ID, Dam, Sire. For use
#'   with \code{makeSdsim}, a fourth column indicates the sex of each individual
#'   in the pedigree.
#' @param N The number of times to iteratively trace alleles through the
#'   pedigree
#' @param heterogametic Character indicating the label corresponding to the
#'   heterogametic sex used in the "Sex" column of the pedigree
#' @param DosageComp A character indicating which model of dosage compensation.
#'   If \code{NULL} then the \dQuote{ngdc} model is assumed.
#' @param parallel A logical indicating whether or not to use parallel
#'   processing. Note, this may only be available for Mac and Linux operating
#'   systems.
#' @param ncores The number of cpus to use when constructing the dominance
#'   relatedness matrix. Default is all available.
#' @param invertD,invertSd A logical indicating whether or not to invert the D
#'   or Sd matrix
#' @param calcSE A logical indicating whether or not the standard errors for
#'   each coefficient of fraternity should be calculated
#' @param returnA,returnS Logical, indicating if the numerator relationship
#'   matrix (A or S) should be stored and returned.
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{A,S }{the A or S matrix in sparse matrix form}
#'     \item{D,Sd }{the approximate D or Sd matrix in sparse matrix form}
#'     \item{logDetD,logDetSd }{the log determinant of the approximate D or
#'       approximate Sd matrix}
#'    \item{Dinv,Sdinv }{the inverse of the approximate D or approximate Sd
#'      matrix in sparse matrix form}
#'    \item{listDinv,listSdinv }{the three column form of the non-zero elements 
#'      for the inverse of the approximate D matrix or the inverse of the
#'      approximate Sd matrix}
#'    \item{Dsim,Sdsim }{the simulated D or Sd matrix in sparse matrix form}
#'    \item{logDetDsim,logDetSdsim }{the log determinant of the simulated D or
#'      simulated Sd matrix}
#'    \item{Dsiminv,Sdsiminv }{the inverse of the simulated D or simulated Sd
#'      matrix in sparse matrix form}
#'    \item{listDsim,listSdsim }{the three column form of the non-zero and
#'      non-self elements for the simulated D or simulated Sd matrix}
#'    \item{listDsiminv,listSdsiminv }{the three column form of the non-zero
#'      elements for the inverse of the simulated D or the inverse of the 
#'      simulated Sd matrix}
#'   }
#' @note This simulation can take a long time for large pedigrees (a few
#'   thousand and higher) and large values of \code{N} (one thousand and 
#'   higher). If unsure, it is advisable to start with a lower \code{N} and 
#'   gradually increase to obtain a sense of the time required to execute a 
#'   desired \code{N}.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeD}}, \code{\link{makeSd}}
#' @references Ovaskainen, O., Cano, J.M., & Merila, J. 2008. A Bayesian
#' framework for comparative quantitative genetics. Proceedings of the Royal
#' Society B 275, 669-678.
#' @examples
#' 
#'   simD <- makeDsim(Mrode9, N = 1000, parallel = FALSE,
#' 		invertD = TRUE, calcSE = TRUE)$listDsim
#' 
#'   simSd <- makeSdsim(FG90, heterogametic = "0", N = 1000, parallel = FALSE,
#' 		invertSd = TRUE, calcSE = TRUE)$listSdsim
#' @export
makeDsim <- function(pedigree, N, parallel = FALSE, ncores = getOption("mc.cores", 2L), invertD = TRUE, calcSE = FALSE, returnA = FALSE){

  approxD <- makeD(pedigree, parallel = parallel, ncores = ncores, invertD = invertD, returnA = returnA)
  lapproxD <- summary(approxD$D)

  nPed <- numPed(pedigree)
  n <- dim(pedigree)[1]
  alleles <- matrix(as.integer(-998), nrow = n, ncol=2) 
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  ndfounders <- length(dfounders)
  
  alleles[dfounders, 1] <- as.integer(seq(1, ndfounders, 1)) 
  alleles[sfounders, 2] <- as.integer(seq(ndfounders+1, (ndfounders + length(sfounders)), 1))
  dalleles <- rep(alleles[, 1], each = N)
  salleles <- rep(alleles[, 2], each = N)
  
  cat("making Dsim ...")

  Cout <- .C("dsim", PACKAGE = "nadiv",
	as.integer(dalleles),
	as.integer(salleles),
	as.integer(N),
	as.integer(n),
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1),
	as.integer(approxD$D@i),
	as.integer(approxD$D@p),
	as.integer(rep(0, length(approxD$D@i))))

  lapproxD$simD <- Cout[[9]] / N
  lapproxD <- lapproxD[which(lapproxD[, 4] != 0), ]
  listDsim <- NULL
  if(calcSE) {
     lapproxD$Dse <- vapply(lapproxD$simD, FUN = function(x, N){(sqrt(x * (1 - x))) / sqrt(N)}, FUN.VALUE = vector("numeric", 1), N)
     listDsim <- lapproxD
  } 

  Dsim.row<- lapproxD[,1]
  Dsim.col<- lapproxD[,2]
  Dsim.x<- lapproxD[,4]
  order.index<-order(Dsim.col + Dsim.row/(n+1), decreasing=FALSE)

  Dsim<-Matrix(0, n, n, dimnames = list(as.character(pedigree[, 1]), NULL))
  Dsim@uplo<-"U"
  Dsim@i<-as.integer(Dsim.row[order.index]-1)
  Dsim@p<-as.integer(c(match(1:n, Dsim.col[order.index]), length(order.index)+1)-1)
  Dsim@x<-Dsim.x[order.index]
  diag(Dsim) <- diag(approxD$D)
  cat(".done", "\n")
  logDetDsim <- determinant(Dsim, logarithm = TRUE)$modulus[1]
  
  if(invertD){
    cat("inverting Dsim ...")
    Dsiminv <- solve(Dsim)
      Dsiminv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    cat(".done", "\n")
    listDsiminv <- sm2list(Dsiminv, rownames = pedigree[,1], colnames = c("row", "column", "simDinverse"))
    Dsim <- as(Dsim, "dgCMatrix")
    return(list(A = approxD$A, D = approxD$D, logDetD = approxD$logDet, Dinv = approxD$Dinv, listDinv = approxD$listDinv, Dsim = Dsim, logDetDsim = logDetDsim, Dsiminv = Dsiminv, listDsim = listDsim, listDsiminv = listDsiminv))
  } else{
      return(list(A = approxD$A, D = approxD$D, logDetD = approxD$logDet, Dsim = Dsim, logDetDsim = logDetDsim, listDsim = listDsim))
    } 

}


