#' Create the dominance genetic relationship matrix
#' 
#' Given a pedigree, the matrix of coefficients of fraternity are returned -
#' the D matrix for autosomes and the Sd matrix for sex chromosomes. Note,
#' inbreeding is not directly incorporated into the calculation of the
#' coefficients (see Details). Functions will return the inverses of the D and
#' Sd matrices by default, otherwise this operation can be skipped if desired.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' There exists no convenient method of obtaining the inverse of the dominance
#' genetic relatedness matrix (or the D matrix itself) directly from a pedigree
#' (such as for the inverse of A, i.e., Quaas (1995)). Therefore, these
#' functions computes the coefficient of fraternity (Lynch and Walsh, 1998) for
#' every individual in the pedigree with a non-zero additive genetic
#' relatedness in the case of autosomes (\code{makeD}) or for the homogametic
#' sex only in the case of sex chromosomes (\code{makeSd}, because the
#' heterogametic sex has only one copy of the shared sex chromosome and
#' therefore cannot express dominance allelic interactions).
#' 
#' The coefficients of fraternity are only approximations that assume no
#' inbreeding. The algorithm used here, however, incorporates inbreeding into
#' the calculation of coefficients of coancestry (using `makeA()`) that are
#' used to calculate coefficients of fraternity. Similarly, the diagonals of
#' the D and Sd matrices are corrected for inbreeding. Meaning, the diagonals
#' of D and Sd are (1-f) so that the overall dominance genetic variance is
#' equal to (1-f)V_D, where f is the coefficient of inbreeding and V_D is
#' dominance genetic variance. This is interpreted as the amount of dominance
#' genetic variance that would be expected if the allele frequencies in the
#' inbred population were representative of a non-inbred, randomly mating
#' population (Shaw et al. 1998; Wolak and Keller 2014). Note, the construction
#' of the D matrix is more computationally demanding (in computing time and
#' space requirements) than is the construction of A. This is possibly also the
#' case for construction of Sd in comparison to the S matrix.
#' 
#' To overcome the computational difficulties, this function can enable
#' parallel processing (see package \code{parallel} included in the R
#' distribution) to speed up the execution. Note this is not be possible on
#' Windows (See \code{parallel} documentation for further information),
#' therefore \code{parallel} = TRUE should only be used on Linux or Mac
#' operating systems (i.e., not Windows). The default is to use the maximum
#' number of cpus available to the machine, but this can be restricted by
#' indicating the number desired in the argument \code{ncores}. Setting up the
#' multi-processing takes some overhead, so no real advantage is gained for
#' small pedigrees. Also, since all processes are sharing a fixed amount of
#' RAM, very large pedigrees using many processes in parallel may not be
#' feasible due to RAM restrictions (i.e., if each process needs "n" amount of
#' RAM to run, then \code{ncores} should be set to = total RAM/n). Otherwise
#' the machine can become overworked.
#' 
#' Note, for very large pedigrees \code{returnA} or \code{returnS} should be
#' set to FALSE to avoid drastically increasing the memory requirements while
#' making D or Sd, respectively. When this occurs, 'NULL' is returned for the
#' element of 'A' in the output of \code{makeD} or for the element of 'S' in
#' the output of \code{makeSd}.
#' 
#' @aliases makeD makeSd
#' @param pedigree A pedigree with columns organized: ID, Dam, Sire. For use
#'   with \code{makeSd}, a fourth column indicates the sex of each individual in
#'   the pedigree.
#' @param heterogametic Character indicating the label corresponding to the
#'   heterogametic sex used in the "Sex" column of the pedigree
#' @param DosageComp A character indicating which model of dosage compensation.
#'   If \code{NULL} then the \dQuote{ngdc} model is assumed.
#' @param parallel Logical, indicating whether computation should be run on
#'   multiple processors at once. See details for considerations.
#' @param ncores Number of cores to use when parallel = TRUE.  Default is
#'   maximum available.  Otherwise, set with an integer. See details for
#'   considerations.
#' @param invertD,invertSd A logical indicating whether or not to invert the D
#'   or S matrix
#' @param returnA,returnS Logical, indicating if the numerator relationship
#'   matrix (A or S) should be stored and returned.
#' @param det Logical, indicating if the determinant of the D or Sd matrix
#'   should be returned.
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{A,S }{the A or S matrix in sparse matrix form}
#'     \item{D,Sd }{the D or Sd matrix in sparse matrix form}
#'     \item{logDet }{the log determinant of the D or Sd matrix} 
#'     \item{Dinv,Sdinv }{the inverse of the D or inverse of the Sd matrix in
#'       sparse matrix form}
#'     \item{listDinv,listSdinv }{the three column form of the non-zero 
#'       elements for the inverse of the D or the inverse of the Sd matrix}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeDsim}}, \code{\link{makeSdsim}}
#' @references Quaas, R.L. 1995. Fx algorithms. An unpublished note.
#' 
#' Lynch M., & Walsh, B. 1998. Genetics and Analysis of Quantitative Traits.
#' Sinauer, Sunderland, Massachusetts.
#' 
#' Shaw, R.G., D.L. Byers, and F.H. Shaw. 1998. Genetic components of variation
#' in Nemophila menziesii undergoing inbreeding: Morphology and flowering time.
#' Genetics. 150:1649-1661.
#' 
#' Wolak, M.E. and L.F. Keller. 2014. Dominance genetic variance and inbreeding
#' in natural populations. In Quantitative Genetics in the Wild, A.
#' Charmantier, L.E.B. Kruuk, and D. Garant eds. Oxford University Press, pp.
#' 104-127.
#' @examples
#' 
#'   DinvMat <- makeD(Mrode9, parallel = FALSE)$Dinv
#' 
#'   SdinvMat <- makeSd(FG90, heterogametic = "0", parallel = FALSE)$Sdinv
#'   # Check to make sure getting correct elements
#'   ## `simPedDFC()` for pedigree with 4 unique sex-linked dominance relatedness values
#'   uSdx <- unique(makeSd(simPedDFC(3), heterogametic = "M", returnS = FALSE)$Sd@x)
#'   stopifnot(all(uSdx %in% c(1, 0.5, 3/16, 1/16))) #<-- must match one of these 4
#' 
#' @export
makeD <- function(pedigree, parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertD = TRUE, returnA = FALSE, det = TRUE){

  numeric.pedigree <- numPed(pedigree) 
  N <- dim(pedigree)[1]
  A <- makeA(pedigree)
  dA <- diag(A)

  if(parallel){
     if(length(A@x)/ncores < 10){
        warning("pedigree too small - 'parallel' set to FALSE instead")
        parallel <- FALSE
     }
  }

  if(!parallel){
     cat("starting to make D...")
     Cout <- .C("dij",
                as.integer(numeric.pedigree[, 2] - 1), 
		as.integer(numeric.pedigree[, 3] - 1), 
		as.integer(A@i), 			
		as.integer(A@p),                        
		as.double(A@x/2),                       
		as.integer(N),                           
		as.double(rep(0, length(A@x))),         
		as.integer(rep(0, length(A@i))),        
                as.integer(rep(0, N)),                  
		as.integer(0))	                        

     D <- Matrix(0, N, N, sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
     D@uplo <- "U"
     D@i <- Cout[[8]][1:Cout[[10]]]
     D@p <- c(Cout[[9]], Cout[[10]])
     D@x <- Cout[[7]][1:Cout[[10]]]
     diag(D) <- 2 - dA

     if(!returnA) A <- NULL
     rm("Cout")

   } else{
        listA <- data.frame(Row = as.integer(rep(1:length(A@p[-1]), diff(A@p))), Column = as.integer(A@i + 1))
        wrap_dij <- function(x){
           sub_lA <- listA[min(x):max(x), 1:2]
           lA_r <- dim(sub_lA)[1]
           Cout <- .C("dijp",
		as.integer(numeric.pedigree[, 2] - 1),
		as.integer(numeric.pedigree[, 3] - 1), 
                as.integer(lA_r),  
		as.integer(sub_lA[, 1] - 1),  
		as.integer(sub_lA[, 2] - 1),  
		as.integer(A@i),  
		as.integer(A@p), 
		as.double(A@x/2),  
		as.double(rep(0, lA_r)))
         Cout[[9]]
        }

        cat("starting to make D...")
        Dijs <- parallel::pvec(seq(1, dim(listA)[1], 1), FUN = wrap_dij, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
  
        D <- Matrix(0, N, N, sparse = TRUE, dimnames = list(as.character(pedigree[, 1]), NULL))
        D@uplo <- "U"
        D@i <- A@i
        D@p <- A@p
        if(!returnA) A <- NULL
        D@x <- Dijs
        D <- drop0(D)
        diag(D) <- 2 - dA

     }

  cat(".done", "\n")
  
  if(det) logDet <- determinant(D, logarithm = TRUE)$modulus[1] else logDet <- NULL
  if(invertD){
    cat("starting to invert D...")
    Dinv <- as(solve(D), "dgCMatrix")
      Dinv@Dimnames <- D@Dimnames
    cat(".done", "\n")
    listDinv <- sm2list(Dinv, rownames=pedigree[,1], colnames=c("row", "column", "Dinverse"))
 return(list(A = A, D = D, logDet = logDet, Dinv=Dinv, listDinv=listDinv))
  } else{
    return(list(A = A, D = D, logDet = logDet))
    } 
}

