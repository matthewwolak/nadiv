#' Creates the additive by dominance and dominance by dominance epistatic
#' genetic relationship matrices
#' 
#' Given a pedigree, the matrix of additive by dominance (AD) genetic
#' relatedness, dominance by dominance (DD) genetic relatedness, or both are
#' returned.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' Because of the computational demands of constructing the D matrix (see
#' \code{\link{makeD}}), this function allows for the inverses that are derived
#' from the D matrix (i.e., D-inverse, AD-inverse, and DD-inverse)to be
#' constructed at the same time.  This way, the D matrix will only have to be
#' constructed once for use in the three separate genetic relatedness inverse
#' matrices that depend upon it.  However, using the \code{output} and
#' \code{invertD} options in different combinations will ensure that only the
#' desired matrix inverses are constructed.
#' 
#' \code{parallel} = TRUE should only be used on Linux or Mac OSes (i.e., not
#' Windows).
#' 
#' Both the AD and DD matrix are computed from the Hadamard product of the
#' respective matrices (see also, \code{\link{makeAA}}).
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param output Character(s) denoting which matrix and its inverse is to be
#'   constructed.
#' @param parallel A logical indicating whether or not to use parallel
#'   processing. Note, this may only be available on Mac and Linux operating
#'   systems.
#' @param invertD A logical indicating whether or not to invert the D matrix
#' @param det A logical indicating whether or not to return the determinants
#'   for the epistatic relationship matrices
#'
#' @return All of the following will be returned. However, the values of the
#'   \code{output} and \code{invertD} options passed to the function will
#'   determine which of the following are not NULL objects within the list:
#'   \describe{ 
#'     \item{D }{the D matrix in sparse matrix form}
#'     \item{logDetD }{the log determinant of the D matrix}
#'     \item{AD }{the AD matrix in sparse matrix form}
#'     \item{logDetAD }{the log determinant of the AD matrix}
#'     \item{DD }{the DD matrix in sparse matrix form} 
#'     \item{logDetDD }{the log determinant of the DD matrix}
#'     \item{Dinv }{the inverse of the D matrix in sparse matrix form}
#'     \item{ADinv }{the inverse of the AD matrix in sparse matrix form}
#'     \item{DDinv }{the inverse of the DD matrix in sparse matrix form}
#'     \item{listDinv }{the three column form of the non-zero elements for the
#'       inverse of the D matrix}
#'     \item{listADinv }{the three column form of the non-zero elements for the 
#'       inverse of the AD matrix}
#'     \item{listDDinv }{the three column form of the non-zero elements for the 
#'       inverse of the DD matrix}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeA}}, \code{\link{makeD}}, \code{\link{makeAA}}
#' @examples
#' 
#'   Boutput <- makeDomEpi(Mrode9, output = "b", parallel = FALSE, invertD = FALSE)
#'   str(Boutput)
#' 	
#'   DADoutput <- makeDomEpi(Mrode9, output = "AD", parallel = FALSE, invertD = TRUE)
#'   str(DADoutput)
#' 
#' @export
makeDomEpi <- function(pedigree, output = c("AD", "DD", "both"),
	parallel = FALSE, invertD = FALSE, det = TRUE)
{
  type <- match.arg(output)
  Dout <- makeD(pedigree, parallel = parallel, invertD = invertD, returnA = TRUE)

  if(type == "AD"){
    AD <- Dout$A * Dout$D
    if(det) logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1] else logDetAD <- NULL
    ADinv <- as(solve(AD), "dgCMatrix")
      ADinv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    DD <- NULL
    logDetDD <- NULL
    DDinv <- NULL
    listDDinv <- NULL
    }      

  if(type == "DD"){
    DD <- Dout$D * Dout$D
    if(det) logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1] else logDetDD <- NULL
    DDinv <- as(solve(DD), "dgCMatrix")
      DDinv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    listDDinv<-sm2list(DDinv, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
    AD <- NULL
    logDetAD <- NULL
    ADinv <- NULL
    listADinv <- NULL
   }

  if(type == "both"){
    AD <- Dout$A * Dout$D
    ADinv <- Matrix(solve(AD), sparse=TRUE)
    listADinv <-sm2list(ADinv, rownames=pedigree[,1], colnames=c("row", "column", "ADinverse"))
    ADinv <- as(ADinv, "dgCMatrix")
      ADinv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    DD <- Dout$D * Dout$D
    if(det){
      logDetAD <- determinant(AD, logarithm = TRUE)$modulus[1]
      logDetDD <- determinant(DD, logarithm = TRUE)$modulus[1]
    } else{ logDetAD <- logDetDD <- NULL}
    DDinv <- as(solve(DD), "dgCMatrix")
      DDinv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    listDDinv<-sm2list(DDinv, rownames=pedigree[,1], colnames=c("row", "column", "DDinverse"))
     }
return(list(D=Dout$D, logDetD = Dout$logDet, AD=AD, logDetAD = logDetAD, DD=DD, logDetDD = logDetDD, Dinv=Dout$Dinv, ADinv=ADinv, DDinv=DDinv, listDinv=Dout$listDinv, listADinv=listADinv, listDDinv=listDDinv))

}

