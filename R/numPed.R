################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################



#' Integer Format Pedigree
#' 
#' Conversion, checking, and row re-ordering of a pedigree in integer form of
#' class \sQuote{numPed}.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', '-998', or '*'.
#' 
#' Individuals must appear in the ID column in rows preceeding where they
#' appear in either the Dam or Sire column. See the
#' \code{\link[nadiv]{prepPed}} function if this is not the case.
#' 
#' If pedigree inherits the class "numPed" (from a previous call to
#' \code{numPed()}) and \code{check = TRUE}, the checks are skipped. If
#' \code{check = FALSE} any pedigree will be transformed into a pedigree
#' consisting of integers and missing values denoted by '-998'.
#' 
#' Based on code from the \code{MCMCglmm} package
#' 
#' @aliases numPed is.numPed ronPed
#' @param pedigree A three column pedigree object, where the columns correspond 
#'   to: ID, Dam, & Sire
#' @param check A logical argument indicating if checks on the validity of the 
#'   pedigree structure should be made, but see Details
#' @param x A pedigree of class \sQuote{\code{numPed}}
#' @param i,\dots Index specifying elements to extract or replace: see
#'   \code{\link[base]{[}}
#'
#' @return An S3 object of class \dQuote{numPed} representing the pedigree, 
#'   where individuals are now numbered from 1 to \code{n} and unknown parents 
#'   are assigned a value of \sQuote{-998}.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link[nadiv]{prepPed}}, \code{\link[MCMCglmm]{MCMCglmm}},
#' \code{\link[base]{[}}
#' @examples
#' 
#' (nPed <- numPed(Mrode2))
#' is(nPed)
#' 
#' # re-order and retain class 'numPed'
#' ronPed(nPed, order(nPed[, 2], nPed[, 3]))
#' is(nPed)
#' 
#' @export
numPed <- function(pedigree, check = TRUE){
 if(!is.numPed(pedigree) && check){      
  if(any(pedigree[, 2] == 0, na.rm = TRUE)){
    pedigree[which(pedigree[, 2] == 0), 2] <- NA
    warning("Zero in the dam column interpreted as a missing parent")
  }
  if(any(pedigree[, 3] == 0, na.rm = TRUE)){
    pedigree[which(pedigree[, 3] == 0), 3] <- NA
    warning("Zero in the sire column interpreted as a missing parent")
  }
  if(any(pedigree[, 2] == -998, na.rm = TRUE)){
    pedigree[which(pedigree[, 2] == -998), 2] <- NA
    if(!is.numPed(pedigree)) warning("-998 in the dam column interpreted as a missing parent")
  }
  if(any(pedigree[, 3] == -998, na.rm = TRUE)){
    pedigree[which(pedigree[, 3] == -998), 3] <- NA
    if(!is.numPed(pedigree)) warning("-998 in the sire column interpreted as a missing parent")
  }
  if(any(pedigree[,2] == "*", na.rm = TRUE)) pedigree[which(pedigree[, 2] == "*"), 2] <- NA
  if(any(pedigree[,3] == "*", na.rm = TRUE)) pedigree[which(pedigree[, 3] == "*"), 3] <- NA

  if(all(is.na(pedigree[, 2])) & all(is.na(pedigree[, 3]))){
     stop("All dams and sires are missing")
  }
  if(dim(pedigree)[2] != 3){
     stop("pedigree must have three columns: ID, Dam and Sire")
  }
  if(sum((na.omit(pedigree[, 2]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 2]) == FALSE)){
     stop("individuals appearing as dams but not in pedigree: first use the 'prepPed' function")
  }
  if(sum((na.omit(pedigree[, 3]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 3]) == FALSE)){
     stop("individuals appearing as sires but not in pedigree: first use the 'prepPed' function")
  }
  if(any(duplicated(pedigree[, 1]))){
     stop("some individuals appear more than once in the pedigree")
  }
 }
  nPed <- matrix(as.integer(-998), dim(pedigree)[1], dim(pedigree)[2])
  nPed[, 1] <- as.integer(seq(1, dim(pedigree)[1], 1))
  nPed[, 2] <- match(pedigree[, 2], pedigree[, 1], nomatch = -998)
  nPed[, 3] <- match(pedigree[, 3], pedigree[, 1], nomatch = -998)
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
 if(check){
  if(length(intersect(nPed[, 2][dnmiss], nPed[, 3][snmiss])) > 0 & (length(dnmiss) > 0) & (length(snmiss) > 0)){
      warning("Dams appearing as Sires - assumed selfing in pedigree")
  }
  if(any(nPed[, 2][dnmiss] > nPed[, 1][dnmiss]) & (length(dnmiss) > 0)){
     stop("Offspring appearing before their dams: first use the 'prepPed' function")
  }
  if(any(nPed[, 3][snmiss] > nPed[, 1][snmiss]) & (length(snmiss) > 0)){
     stop("Offspring appearing before their Sires: first use the 'prepPed' function")
  }
  if(any((nPed[, 1] - nPed[, 2]) == 0)){
     cat("Individual(s):", nPed[which((nPed[, 1] - nPed[, 2]) == 0), 1], "\n")
     stop("Individual appearing as its own Dam")
  }
  if(any((nPed[, 1] - nPed[, 3]) == 0)){
     cat("Individual(s):", nPed[which((nPed[, 1] - nPed[, 3]) == 0), 1], "\n")
     stop("Individual appearing as its own Sire")
  }
 }
  nPed <- structure(nPed, class = "numPed")
 nPed
}



#' @ method is numPed
#' @rdname numPed
#' @export
is.numPed <- function(x) inherits(x, "numPed")






# re-ordering rows of object with class 'numPed'

#' @rdname numPed
#' @export
ronPed <- function(x, i, ...){
   r <- structure(unclass(x)[i, ,...], class = "numPed")
   r
}


