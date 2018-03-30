# Generic


#' Generation assignment
#' 
#' Given a pedigree, the function assigns the generation number to which each
#' individual belongs.
#' 
#' 0 is the base population.
#' 
#' Migrants, or any individuals where both parens are unknown, are assigned to
#' generation zero.  If parents of an individual are from two different
#' generations (e.g., dam = 0 and sire = 1), the individual is assigned to the
#' generation following the greater of the two parents (e.g., 2 in this
#' example).
#' 
#' @aliases genAssign genAssign.default genAssign.numPed
#' @usage genAssign(pedigree, \dots{}) \method{genAssigndefault}(pedigree,
#' \dots{}) \method{genAssignnumPed}(pedigree, \dots{})
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param \dots Arguments to be passed to methods
#' @return A vector of values is returned.  This vector is in the same order as
#' the ID column of the pedigree.
#' @author \email{matthewwolak@@gmail.com}
#' @export genAssign
genAssign <- function(pedigree, ...){
  UseMethod("genAssign", pedigree)
}

###############################################################################
# Methods:

genAssign.default <- function(pedigree, ...)
{ 
   n <- nrow(pedigree)
   if(any(apply(pedigree[, 1:3], MARGIN = 2, FUN = function(x){min(x, na.rm = TRUE) < 0}))){
      warning("Negative values in pedigree interpreted as missing values")
      pedigree[pedigree < 0] <- -998
   }
   if(!all(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.numeric)) | any(apply(pedigree[, 1:3], MARGIN = 2, FUN = is.na))){
      pedigree[, 1:3] <- numPed(pedigree[, 1:3])
   }

   Cout <- .C("ga",
	as.integer(pedigree[, 2] - 1),
	as.integer(pedigree[, 3] - 1),
        vector("integer", length = n),
	as.integer(n))
  Cout[[3]]
}

######################################

genAssign.numPed <- function(pedigree, ...)
{ 
   n <- nrow(pedigree)
   Cout <- .C("ga",
	as.integer(pedigree[, 2] - 1),
	as.integer(pedigree[, 3] - 1),
        vector("integer", length = n),
	as.integer(n))
  Cout[[3]]
}

