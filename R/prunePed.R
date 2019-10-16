# Generic


#' Prunes a pedigree based on individuals with phenotypes
#' 
#' This function removes individuals who are either not themselves or not
#' ancestors to phenotyped individuals
#' 
#' Often mixed effect models run much faster when extraneous information is
#' removed before running the model. This is particularly so when reducing the
#' number of random effects associated with a relationship matrix constructed
#' from a pedigree.
#' 
#' NOTE: more columns than just a pedigree can be passed in the \code{pedigree}
#' argument.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' This function is very similar to (and the code is heavily borrowed from) a
#' function of the same name in the \code{MCMCglmm} package by Jarrod Hadfield.
#' 
#' @aliases prunePed prunePed.default prunePed.numPed
#' @param pedigree An object, where the first 3 columns correspond to: ID, Dam,
#'   & Sire. See details.
#' @param phenotyped A vector indicating which individuals in the pedigree have
#'   phenotypic information available.
#' @param \dots Arguments to be passed to methods
#'
#' @return The pedigree object (can have more columns than just ID, Dam, and
#'   Sire), where the ID column contains an ID for all individuals who are
#'   actually phenotyped or are an ancestor to an individual with a phenotype
#'   (and are thus informative for estimating parameters in the base 
#'   population).
#' @seealso \code{\link[nadiv]{prepPed}}
#' @examples
#' 
#' 
#' # Make a pedigree (with sex) from the warcolak dataset
#'   warcolak_ped <- warcolak[, 1:4]
#' 
#' # Reduce the number of individuals that have a phenotype for "trait1" in
#'   #the warcolak dataset
#'   t1phenotyped <- warcolak
#'   t1phenotyped[sample(seq.int(nrow(warcolak)), 1500, replace = FALSE), "trait1"] <- NA
#'   t1phenotyped <- t1phenotyped[which(!is.na(t1phenotyped$trait1)), ]
#' 
#' # The following will give a pedigree with only individuals that have a 
#' # phenotype for "trait1" OR are an ancestor to a phenotyped individual.
#'   pruned_warcolak_ped <- prunePed(warcolak_ped, phenotyped = t1phenotyped$ID)
#' 
#' # Now compare the sizes (note, pruned_warcolak_ped retained its column indicating sex.
#' # We could have kept all of the data associated with individuals who had phenotypic
#' # information on "trait1" by instead specifying 
#' #  prunePed(warcolak, phenotyped = t1phenotyped$ID) 
#'   dim(warcolak_ped)
#'   dim(pruned_warcolak_ped)
#' 
#' 
#' @export
prunePed <- function(pedigree, phenotyped, ...){
  UseMethod("prunePed", pedigree)
}

###############################################################################
###############################################################################
#     Methods:
###############################################################################
###############################################################################

################################################
#Borrowed heavily from the 'prunePed' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################

#' @method prunePed default
#' @rdname prunePed
#' @export
prunePed.default <- function (pedigree, phenotyped, ...) {
   nPed <- numPed(pedigree[, 1:3])
   ikeep <- match(phenotyped, pedigree[, 1])
   nind <- length(ikeep) + 1
   while(length(ikeep) != nind){
      nind <- length(ikeep)
      ikeep <- union(c(nPed[ikeep, 2:3]), ikeep)
      ikeep <- ikeep[which(ikeep > 0)]
   }

   pedigree <- pedigree[sort(ikeep), ]
   pedigree[, 1] <- as.factor(as.character(pedigree[, 1]))
   pedigree[, 2] <- as.factor(as.character(pedigree[, 2]))
   pedigree[, 3] <- as.factor(as.character(pedigree[, 3]))
 pedigree
}





#' @method prunePed numPed
#' @rdname prunePed
#' @export
prunePed.numPed <- function (pedigree, phenotyped, ...) {
   ikeep <- match(phenotyped, pedigree[, 1])
   nind <- length(ikeep) + 1
   while(length(ikeep) != nind){
      nind <- length(ikeep)
      ikeep <- union(c(pedigree[ikeep, 2:3]), ikeep)
      ikeep <- ikeep[which(ikeep > 0)]
   }
 pedigree[sort(ikeep), ]
}

