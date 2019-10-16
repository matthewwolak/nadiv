#' Functions to conduct gene dropping through a pedigree
#' 
#' Functions that perform and summarize gene dropping conducted on supplied pedigrees
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA', '0'
#' , or '*'.
#'
#' \code{parallel} = TRUE should only be used on Linux or Mac operating systems
#' (i.e., not Windows).
#'
#' Founder allelic values (the alleles assigned to an individual's maternal,
#' paternal, or both haplotypes when the maternal, paternal, or both parents are
#' missing) are equivalent positive and negative integer values corresponding to
#' the maternal and paternal haplotypes, respectively. For example, if the first
#' individual in the pedigree has two unknown parents it will have the following
#' two allelic values: 1=maternal haplotype and -1=paternal haplotype.  
#'
#' @aliases geneDrop
#' @param pedigree A pedigree with columns organized: ID, Dam, Sire.
#' @param N The number of times to iteratively trace alleles through the
#'   pedigree
#' @param parallel A logical indicating whether or not to use parallel
#'   processing. Note, this may only be available for Mac and Linux operating
#'   systems.
#' @param ncores The number of cpus to use when constructing the dominance
#'   relatedness matrix. Default is all available.
#' @param \dots Other arguments that can be supplied to alter what summaries are
#'   reported.
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{IDs }{Original identities in the pedigree}
#'     \item{maternal }{Simulated maternal haplotypes}
#'     \item{paternal }{Simulated paternal haplotypes} 
#'     \item{numericPedigree }{Pedigree in class \code{numPed} for convenient
#'       post-processing of haplotypes}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeDsim}}
#' @examples
#' 
#'   geneDrop(Mrode2, N = 10
#' 
#' @export

geneDrop <- function(pedigree, N, parallel = FALSE, ncores = getOption("mc.cores", 2L), ...){

  nPed <- if(is.numPed(pedigree)) pedigree else numPed(pedigree)
  n <- nrow(pedigree)
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  #FIXME allow alleles to be specified, 
  ## but associate user supplied alleles with integers (work on integers in c++)
  dalleles <- salleles <- vector("integer", length = n) 
  #TODO allow supplied inbreeding coefficients so founders can be inbred
  dalleles[dfounders] <- as.integer(dfounders)
  salleles[sfounders] <- as.integer(-sfounders)
  Ndalleles <- rep(dalleles, each = N)
  Nsalleles <- rep(salleles, each = N)

  #TODO execute in parallel
  Cout <- .C("genedrop",
	as.integer(Ndalleles),
	as.integer(Nsalleles),
	as.integer(N),
	as.integer(n),
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1))

 return(list(IDs = pedigree[, 1],
	maternal = matrix(Cout[[1]], ncol = N, byrow = TRUE),
	paternal = matrix(Cout[[2]], ncol = N, byrow = TRUE),
	numericPedigree = nPed))
}

