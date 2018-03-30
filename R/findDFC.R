#' Finds the double first cousins in a pedigree
#' 
#' Given a pedigree, all pairs of individuals that are double first cousins are
#' returned.
#' 
#' When exact = TRUE, only those individuals whose grandparents are completely
#' unrelated will be identified as double first cousins.  When exact = FALSE,
#' as long as the parents of individuals i and j are two sets of siblings
#' (i.e., either sires full brothers/dams full sisters or two pairs of opposite
#' sex full sibs) then i and j will be considered double first cousins.  In the
#' event where the grandparents of i and j are also related, exact = FALSE will
#' still consider i and j full sibs, even though genetically they will be more
#' related than exact = TRUE double first cousins.
#' 
#' \code{parallel} = TRUE should only be used on Linux or Mac OSes (i.e., not
#' Windows).
#' 
#' @param pedigree A pedigree with columns organized: ID, Dam, Sire
#' @param exact A logical statement indicating if individuals who are exactly
#'   double first cousins are to be identified
#' @param parallel A logical statement indicating if parallelization should be
#'   attempted.  Note, only reliable for Mac and Linux operating systems.
#' @param ncores Number of cpus to use, default is maximum available
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{PedPositionList }{gives the list of row numbers for all the
#'       pairs of indidivuals that are related as double first cousins.}
#'     \item{DFC }{gives the list of IDs, as characters, for all the pairs of 
#'       individuals that are related as double first cousins.}
#'     \item{FamilyCnt }{If two individuals, i and j, are double first cousins,
#'       then i's siblings will also be double first cousins with j's siblings.
#'       Therefore, this is the total number of family pairs where offspring
#'       are related as double first cousins.}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @export
findDFC <- function(pedigree, exact = FALSE, parallel = FALSE, ncores = getOption("mc.cores", 2L))
{
  numeric.pedigree <- numPed(pedigree)
  ped <- cbind(numeric.pedigree, genAssign(numeric.pedigree), rep(0, dim(numeric.pedigree)[1]))
  num.out <- ped[ped[,4] >= 2, ]
  ni <- dim(num.out)[1]
  maxid <- max(num.out[,1]) 

  i <- unlist(mapply(rep, num.out[-ni, 1], each = seq((ni-1), 1)))
  j <- unlist(lapply(seq(2,ni), FUN = function(x) num.out[x:ni, 1]))
  if(exact) exct <- 1 else exct <- 0

  if(parallel) {
     wrap_DFC <- function(x){
         i.tmp <- i[min(x):max(x)]
         j.tmp <- j[min(x):max(x)]
         Cout <- .C("dfc",
	    as.integer(numeric.pedigree[, 2] - 1),
            as.integer(numeric.pedigree[, 3] - 1),
	    as.integer(i.tmp - 1),
	    as.integer(j.tmp - 1),
	    as.integer(length(i.tmp)),
	    as.integer(exct))
        Cout[[3]]
     }
     dfcs.vec <- parallel::pvec(seq.int(length(i)), FUN = wrap_DFC, mc.set.seed = FALSE, mc.silent = TRUE, mc.cores = ncores, mc.cleanup = TRUE)
     } else{ 
          Cout <- .C("dfc",
	     as.integer(numeric.pedigree[, 2] - 1),
             as.integer(numeric.pedigree[, 3] - 1),
	     as.integer(i - 1),
	     as.integer(j - 1),
	     as.integer(length(i)),
	     as.integer(exct))
          dfcs.vec <- Cout[[3]]
       }


  yes.dfcs <- which(dfcs.vec == 1)

return(list(PedPositionList = data.frame(i = i[yes.dfcs], j = j[yes.dfcs]), DFC = data.frame(i = pedigree[i[yes.dfcs], 1], j = pedigree[j[yes.dfcs], 1]), FamilyCnt = dim(unique(cbind(pedigree[i[yes.dfcs], 2:3], pedigree[j[yes.dfcs], 2:3])))[1]))
}

