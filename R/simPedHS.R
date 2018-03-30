#' Half-sib pedigree construction
#' 
#' Simulates a pedigree for a half-sib mating design (sometimes also called the
#' North Carolina Design 1).
#' 
#' \code{n} must be greater than or equal to 2, because one male and one female
#' offspring are produced from each mating
#' 
#' Some functions/calculations get bogged down if no two dams have the same ID
#' in the entire pedigree (e.g., \code{aov}).  However, other functions must
#' have unique identifiers for every individual.
#' 
#' @usage simPedHS(s, d, n, uniqueDname = TRUE, prefix = NULL)
#' @param s Number of sires
#' @param d Number of dams per sire
#' @param n Number of offspring per mating (must be > or = 2)
#' @param uniqueDname Logical indicating if dams should have unique names
#' within sire families or throughout the entire pedigree
#' @param prefix Optional prefix to every identity
#' @return A \code{data.frame} with columns corresponding to: id, dam, sire,
#' and sex.  Sex is "M" for males and "F" for females.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{simPedDFC}}
#' @examples
#' 
#'   simPedHS(s = 1, d = 3, n = 2)
#' 
#' @export simPedHS
simPedHS <- function(s, d, n, uniqueDname = TRUE, prefix = NULL){
  if(n<2){ stop("must have more than 1 offspring per family (n > or = 2)")
   }
  sires <- paste("s", seq(1, s, 1), sep="")
  if(uniqueDname) {
    dams <- paste("d", seq(1, (s*d), 1), sep="")
    } else{
      dams <- rep(paste("d", seq(1,d, 1), sep=""), s)
      }
  offspring <- paste("o", seq(1, (s*d*n), 1), sep="")

  if(is.null(prefix)){
    ped <- data.frame(id = c(sires, dams, offspring),
	dam = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(dams, each=n)),
	sire = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(sires, each=d*n)),
	sex = c(rep("M", length(sires)), rep("F", length(dams)), rep(c("M","F"), length.out = length(offspring))))
  } else{
       p <- as.character(prefix)
       ped <- data.frame(id = paste0(p, c(sires, dams, offspring)),
	dam = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(paste0(p, dams), each=n)),
	sire = c(rep(NA, length(sires)), rep(NA, length(dams)), rep(paste0(p, sires), each=d*n)),
	sex = c(rep("M", length(sires)), rep("F", length(dams)), rep(c("M","F"), length.out = length(offspring))))     
    }
ped
}

