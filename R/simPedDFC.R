#' Double first cousin pedigree construction
#' 
#' Simulates a pedigree for the \dQuote{double first cousin} mating design
#' (Fairbairn and Roff 2006).
#' 
#' This is an adaption to a half-sib breeding design which also produces first
#' cousins and double first cousins.  Double first cousins are produced by
#' mating two brothers to two sisters (the offspring of the resulting two
#' families are double first cousins with one another).  This is described in
#' Fairbairn and Roff (2006) as being particularly effective for separating
#' autosomal additive genetic variance from sex chromosomal additive genetic
#' variance.  It is also amenable to estimating dominance variance, however, it
#' still has difficulty separating dominance variance from common maternal
#' environmental variance (Meyer 2008).
#' 
#' @usage simPedDFC(F, gpn = 4, fsn = 4, s = 2, prefix = NULL)
#' @param F Number of blocks for the design
#' @param gpn Number of grandparents in the first/GP generation (must be >= 2)
#' @param fsn Number of offspring in the full-sib families of the second/P
#' generation (must be an even number >= 4)
#' @param s Number of sires per full-sib family in the second/P generation
#' (must be >=2)
#' @param prefix Optional prefix to add to every identity
#' @return A \code{data.frame} with columns corresponding to: id, dam, sire,
#' and sex.  Sex is \code{M} for males and \code{F} for females.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{simPedHS}}, \code{\link{warcolak}}
#' @references Fairbairn, D.J. and D.A. Roff. 2006.  The quantitative genetics
#' of sexual dimorphism: assessing the importance of sex-linkage. Heredity
#' 97:319-328.
#' 
#' Meyer, K. 2008.  Likelihood calculations to evaluate experimental designs to
#' estimate genetic variances. Heredity 101:212-221.
#' @examples
#' 
#'   DFC1 <- simPedDFC(F = 1, gpn = 2, fsn = 4, s = 2)
#' 
#' @export simPedDFC
simPedDFC <- function(F, gpn = 4, fsn = 4, s = 2, prefix = NULL)
{

if(gpn < 2) stop("Number of founding grand-parents ('gpn') must be greater than or equal to 2")
if(fsn < 4) stop("Full-sib family size ('fsn') must be greater than or equal to 4")
if(s < 2) stop("Number of sires per full-sib family in P generation ('s') must be greater than or equal to 2")
if(floor(fsn/2) != (fsn/2)) stop("Full-sib family size ('fsn') must be an even number")

unitFun <- function(Fx){
   design <- matrix(NA, nrow = fsn, ncol = gpn)
   rc <- cbind(c(1:(floor(fsn/s)*s)), rep(1:floor(fsn/s), each = s))
   sires <- paste(rep(paste("u", Fx, "_s", seq.int(fsn/s), sep = ""), each = s), letters[1:(floor(fsn/2)*s)], sep = "")

   for(x in 1:dim(rc)[1]){
      design[x, rc[x,2]] <- sires[x]
      damNumb <- which(is.na(design[x, ]))
      design[x, damNumb] <- paste("u", Fx, "_d", damNumb, letters[x], sep = "")
   }
   sexDesign <- grepl(paste("^u", Fx, "_s", sep = ""), design)

   if(is.null(prefix)){
     unitPed <- data.frame(id = as.character(c(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), paste("u", Fx, "_gd", seq.int(gpn), sep = ""), c(design), paste("u", Fx, "_m", rep(seq.int((fsn*(gpn-1))), each = fsn), rep(c("m", "f"), each = (fsn/2)), rep(seq.int(fsn/2), (fsn*(gpn-1))), sep = ""))),
	dam = as.character(c(rep(NA, (2*gpn)), rep(paste("u", Fx, "_gd", seq.int(gpn), sep = ""), each = fsn), rep(unlist(sapply(t(design), FUN = function(x) {x[grepl(paste("^u", Fx, "_d", sep = ""), x)]})), each = fsn))),
	sire = as.character(c(rep(NA, (2*gpn)), rep(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), each = fsn), rep(sires, each = ((gpn-1)*fsn)))),
	sex = c(rep("M", gpn), rep("F", gpn), sapply(sexDesign, FUN = function(x) {if(x) "M" else "F"}), rep(rep(c("M", "F"), each = (fsn/2)), ((gpn-1)*fsn))))
   } else{
       p <- as.character(prefix)   
       unitPed <- data.frame(id = as.character(paste0(p, c(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), paste("u", Fx, "_gd", seq.int(gpn), sep = ""), c(design), paste("u", Fx, "_m", rep(seq.int((fsn*(gpn-1))), each = fsn), rep(c("m", "f"), each = (fsn/2)), rep(seq.int(fsn/2), (fsn*(gpn-1))), sep = "")))),
	dam = as.character(c(rep(NA, (2*gpn)), paste0(p, c(rep(paste("u", Fx, "_gd", seq.int(gpn), sep = ""), each = fsn), rep(unlist(sapply(t(design), FUN = function(x) {x[grepl(paste("^u", Fx, "_d", sep = ""), x)]})), each = fsn))))),
	sire = as.character(c(rep(NA, (2*gpn)), paste0(p, c(rep(paste("u", Fx, "_gs", seq.int(gpn), sep = ""), each = fsn), rep(sires, each = ((gpn-1)*fsn)))))),
	sex = c(rep("M", gpn), rep("F", gpn), sapply(sexDesign, FUN = function(x) {if(x) "M" else "F"}), rep(rep(c("M", "F"), each = (fsn/2)), ((gpn-1)*fsn))))
     }


   unitPed
   }


 ped_out <- do.call(rbind, lapply(seq.int(F), FUN = unitFun))


return(ped_out)
}

