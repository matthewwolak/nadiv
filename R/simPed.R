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
#' @param s Number of sires
#' @param d Number of dams per sire
#' @param n Number of offspring per mating (must be > or = 2)
#' @param uniqueDname Logical indicating if dams should have unique names
#'   within sire families or throughout the entire pedigree
#' @param prefix Optional prefix to every identity
#'
#' @return A \code{data.frame} with columns corresponding to: id, dam, sire,
#'   and sex. Sex is "M" for males and "F" for females.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{simPedDFC}}
#' @examples
#' 
#'   simPedHS(s = 1, d = 3, n = 2)
#' 
#' @export
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






################################################################################
#' Double first cousin pedigree construction
#' 
#' Simulates a pedigree for the \dQuote{double first cousin} mating design
#' (Fairbairn and Roff 2006).
#' 
#' This is an adaption to a half-sib breeding design which also produces first
#' cousins and double first cousins. Double first cousins are produced by
#' mating two brothers to two sisters (the offspring of the resulting two
#' families are double first cousins with one another). This is described in
#' Fairbairn and Roff (2006) as being particularly effective for separating
#' autosomal additive genetic variance from sex chromosomal additive genetic
#' variance. It is also amenable to estimating dominance variance, however, it
#' still has difficulty separating dominance variance from common maternal
#' environmental variance (Meyer 2008).
#'
#' For a given unit of the design (\code{F} total), 2*\code{gpn} 0-generation 
#' (grandparental) individuals are created and paired to make \code{gpn} full-
#' sib families. Then the first \code{fws} families are each allocated \code{s}
#' males/sires and \code{s*(fws-1)} females/dams in the 1/P generation. The
#' remaining (\code{gpn-fws}) families (only when: \code{gpn > fws}) are assigned
#' \code{s*fws} females/dams. If \code{fsn > (s*fws)}, the remaining 1/P 
#' generation individuals in each full-sib family (\code{fsn - (s*fws)}) are
#' allocated to each family with equal numbers of females and males (this allows
#' for more individuals to be phenotyped in the 1/P generation than are used to
#' produce the 2/F1 generation). The 2/F1 generation is then assigned, based on
#' the mating design in Fairbairn and Roff (2006) - essentially each sire (of
#' the \code{s} per full-sib family in the 1/P generation) is mated to a female
#' from each of the other \code{gpn-1} full-sib families to produce \code{fsn}
#' offspring (with equal numbers of females and males).
#' 
#' @param F Number of blocks for the design
#' @param gpn Number of grandparent pairs in the 0/GP generation
#'   (must be >= 2). Equals the number of full-sib families in the
#'   1/P generation.
#' @param fsn Number of offspring in each full-sib family of the 1/P
#'   and 2/F1 generations (must be an even number >= 4).
#' @param s Number of sires per full-sib family in the 1/P generation
#'   (must be >=2)
#' @param fws Number of 1/P generation families with sires. Together, with
#'   \code{s}, sets up how cousins and double first cousins are produced 
#' @param prefix Optional prefix to add to every identity
#'
#' @return A \code{data.frame} with columns corresponding to: id, dam, sire,
#'   and sex. Sex is \code{M} for males and \code{F} for females.
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
#' @export
simPedDFC <- function(F, gpn = 4, fsn = 4, s = 2, fws = 2, prefix = NULL)
{

  if(gpn < 2){
    stop("Number of founding grand-parents ('2*gpn') must be greater than or equal to 4")
  }
  if(fsn < 4){
    stop("Full-sib family size ('fsn') must be greater than or equal to 4")
  }
  if(s < 2){
    stop("Number of sires per full-sib family in 1/P generation ('s') must be greater than or equal to 2")
  }
  #FIXME let fsn != even number (then incorporate this below for P and F generations)
  if(floor(fsn/2) != (fsn/2)){
    stop("Full-sib family size ('fsn') must be an even number")
  }
  if(fsn < (s*fws)){
    stop("Full-sib family size ('fsn') must be greater than or equal to the product: number of families with sires ('fws') * number of sires per family ('s')")
  }

  if(is.null(prefix)) p <- "" else p <- as.character(prefix)

  # Calculate numbers of types of individuals *PER UNIT*
  ## 0/GP generation (males and females)
  gppu <- 2*gpn
  ## 1/P sires (can be << 0.5*fsn if fsn large - these are 1/P males that do mate)
  pspu <- s*fws
  ## 1/P dams (can be << 0.5*fsn if fsn large - these are 1/P females that do mate)
  pdpu <- (pspu*gpn) - pspu
  ## 1/P individuals that don't contribute to F1 (alternate male/female)
  pnmpu <- ifelse(((fsn*gpn) - (pspu + pdpu)) <= 0, 0, (fsn*gpn)-(pspu+pdpu))
  ## 2/F1 individuals per unit
  f1pu <- pdpu*fsn

  # total individuals in a unit
  n <- gppu + pspu + pdpu + pnmpu + f1pu
  # total individuals in pedigree
  N <- n*F
  # allocate output data.frame
  ped_out <- data.frame(id = rep(NA, N), dam = rep(NA, N), sire = rep(NA, N),
	sex = rep(NA, N))

  # Calculate starting/ending indices for each group
  gpsi <- seq(from = 1, by = n, length.out = F)
    gpei <- gpsi + gppu - 1
  psi <- gpei + 1
    pei <- psi + pspu + pdpu + pnmpu - 1
  f1si <- pei + 1
    f1ei <- f1si + f1pu - 1

  # designate the letters for distinguishing among full-sibs in same family of
  ## the 1/P generation
  if(fsn > 26){
    pfsletters <- unlist(Reduce(paste0, 
       replicate(fsn %/% length(letters), letters, simplify = FALSE),
       init = letters,
       accumulate = TRUE))[1:fsn]
  } else pfsletters <- letters[seq(fsn)]


  # make a generic 1/P generation block
  design <- matrix(NA, nrow = fsn, ncol = gpn)
  sirepos <- cbind(seq(s*fws), rep(seq(fws), each = s))

  # conduct the following on the `Fx`th unit out of `F` total
  unitFun <- function(Fx){
    # setup mating females and males (dams/sires) for a unit
    ## 1/P names are: unit | sire/dam/non-mating male/nom-mating female | gp family number | full-sib letter (order of full-sib within gp family)
    sires <- paste0(rep(paste0(p, "u", Fx, "_s", seq(fws)), each = s), pfsletters[seq(pspu)])

    for(x in 1:pspu){
      design[x, sirepos[x,2]] <- sires[x]
      damNumb <- which(is.na(design[x, ]))
      design[x, damNumb] <- paste0(p, "u", Fx, "_d", damNumb, pfsletters[x])
    }
    if(pnmpu != 0){
      for(y in seq(pspu + 1, fsn, by = 1)){
        # below so that alternate sexes each row
        if((y/2) != ceiling(y/2)){
          design[y, ] <- paste0(p, "u", Fx, "_", "m", seq(gpn), pfsletters[y])  
        } else design[y, ] <- paste0(p, "u", Fx, "_", "f", seq(gpn), pfsletters[y]) 
      }
    }

    # Record male==1==TRUE or female==0==FALSE for a unit of the 1/P generation
    sexDesign <- grepl(paste("^", p, "u", Fx, "_s", sep = ""), design) | grepl(paste("^", p, "u", Fx, "_m", sep = ""), design)

    #########################
    # Fill-in ped_out
    ## Fill-in 0/GP positions for `Fx`th unit of the pedigree
    ped_out[gpsi[Fx]:gpei[Fx], ] <<- cbind(c(paste0(p, "u", Fx, "_gs", seq(gpn)), paste0(p, "u", Fx, "_gd", seq(gpn))),
	rep(NA, 2*gpn),
	rep(NA, 2*gpn),
	rep(c("M", "F"), each = gpn))

    ## Fill-in 1/P positions for `Fx`th unit of the pedigree
    ped_out[psi[Fx]:pei[Fx], ] <<- cbind(c(design),
	rep(paste0(p, "u", Fx, "_gd", seq(gpn), sep = ""), each = fsn),
	rep(paste0(p, "u", Fx, "_gs", seq(gpn), sep = ""), each = fsn),
	c("F", "M")[sexDesign+1])

    ## Fill-in 2/F1 positions for `Fx`th unit of the pedigree
    ped_out[f1si[Fx]:f1ei[Fx], ] <<- cbind(paste0(p, "u", Fx, "_m", rep(seq(pdpu), each = fsn), rep(c("m", "f"), each = (fsn/2)), rep(seq(fsn/2), (2*pdpu))),
	rep(unlist(sapply(t(design), FUN = function(x) {x[grepl(paste("^", p, "u", Fx, "_d", sep = ""), x)]})), each = fsn),
	rep(sires, each = ((gpn-1)*fsn)),
	rep(c("M", "F"), each = (fsn/2)))


   invisible(Fx)
   }


  sapply(seq(F), FUN = unitFun)

 return(ped_out)
}

