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
#' @param prefix Optional prefix to add to every identity
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
#' For a given unit of the design (\code{U} total), \code{2*gpn} 0-generation 
#' (grandparental or GP) individuals are created and paired to make \code{gpn}
#' full-sib families. Then the first \code{fws} families are each allocated
#' \code{s} males/sires and \code{s*(fws-1)} females/dams in the 1 (parental or P)
#' generation. The remaining (\code{gpn-fws}) families (only when: 
#' \code{gpn > fws}) are assigned \code{s*fws} females/dams. If
#' \code{fsn > (s*fws)}, the remaining generation 1 (P) individuals in each 
#' full-sib family (\code{fsn - (s*fws)}) are allocated to each family with
#' equal numbers of females and males [this allows for more individuals to be
#' phenotyped in generation 1 (P) than are used to produce generation 2 (F1)].
#' Generation 2 (F1) is then assigned, based on the mating design in Fairbairn
#' and Roff (2006) - essentially each sire [of the \code{s} per full-sib family
#' in generation 1 (P)] is mated to a female from each of the other \code{gpn-1}
#' full-sib families to produce \code{fsn} offspring (with equal numbers of
#' females and males).
#' 
#' @param U An integer number of units or blocks for the design
#' @param gpn Number of grandparent pairs in the generation 0 (GP)
#'   (must be >= 2). Equals the number of full-sib families in generation 1 (P).
#' @param fsn Number of offspring in each full-sib family of generations 1 and 2
#'   (P and F1 - must be an even number >= 4).
#' @param s Number of sires per full-sib family in generation 1 (P - must be >=2)
#' @param fws Number of generation 1 (P) families with sires. Together, with
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
#'   DFC1 <- simPedDFC(U = 1, gpn = 2, fsn = 4, s = 2, fws = 2)
#' 
#' @export
simPedDFC <- function(U, gpn = 4, fsn = 4, s = 2, fws = 2, prefix = NULL)
{

  if(gpn < 2){
    stop("Number of founding grand-parents ('2*gpn') must be greater than or equal to 4")
  }
  if(fsn < 4){
    stop("Full-sib family size ('fsn') must be greater than or equal to 4")
  }
  if(s < 2){
    stop("Number of sires per full-sib family in generation 1 or P ('s') must be greater than or equal to 2")
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
  N <- n*U
  # allocate output data.frame
  ped_out <- data.frame(id = rep(NA, N), dam = rep(NA, N), sire = rep(NA, N),
	sex = rep(NA, N))

  # Calculate starting/ending indices for each group
  gpsi <- seq(from = 1, by = n, length.out = U)
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


  # locate sires in a generic 1/P generation block (created in `design` below)
  sirepos <- cbind(seq(s*fws), rep(seq(fws), each = s))


  ####################################################################
  # Begin FOR LOOP
  # conduct the following on the `Fx`th unit out of `U` total
  for(Fx in seq.int(U)){
    # make a generic 1/P generation block
    design <- matrix(NA, nrow = fsn, ncol = gpn)
    # setup mating females and males (dams/sires) for a unit
    ## 1/P names are: unit | sire/dam/non-mating male/nom-mating female | gp
    ### family number | full-sib letter (order of full-sib within gp family)
    sires <- paste0(rep(paste0(p, "u", Fx, "_s", seq(fws)), each = s),
                                                         pfsletters[seq(pspu)])

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
        } else design[y, ] <- paste0(p, "u", Fx, "_", "f", seq(gpn),
                                                                  pfsletters[y]) 
      }
    }

    # Record male==1==TRUE or female==0==FALSE for a unit of the 1/P generation
    sexDesign <- grepl(paste("^", p, "u", Fx, "_s", sep = ""), design) |
                           grepl(paste("^", p, "u", Fx, "_m", sep = ""), design)

    #########################
    # Fill-in ped_out
    ## Fill-in 0/GP positions for `Fx`th unit of the pedigree
    ped_out[gpsi[Fx]:gpei[Fx], ] <- cbind(c(paste0(p, "u", Fx, "_gs", seq(gpn)),
                                           paste0(p, "u", Fx, "_gd", seq(gpn))),
	rep(NA, 2*gpn),
	rep(NA, 2*gpn),
	rep(c("M", "F"), each = gpn))

    ## Fill-in 1/P positions for `Fx`th unit of the pedigree
    ped_out[psi[Fx]:pei[Fx], ] <- cbind(c(design),
	rep(paste0(p, "u", Fx, "_gd", seq(gpn), sep = ""), each = fsn),
	rep(paste0(p, "u", Fx, "_gs", seq(gpn), sep = ""), each = fsn),
	c("F", "M")[sexDesign+1])

    ## Fill-in 2/F1 positions for `Fx`th unit of the pedigree
    ped_out[f1si[Fx]:f1ei[Fx], ] <- cbind(paste0(p, "u", Fx, "_m",
          rep(seq(pdpu), each = fsn),
          rep(c("m", "f"), each = (fsn/2)),
          rep(seq(fsn/2), (2*pdpu))),
	rep(unlist(sapply(t(design), FUN = function(x) {
	   x[grepl(paste("^", p, "u", Fx, "_d", sep = ""), x)]})), each = fsn),
	rep(sires, each = ((gpn-1)*fsn)),
	rep(c("M", "F"), each = (fsn/2)))
  }  #<-- end for Fx LOOP
  ####################################################################

  ped_out$id <- as.factor(ped_out$id)
  ped_out$dam <- as.factor(ped_out$dam)
  ped_out$sire <- as.factor(ped_out$sire)
  ped_out$sex <- as.factor(ped_out$sex)

 return(ped_out)
}









################################################################################
#' Middle Class Neighborhood pedigree construction
#' 
#' Simulates a pedigree for the \dQuote{middle class neighborhood} mating design
#' (Shabalina, Yampolsky, and Kondrashov 1997).
#' 
#' This creates a pedigree following a breeding design which maintains equal
#' contributions to the next generation by each family in the design. It
#' effectively removes the effect of natural selection which makes it amenable
#' to quantify the contribution of mutations to phenotypic variance over the
#' course of the breeding design.
#'
#' For a starting pedigree template (\code{pedTemp}), the last generation is used
#' as parents to begin the breeding design for the next \code{g} generations.
#' The number of families in the last generation of the template pedigree
#' (\code{pedTemp}) will be the number of families in each generation.
#'
#' Alternatively, if no template pedigree is provided (\code{pedTemp=NULL}),
#' \code{Nfam} number of families will be produced in the first generation from
#' \code{Nfam} unique sire and \code{Nfam} unique dams.
#'
#' Either \code{pedTemp} or \code{Nfam} must be \code{NULL}, but not both.
#' 
#' @param pedTemp A \code{data.frame} pedigree of a template pedigree from which
#'   the middle class neighborhood design should continue. If \code{NULL}, a new
#'   pedigree will be created with \code{Nfam} families.
#' @param g Integer number of generations to produce from the middle class
#'   neighborhood design
#' @param Nfam Integer number of families with which to start a new pedigree
#'   following the middle class neighborhood design. 
#' @param noff Integer number of full-sib offspring produced by each family
#'   (must be >=2).
#'
#' @return A \code{data.frame} with columns corresponding to: id, dam, sire, sex,
#'   and generation. Sex is \code{M} for males and \code{F} for females. The
#'   first generation produced in the middle class neighborhood scheme is assigned
#'   a value of \dQuote{1}, with their parents being assigned to generation
#'   \code{0}. If \code{pedTemp} was provided, the generations from this pedigree
#'   will be denoted with negative integers.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{simPedHS}}, \code{\link{simPedDFC}}
#' @references Shabalina, S.A, L.Y. Yampolsky, and A.S. Kondrashov. 1997. Rapid
#' decline of fitness in panmictic populations of Drosophila melanogaster
#' maintained under relaxed natural selection. Proc. Natl. Acad. Sci. USA.
#' 94:13034-13039.
#' @examples
#'  # No template pedigree provided - start from scrtach
#'   mcn1 <- simPedMCN(pedTemp = NULL, g = 3, Nfam = 4, noff = 2)
#' 
#'  # Provide a template pedigree (half-sib design)
#'   hsped <- simPedHS(s = 2, d = 2, n = 4)
#'   mcnHS <- simPedMCN(pedTemp = hsped, g = 3)
#' @export
simPedMCN <- function(pedTemp, g, Nfam = NULL, noff = 2)
{

  temp <- !is.null(pedTemp)
  if(temp && !is.null(Nfam)) warning("When `pedTemp` provided, `Nfam` is ignored")
  if(!temp && is.null(Nfam)){
    stop("Must provide argument to either `pedTemp` or `Nfam` (but not both)")
  }
  if(noff < 2) stop("'noff' must be >= 2")
  
  if(temp){
    #Find last generation from `pedTemp`
    pedTemp$gen <- genAssign(pedTemp)
    # Create head of a new pedigree for this generation
    pedSt <- pedTemp[which(pedTemp$gen == max(pedTemp$gen)), ]
      names(pedSt)[1:3] <- c("id", "dam", "sire")  #<-- standardize within function
    # Determine number of unique families
    Nfam <- sum(!is.na(unique(paste(as.character(pedSt[, 2]),
      as.character(pedSt[, 3]), sep = "_"))))
    # Figure out if a "sex" column is present in pedTemp
    sexCol <- NULL
    if(ncol(pedTemp) > 3){
      dmatch <- match(pedSt$dam, pedTemp[, 1])
      smatch <- match(pedSt$sire, pedTemp[, 1])
      for(c in 4:ncol(pedTemp)){
        nUniqDcVals <- length(unique(pedTemp[dmatch, c]))
        nUniqScVals <- length(unique(pedTemp[smatch, c]))
        if((nUniqDcVals == 1) & (nUniqScVals == 1)){
          if(unique(pedTemp[dmatch, c]) != unique(pedTemp[smatch, c])){
            sexCol <- c
            break
          }
        }  # end if 1 unique value each
      }  # end for c 

      if(is.null(sexCol)) pedSt$sex <- NA else{
        pedSt <- pedSt[, c(1:3, sexCol)]
          names(pedSt[, 4]) <- "sex"  #<-- standardize within function
      }  
    } else{ 
        pedSt$sex <- NA
        # Arbitrarily assign first ~half Male and remainder Female
        tmpi <- floor(0.5 * nrow(pedSt))
        pedSt$sex[1:tmpi] <- "M"
        pedSt$sex[(1+tmpi):nrow(pedSt)] <- "F"        
      }  
  } else{
      pedSt <- data.frame(id = c(paste0("s", seq(Nfam)), paste0("d", seq(Nfam))),
        dam = NA, sire = NA,
        sex = rep(c("M", "F"), each = Nfam))
    }  #<-- end if/else temp

  ######################################################
  # Figure out total number of new individuals to create
  N <- noff * Nfam * g  
  ped <- data.frame(id = c(as.character(pedSt$id),
      paste(paste0("mcn", rep(seq(Nfam * g), each = noff)),
      rep(seq(noff), Nfam * g), sep = "_")),
    dam = c(as.character(pedSt$dam), rep(NA, N)),
    sire = c(as.character(pedSt$sire), rep(NA, N)),
    sex = c(as.character(pedSt$sex), rep(c("M", "F"), length.out = N)),
    gen = c(rep(0, nrow(pedSt)), rep(seq(g), each = noff * Nfam)))
  
  parentDS <- rep(NA, nrow(ped))
    # find which missing both parents
    bthParNA <- (is.na(ped$dam) & is.na(ped$sire))
    # assign unique parent pair ID to all individuals with atleast 1 parent known
    parentDS[!bthParNA] <- paste0(as.character(ped$dam[!bthParNA]),
      as.character(ped$sire[!bthParNA]))
    # CHECK: if ALL 0 generation (base of `ped`) parents NA
    ## if no pedigree given to function, give phantom parents, else STOP
    if(all(is.na(parentDS[which(ped$gen == 0)]))){
      if(temp){
        stop("If `pedTemp` provided, all parents of the last generation cannot be NA")
      }  
      # `pedSt` created because !temp, sires and dams from 2-offspring families   
      parentDS[which(ped$gen == 0)] <- paste0(rep(paste0("phtmD", seq(Nfam)), 2),
        rep(paste0("phtmS", seq(Nfam)), 2))
    }
  dpool <- spool <- rep(NA, Nfam) 
  iparents <- matrix(NA, nrow = Nfam, ncol = 2)
      
      
  ###### For LOOP through generations #########
  st <- nrow(pedSt) + 1; end <- st - 1 + noff * Nfam  
  for(i in seq(g)){
    # pick individuals to be dams or sires (`dpool` or `spool`) for generation i
    cnt <- 1
    for(u in unique(parentDS[which(ped$gen == (i-1))])){
      uM <- ped$id[which((pDSinU <- parentDS %in% u) & ped$sex == "M")]
        if((luM <- length(uM)) == 0){
          stop(paste("Not every family in generation", i-1,
            "has a male"))                 
        }
      uF <- ped$id[which(pDSinU & ped$sex == "F")]
        if((luF <- length(uF)) == 0){
          stop(paste("Not every family in generation", i-1,
            "has a female"))                 
        }
      # Funny thing with sample below, because sample(1 integer) invokes sample.int  
      spool[cnt] <- uM[sample(x = luM, size = 1)]
      dpool[cnt] <- uF[sample(x = luF, size = 1)]
      cnt <- cnt + 1
    }  #<-- end for u
    
    # create mating pairs
    ## Since went by unique family in previous generation,
    ### brother and sister should reside at same index of dpool and spool
    usedSindex <- c()
    if(Nfam > 2){  
      for(u in 1:(Nfam-2)){  # leave 2 so don't end with sibling as only sire left
        # randomly choose male (avoiding brother of the dam)
        usire <- sample(spool[-c(u, usedSindex)], size = 1)
          usedSindex <- c(usedSindex, which(spool == usire))
        iparents[u, 1:2] <- c(dpool[u], usire)
      }
      uleft <- seq(Nfam)[-usedSindex]
      u2 <- c(Nfam - 1, Nfam)
      if(uleft[1] == u2[1] | uleft[2] == u2[2]) uleft <- rev(uleft)
      iparents[u2, 1:2] <- cbind(dpool[u2], spool[uleft])
    } else iparents[, 1:2] <- cbind(dpool, rev(spool))

    # stick mating pairs into pedigree (repeat for noff per family) for i generation
    ped[st:end, c("dam", "sire")] <- matrix(rep(c(iparents), each = noff), ncol = 2)
      # update `parentDS` too
      parentDS[st:end] <- paste0(as.character(ped$dam[st:end]),
        as.character(ped$sire[st:end]))
        
    # move indices
    st <- end + 1
    end <- st - 1 + noff * Nfam 
  }  #<-- end for i
  
  #############################################
  # stitch graft ped onto remainder of pedTemp
  if(temp){
    pedHd <- pedTemp[which(pedTemp$gen < max(pedTemp$gen)), 1:3]
      names(pedSt)[1:3] <- c("id", "dam", "sire")  #<-- standardize within function
    if(is.null(sexCol)){
      pedHd$sex <- NA
      pedHd$sex[which(pedHd$id %in% pedTemp[, 2])] <- "F"
      pedHd$sex[which(pedHd$id %in% pedTemp[, 3])] <- "M"
    } else{
        pedHd$sex <- pedTemp[which(pedTemp$gen < max(pedTemp$gen)), sexCol]
      }  #<-- end if/else 3 columns in pedTemp
    pedHd$gen <- pedTemp$gen[which(pedTemp$gen < max(pedTemp$gen))] - max(pedTemp$gen)
    # Add head of pedTemp to top of ped
    ped <- rbind(pedHd, ped)
  }  #<-- end if a template pedigree supplied
  
 return(ped)
} 

