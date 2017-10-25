makeSdsim <- function(pedigree, heterogametic, N,
	DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"),
	parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertSd = TRUE, calcSE = FALSE, returnS = FALSE){

  if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")

  dc.model <- match.arg(DosageComp)
  if(is.null(dc.model)){
    warning("Assuming 'ngdc' dosage compensation model")
    dc.model <- "ngdc"
  }
  if(dc.model == "hopi" | dc.model == "hori"){
    cat("Assume sex chromosomal dominance allelic interactions do not occur under 'hopi' or 'hori'\n")
    return(NULL)
  }

  n <- dim(pedigree)[1L]
  nPed <- numPed(pedigree[, 1:3])
  damsex <- pedigree[unique(nPed[, 2])[-1], 4]
  if(any(damsex == heterogametic)){
    pedname <- names(pedigree)
    pedigree <- pedigree[, c(1,3,2,4)]
    names(pedigree) <- pedname
    nPed <- numPed(pedigree[, 1:3])
   cat("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system\n")
  } else cat("Assuming male heterogametic (e.g., XX/XY) sex chromosome system\n")
  sex <- rep(-998, n)
  sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
  sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
  nhom <- sum(sex)  # Number of individuals with homogametic sex chromosomes

  #TODO delete next note once consolidated 'gene dropping' functions/code
  # diverges from `makeDsim()` and follows simplifications in `geneDrop()`
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998 & sex == 1)
  dalleles <- salleles <- vector("integer", length = n)

  dalleles[dfounders] <- as.integer(dfounders)
  salleles[sfounders] <- as.integer(-sfounders)
  Ndalleles <- rep(dalleles, each = N)
  Nsalleles <- rep(salleles, each = N)
  
  cat("making Sdsim ...")

  # diversion to calculate maximum expected entries in sex-chromosome D matrix
  ## based on calculation for sex-chromosome S matrix (additive)
  dnmiss <- which(nPed[, 2] != -998 & sex == 1)
  snmiss <- which(nPed[, 3] != -998 & sex == 1)
  bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998 & sex == 1)
  nSd <- nhom + 2 * length(dnmiss) + 2 * length(snmiss)
  nSd <- nSd + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)

  Cout <- .C("sdsim",
	as.integer(Ndalleles),        # [[1]] N dam alleles (or homogametic sex if ZZ/ZW)
	as.integer(Nsalleles),	      # [[2]] N sire alleles (or heterogametic sex if ZZ/ZW)
	as.integer(N),  	      # [[3]] N (number of replications)
	as.integer(n), 		      # [[4]] n pedigree size
	as.integer(nPed[, 2] - 1),    # [[5]] dam number IDs
	as.integer(nPed[, 3] - 1),    # [[6]] sire number IDs
	as.integer(sex),	      # [[7]] sex or number of homogametic sex chromosomes
	as.integer(rep(0, nSd)),      # [[8]] i slot of sex-chrom. dom. relatedness matrix
	as.integer(rep(0, n+1)),      # [[9]] p slot of matrix
	as.integer(rep(0, nSd)))      # [[10]] x slot of matrix

  nSd <- Cout[[9]][nhom+1]  # change to reflect actual number of non-zeroes
  Sdsim <- Matrix(0, nhom, nhom, sparse = TRUE, dimnames = list(as.character(pedigree[homs, 1]), NULL))
  Sdsim@uplo <- "U"
  Sdsim@i <- Cout[[8]][1:nSd]
  Sdsim@p <- Cout[[9]][1:(nhom+1)]
  Sdsim@x <- Cout[[10]][1:nSd] / N
#TODO in `makeDsim()` stick in the diagonal the values from the analytical D?
## If do, don't need to calculate the diagonals in c++ code ('sdsim.cc')

#FIXME decide what to do (and below where return)
#  approxSd <- makeSd(pedigree, homogametic = homogametic, DosageComp = DosageComp,
#	parallel = parallel, ncores = ncores,
#	invertSd = invertSd, returnS = returnS, det = det)

  lapproxSd <- summary(Sdsim)
  names(lapproxSd)[3L] <- "simSd"
  listSdsim <- NULL
  if(calcSE) {
#TODO check that this is the same for sex-chromosomal case as it is for autosomes
## Could see how it would differ, because difference chances of inheriting certain alleles
## Not strictly binomial sampling
     lapproxSd$Sdse <- vapply(lapproxSd$simSd, FUN = function(x, N){(sqrt(x * (1 - x))) / sqrt(N)}, FUN.VALUE = vector("numeric", 1), N)
     listSdsim <- lapproxSd
  } 

  cat(".done", "\n")
  logDetSdsim <- determinant(Sdsim, logarithm = TRUE)$modulus[1]
  
  if(invertSd){
    cat("inverting Sdsim ...")
    Sdsiminv <- solve(Sdsim)
      Sdsiminv@Dimnames <- Sdsim@Dimnames
    cat(".done", "\n")
    listSdsiminv <- sm2list(Sdsiminv, rownames = Sdsim@Dimnames[[1L]], colnames = c("row", "column", "simSdinverse"))
    Sdsim <- as(Sdsim, "dgCMatrix")
#FIXME decide what to do with approximate objects
#    return(list(S = approxSd$S, Sd = approxSd$Sd, logDetSd = approxSd$logDet, Sdinv = approxSd$Sdinv, listSdinv = approxSd$listSdinv, Sdsim = Sdsim, logDetSdsim = logDetSdsim, Sdsiminv = Sdsiminv, listSdsim = listSdsim, listSdsiminv = listSdsiminv))
#  } else{
#      return(list(S = approxSd$S, Sd = approxSd$Sd, logDetSd = approxSd$logDet, Sdsim = Sdsim, logDetSdsim = logDetSdsim, listSdsim = listSdsim))
#    } 
    return(list(Sdsim = Sdsim, logDetSdsim = logDetSdsim, Sdsiminv = Sdsiminv, listSdsim = listSdsim, listSdsiminv = listSdsiminv))
  } else{
      return(list(Sdsim = Sdsim, logDetSdsim = logDetSdsim, listSdsim = listSdsim))
    } 


}




################################################################################

#pedFS <- simPedHS(1, 1, 4)
#makeSdsim(pedFS, heterogametic = "M", N = 10000, returnS = TRUE)
#SdOut <- makeSd(pedigree, heterogametic = "M", returnS = TRUE)

#pedDFC <- simPedDFC(1)
#Sd <- makeSdsim(pedDFC, heterogametic = "M", N = 1000000, returnS = TRUE)$Sdsim

#makeSdsim(FG90, heterogametic = "0", N = 100000, returnS = TRUE)$Sdsim




