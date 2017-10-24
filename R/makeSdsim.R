makeDsim <- function(pedigree, heterogametic, N,
	parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertSd = TRUE, calcSE = FALSE, returnS = FALSE){

  if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")

  dc.model <- match.arg(DosageComp)
  if(is.null(dc.model)){
    warning("Assuming 'ngdc' dosage compensation model")
    dc.model <- "ngdc"
  }
  #TODO check if can have dominance under hori
  #TODO NO inbreeding effects in hets for hedo
  #TODO NO inbreeding effects in hoha
  if(dc.model == "hopi"){
    cat("Assume that sex chromosomal dominance allelic interactions do not occur under 'hopi'\n")
    return(NULL)
  }
if(dc.model != "ngdc"){ #FIXME temporarily only allow ngdc for now
  stop("Currently only supporting 'ngdc' model")
}

  n <- dim(pedigree)[1L]
  nPed <- numPed(pedigree)
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

  # diverges from `makeDsim()` and follows simplifications in `geneDrop()`
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  dalleles <- salleles <- vector("integer", length = n)

  dalleles[dfounders] <- as.integer(dfounders)
  salleles[sfounders] <- as.integer(sfounders)
  Ndalleles <- rep(dalleles, each = N)
  Nsalleles <- rep(salleles, each = N)
  
  cat("making Sdsim ...")

  # diversion to calculate maximum expected entries in sex-chromosome D matrix
  ## based on calculation for sex-chromosome S matrix (additive)
  dnmiss <- which(nPed[, 2] != -998 & sex == 1)
  snmiss <- which(nPed[, 3] != -998 & sex == 1)
  bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998 & sex == 1)
  nSd <- sum(sex) + 2 * length(dnmiss) + 2 * length(snmiss)
  nSd <- nSd + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)


#TODO LEFT OFF HERE 2017-10-24



  Cout <- .C("sdsim",
	as.integer(dalleles),
	as.integer(salleles),
	as.integer(N),
	as.integer(n),
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1),
	as.integer(approxSd$Sd@i),
	as.integer(approxSd$Sd@p),
	as.integer(rep(0, length(approxSd$Sd@i))))

  lapproxSd$simSd <- Cout[[9]] / N
  lapproxSd <- lapproxSd[which(lapproxSd[, 4] != 0), ]
  listSdsim <- NULL
  if(calcSE) {
     lapproxSd$Sdse <- vapply(lapproxSd$simSd, FUN = function(x, N){(sqrt(x * (1 - x))) / sqrt(N)}, FUN.VALUE = vector("numeric", 1), N)
     listSdsim <- lapproxSd
  } 

  Sdsim.row <- lapproxSd[,1]
  Sdsim.col <- lapproxSd[,2]
  Sdsim.x <- lapproxSd[,4]
  order.index <- order(Sdsim.col + Sdsim.row/(n+1), decreasing=FALSE)
  Sdsim <- Matrix(0, n, n, dimnames = list(as.character(pedigree[, 1]), NULL))
  Sdsim@uplo <- "U"
  Sdsim@i <- as.integer(Sdsim.row[order.index]-1)
  Sdsim@p <- as.integer(c(match(1:n, Sdsim.col[order.index]), length(order.index)+1)-1)
  Sdsim@x <- Sdsim.x[order.index]
  diag(Sdsim) <- diag(approxSd$Sd)
  cat(".done", "\n")
  logDetSdsim <- determinant(Sdsim, logarithm = TRUE)$modulus[1]
  
  if(invertSd){
    cat("inverting Sdsim ...")
    Sdsiminv <- solve(Sdsim)
      Sdsiminv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    cat(".done", "\n")
    listSdsiminv <- sm2list(Sdsiminv, rownames = pedigree[,1], colnames = c("row", "column", "simDinverse"))
    Sdsim <- as(Sdsim, "dgCMatrix")
    return(list(S = approxSd$S, Sd = approxSd$Sd, logDetSd = approxSd$logDet, Sdinv = approxSd$Sdinv, listSdinv = approxSd$listSdinv, Sdsim = Sdsim, logDetSdsim = logDetSdsim, Sdsiminv = Sdsiminv, listSdsim = listSdsim, listSdsiminv = listSdsiminv))
  } else{
      return(list(S = approxSd$S, Sd = approxSd$Sd, logDetSd = approxSd$logDet, Sdsim = Sdsim, logDetSdsim = logDetSdsim, listSdsim = listSdsim))
    } 

}


