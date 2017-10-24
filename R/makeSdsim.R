makeDsim <- function(pedigree, heterogametic, N,
	parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertSd = TRUE, calcSE = FALSE, returnS = FALSE){

  approxSd <- makeSd(pedigree, heterogametic = heterogametic, DosageComp = DosageComp,
	parallel = parallel, ncores = ncores,
	invertSd = invertSd, returnS = returnS, det = TRUE)
  lapproxSd <- summary(approxSd$Sd)

  nPed <- numPed(pedigree)
  n <- dim(pedigree)[1]
  alleles <- matrix(as.integer(-998), nrow = n, ncol=2) 
  dfounders <- which(nPed[, 2] == -998)
  sfounders <- which(nPed[, 3] == -998)
  uniqp <- c(unique(nPed[, 2])[-1], unique(nPed[, 3])[-1])
  ndfounders <- length(dfounders)
  
  alleles[dfounders, 1] <- as.integer(seq(1, ndfounders, 1)) 
  alleles[sfounders, 2] <- as.integer(seq(ndfounders+1, (ndfounders + length(sfounders)), 1))
  dalleles <- rep(alleles[, 1], each = N)
  salleles <- rep(alleles[, 2], each = N)
  
  cat("making Sdsim ...")

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


