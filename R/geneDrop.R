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

