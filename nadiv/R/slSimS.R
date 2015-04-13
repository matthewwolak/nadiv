#This code simulates genotypes from single locus within a pedigree for just the X chromosomes
slSimS <- function(pedigree, nloci, xvar, heterogametic = "0"){
  numped <- numPed(pedigree[, 1:3])
  if(length(which(numped[, 2] == -998 & numped[, 3] != -998)) > 0){
     stop("Individuals must either have both parents known or no known parents")
  }
  if(length(which(numped[, 2] != -998 & numped[, 3] == -998)) > 0){
     stop("Individuals must either have both parents known or no known parents")
  }
  n <- dim(numped)[1]
  sex <- vector("numeric", length = n)
  sex[which(pedigree[, 4] == heterogametic)] <- 1
  sex[which(pedigree[, 4] != heterogametic)] <- 2
  hap1 <- matrix(0, nrow = nloci, ncol = n)
  hap2 <- matrix(0, nrow = nloci, ncol = n)
  founders <- which(numped[, 2] == -998 & numped[, 3] == -998)
  hap1[, founders] <- rnorm(length(founders)*nloci, 0, sqrt(0.5*(xvar/nloci)))
  hap2[, founders[which(sex[founders] == 2)]] <- rnorm(length(founders[which(sex[founders] == 2)])*nloci, 0, sqrt(0.5*(xvar/nloci)))
  segregation <- ceiling(runif((nloci*(n - length(founders))), min = 0, max = 2))

  Cout <- .C("slsims",
		as.integer(numped[, 2] - 1),
		as.integer(numped[, 3] - 1),
		as.integer(sex),
		as.integer(n),
		as.integer(nloci),
		as.integer(segregation),
		as.double(hap1),
		as.double(hap2))
  hap1 <- matrix(Cout[[7]], nrow = nloci, ncol = n)
  hap2 <- matrix(Cout[[8]], nrow = nloci, ncol = n)

return(list(hap1 = hap1, hap2 = hap2)) 
}


