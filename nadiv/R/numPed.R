################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################

numPed<-function(pedigree, check = TRUE)
{
 if(check){      
  if(length(which(pedigree[, 2] == 0)) > 0){
    pedigree[which(pedigree[, 2] == 0), 2] <- NA
    warning("Zero in the dam column interpreted as a missing parent")
  }
  if(length(which(pedigree[, 3] == 0)) > 0){
    pedigree[which(pedigree[, 3] == 0), 3] <- NA
    warning("Zero in the sire column interpreted as a missing parent")
  }
  if(length(which(pedigree[,2] == "*")) > 0) pedigree[which(pedigree[, 2] == "*"), 2] <- NA
  if(length(which(pedigree[,3] == "*")) > 0) pedigree[which(pedigree[, 3] == "*"), 3] <- NA

  if(all(is.na(pedigree[, 2])) & all(is.na(pedigree[, 3]))){
     stop("All dams and sires are missing")
  }
  if(dim(pedigree)[2] != 3){
     stop("pedigree must have three columns: ID, Dam and Sire")
  }
  if(sum((na.omit(pedigree[, 2]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 2]) == FALSE)){
     stop("individuals appearing as dams but not in pedigree: first use the 'prepPed' function")
  }
  if(sum((na.omit(pedigree[, 3]) %in% pedigree[, 1]) == FALSE) > 0 & any(is.na(pedigree[, 3]) == FALSE)){
     stop("individuals appearing as sires but not in pedigree: first use the 'prepPed' function")
  }
  if(sum(duplicated(pedigree[, 1])) > 0){
     stop("some individuals appear more than once in the pedigree")
  }
 }
  numped <- matrix(as.integer(-998), dim(pedigree)[1], dim(pedigree)[2])
  numped[, 1] <- as.integer(seq(1, dim(pedigree)[1], 1))
  numped[, 2] <- match(pedigree[, 2], pedigree[, 1], nomatch = -998)
  numped[, 3] <- match(pedigree[, 3], pedigree[, 1], nomatch = -998)
  dnmiss <- which(numped[, 2] != -998)
  snmiss <- which(numped[, 3] != -998)
  bnmiss <- which(numped[, 2] != -998 & numped[, 3] != -998)
 if(check){
  if(length(intersect(numped[, 2][dnmiss], numped[, 3][snmiss])) > 0 & (length(dnmiss) > 0) & (length(snmiss) > 0)){
      warning("Dams appearing as Sires - assumed selfing in pedigree")
  }
  if(any(numped[, 2][dnmiss] > numped[, 1][dnmiss]) & (length(dnmiss) > 0)){
     stop("Offspring appearing before their dams: first use the 'prepPed' function")
  }
  if(any(numped[, 3][snmiss] > numped[, 1][snmiss]) & (length(snmiss) > 0)){
     stop("Offspring appearing before their Sires: first use the 'prepPed' function")
  }
 }
 numped
}


