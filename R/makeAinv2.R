################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################
#rm(list = ls()); pedigree <- Mrode2; algorithm <- "MLb"; det <- FALSE
makeAinv2 <- function(pedigree, algorithm = c("H", "MLa", "MLb"), det = FALSE){
  al2use <- match.arg(algorithm, choices = c("H", "MLa", "MLb"))
  if(al2use == "MLb"){
     renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
     nPed <- numPed(pedigree[renPed, ])
  } else nPed <- numPed(pedigree)
  N <- dim(nPed)[1]
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  Tinv.x <- c(rep(-0.5, length(dnmiss) + length(snmiss)), rep(1, N))
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  inbreeding <- c(rep(0, N), -1)
  dii <- rep(0, N)


  if(al2use == "H"){
     bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
     Tinv <- Matrix(0, N, N, sparse = TRUE)
     Tinv[1, 2] <- 1
     Tinv@i <- as.integer(Tinv.row[el.order] - 1)
     Tinv@p <- as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1)
     Tinv@x <- as.double(Tinv.x[el.order])
     nA <- N + length(dnmiss) + length(snmiss)
     nA <- nA + sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
     nPed[nPed == -998] <- N + 1
     Cout <- .C("acinv",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
            as.integer(Tinv@i),  #iTinvP
	    as.integer(c(Tinv@p, length(Tinv@x))),  #pTinvP
            as.double(Tinv@x),  #xTinvP
            as.integer(N),   #nTinvP
	    as.integer(length(Tinv@x)), #nzmaxTinvP
	    as.integer(rep(0, nA)), #iAP
	    as.integer(rep(0, N + 1)), #pAP
	    as.double(rep(0, nA)), #xAP
	    as.integer(nA)) #nzmaxAP

     Ainv <- Matrix(0, N, N, dimnames = list(as.character(pedigree[, 1]), NULL))
     Ainv[1, 2] <- 1
     Ainv@i <- Cout[[9]][1:Cout[[12]]]
     Ainv@p <- Cout[[10]]
     Ainv@x <- Cout[[11]][1:Cout[[12]]]
     Ainv <- as(Matrix::tcrossprod(Ainv), "dgCMatrix")
   }



   if(al2use == "MLa" | al2use == "MLb"){
      sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	x = TRUE, index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))
      nPed[nPed == -998] <- N + 1
      # transpose gives lower triangle
      Ainv <- t(crossprod(crossprod(Diagonal(n = N, x = TRUE), sTinv)))
      #system.time(Ainv2 <- t(as(crossprod(as(crossprod(Diagonal(n = N, x = TRUE), sTinv), "lgCMatrix")), "lsCMatrix")))
#identical(Ainv@i, Ainv2@i)
#identical(Ainv@p, Ainv2@p)
   }
if(al2use == "MLa"){
# differs from H algorithm:
  # 1: doesn't use matrix multiplication (Tinv %*% Dinv %*% tTinv) to create Ainv, instead assigns element-by-element the non-zero Ainv entries (will be useful when adding genetic groups to the algorithm)
     Cout <- .C("acinvMLa",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
            as.double(rep(0, N)),  #dii
            as.integer(N),   #n
            as.double(rep(0, length(Ainv@x))),  # xA
	    as.integer(Ainv@i), #iA
	    as.integer(Ainv@p), #pA
	    as.integer(length(Ainv@x))) #nzmaxA
     Ainv@x <- Cout[[6]]
     Ainv <- as(Ainv, "dgCMatrix")
}

if(al2use == "MLb"){
# differs from MLa: 
  # 1: Adds Ainv elements in same for loop as calculation of f
  # 2: First checks to see if individual k has same dam and sire as k-1, if so then just assigns k-1's f 
  # 3: simplifies the calculation of the addition to the Ainv element (instead of alphai * 0.25 - defines alphai=alphai*0.25).

#warning("Not producing correct Ainv as of 2pm on 8 June 2015")
cat("10:24am\n")
     Cout <- .C("acinvMLb",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
            as.double(rep(0, N)),  #dii
            as.integer(N),   #n
            as.double(rep(0, length(Ainv@x))),  # xA
	    as.integer(Ainv@i), #iA
	    as.integer(Ainv@p), #pA
	    as.integer(length(Ainv@x))) #nzmaxA
     Ainv@x <- Cout[[6]]
     fsOrd <- as(as.integer(renPed), "pMatrix")
     Ainv <- as(t(fsOrd) %*% Ainv %*% fsOrd, "dgCMatrix")


}

   if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL
 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedigree[, 1], colnames = c("row", "column", "Ainv")), f = Cout[[3]][-(N+1)], dii = dii, logDet = logDet))
}


#library(devtools); dev_mode(TRUE); install.packages("~/Dropbox/nadiv_parent/nadiv_2.13.4.tar.gz", repos = NULL); library(nadiv); library(asreml)

#ped <- Mrode2
#ped <- warcolak[, 1:3]
#ped <- read.table("~/Dropbox/GarlandMice/pedigree_fixed_switched.txt", header = TRUE)
#ped <- read.table("~/Dropbox/DataMenagerie/Liedvogel2012_TitPeds/BlueTit.txt", header = TRUE)
#ped <- read.table("~/Dropbox/DataMenagerie/Liedvogel2012_TitPeds/GreatTit.txt", header = TRUE)

#ped <- data.frame(id = seq(10), dam = c(NA, NA, 2, 2, 4, 5, 2, 2, 5, NA), sire = c(NA, NA, 1, 1, 1, 3, 3, 3, 6, 7)) 

#system.time(Ainv <- makeAinv(ped)$Ainv)
#system.time(AinvH <- nadiv:::makeAinv2(ped, "H")$Ainv)
#system.time(AinvMLa <- nadiv:::makeAinv2(ped, "MLa")$Ainv)
#system.time(AinvMLb <- nadiv:::makeAinv2(ped, "MLb")$Ainv)
#system.time(AinvA <- asreml.Ainverse(ped)$ginv)


#sum(Ainv - AinvH); range(Ainv - AinvH)

#sum(Ainv - AinvMLa); range(Ainv - AinvMLa)

#sum(Ainv - AinvMLb); range(Ainv - AinvMLb)

#sum(AinvMLa - AinvMLb); range(AinvMLa - AinvMLb)



