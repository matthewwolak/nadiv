################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################
pedigree <- Mrode2; algorithm <- "ML"; det <- FALSE
#makeAinv2 <- function(pedigree, algorithm = c("H", "ML"), det = FALSE){
  numped <- numPed(pedigree)
  N <- dim(numped)[1]
  dnmiss <- which(numped[, 2] != -998)
  snmiss <- which(numped[, 3] != -998)
  bnmiss <- which(numped[, 2] != -998 & numped[, 3] != -998)
  Tinv.row <- c(numped[, 1][dnmiss], numped[, 1][snmiss], 1:N)
  Tinv.col <- c(numped[, 2][dnmiss], numped[, 3][snmiss], 1:N)
  Tinv.x <- c(rep(-0.5, length(dnmiss) + length(snmiss)), rep(1, N))
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  Tinv <- Matrix(0, N, N, sparse = TRUE)
  Tinv[1, 2] <- 1
  Tinv@i <- as.integer(Tinv.row[el.order] - 1)
  Tinv@p <- as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1)
  Tinv@x <- as.double(Tinv.x[el.order])
  nA <- N + length(dnmiss) + length(snmiss)
  nA <- nA + sum(duplicated(paste(numped[, 2], numped[, 3])[bnmiss]) == FALSE)
  inbreeding <- c(rep(0, N), -1)
  numped[numped == -998] <- N + 1
  al2use <- match.arg(algorithm, choices = c("H", "ML"))
  if(al2use == "H"){
    Cout <- .C("acinv",
	    as.integer(numped[, 2] - 1), #dam
	    as.integer(numped[, 3] - 1),  #sire
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



   if(al2use == "ML"){
      Ainv.row <- c(seq(N), dnmiss, snmiss, ) 
#     dii <- rep(0, N)
     Cout <- .C("acinvML",
	    as.integer(numped[, 2] - 1), #dam
	    as.integer(numped[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
            as.double(rep(0, N)),  #dii
#            as.integer(Tinv@i),  #iTinvP
#	    as.integer(c(Tinv@p, length(Tinv@x))),  #pTinvP
#            as.double(Tinv@x),  #xTinvP
            as.integer(N),   #n
#	    as.integer(length(Tinv@x)), #nzmaxTinvP

            as.double(rep(0, 6*N)))  # xA
	    as.integer(rep(0, nA)), #iA
	    as.integer(rep(0, N + 1)), #pA
	    as.integer(nA)) #nzmaxA
# Triplets are lower triangle!
#order of triplets in xA: i,i; i,sire; i,dam; sire,sire; sire,dam/dam,sire; dam,dam;  
#     Ainv <- Matrix(0, N, N, dimnames = list(as.character(pedigree[, 1]), NULL))
#     Ainv[1, 2] <- 1
#     Ainv@i <- Cout[[9]][1:Cout[[12]]]
#     Ainv@p <- Cout[[10]]
#     Ainv@x <- Cout[[11]][1:Cout[[12]]]

tmp <- cbind(i = c(rep(seq(N), 3), numped[, 3], apply(numped[, 2:3], MARGIN = 1, FUN = max), numped[, 2]),
	j = c(seq(N), numped[, 3], numped[, 2], numped[, 3], apply(numped[, 2:3], MARGIN = 1, FUN = min), numped[, 2]),
	x = Cout[[6]])[tmp[, 1] != N+1 & tmp[, 2] != N+1, ]

      Ainv <- sparseMatrix(i = tmp[, 1], j = tmp[, 2], x = tmp[, 3])
#     Ainv <- as(Ainv, "dgCMatrix")
   }


   if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedigree[, 1], colnames = c("row", "column", "Ainv")), f = Cout[[3]][-(N+1)], logDet = logDet))
}

