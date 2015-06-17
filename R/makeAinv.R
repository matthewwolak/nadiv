makeAinv <- function(pedigree, det = FALSE){
  renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
  nPed <- numPed(pedigree[renPed, ])
  N <- dim(nPed)[1]
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  inbreeding <- c(rep(0, N), -1)
  dii <- rep(0, N)

  sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	x = as.integer(1, length(Tinv.row)), index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))
  Ainv <- t(crossprod(crossprod(Diagonal(n = N, x = TRUE), sTinv))) # transpose gives lower triangle
  # 1: Adds Ainv elements in same for loop as calculation of f
  # 2: First checks to see if individual k has same dam and sire as k-1, if so then just assigns k-1's f 
  # 3: simplifies the calculation of the addition to the Ainv element (instead of alphai * 0.25 - defines alphai=alphai*0.25).
  nPed[nPed == -998] <- N + 1
  Cout <- .C("ainvML",
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
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedigree[, 1], colnames = c("row", "column", "Ainv")), f = Cout[[3]][-(N+1)], logDet = logDet))
}

