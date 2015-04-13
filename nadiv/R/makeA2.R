makeA2 <- function(pedigree, f = NULL)
{
  numped <- numPed(pedigree)
  N <- dim(numped)[1]
####
#  nA <- N + 2 * length(dnmiss) + 2 * length(snmiss)
#  nA <- nA + 2 * sum(duplicated(paste(numped[, 2], numped[, 3])[bnmiss]) == FALSE)
nA <- (N^2)#/ 2
####
  inbreeding <- if(is.null(f)) c(rep(0, N), -1) else c(f, -1)
  dii <- rep(1, N)
  numped[numped == -998] <- N + 1
    Cout <- .C("makea",
	    as.integer(numped[, 2] - 1), #dam
	    as.integer(numped[, 3] - 1),  #sire
	    as.double(inbreeding),  #f
	    as.double(dii),  #dii
            as.integer(N),   #N
	    as.integer(rep(0, nA)), #iAP
	    as.integer(rep(0, N + 1)), #pAP
	    as.double(rep(0, nA)), #xAP
	    as.integer(nA)) #nzmaxAP

   A <- Matrix(0, N, N)
   A[1, 2] <- 1
   A@i <- Cout[[6]][1:Cout[[9]]]
   A@p <- Cout[[7]]
   A@x <- Cout[[8]][1:Cout[[9]]]
as(A, "symmetricMatrix")
}

