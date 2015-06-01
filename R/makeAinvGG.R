################################################
#Adapted from part of the 'inverseA' function
# written by Jarrod Hadfield
#in the 'MCMCglmm' package
################################################
# First, entirely R code version, followed by implementation in C++
## Note: this is from Quaas 1988 which extends Henderson's method of obtaining A-inverse
## Note: also, this is 'backwards' from Quaas - where I put the groups on top of the list of equations and the normal pedigree below (this mixes up the modified MME).
#TODO: Could consider reducing steps to be unique for each unique family (i.e., do all Full-sibs at once)

library(nadiv)
#rm(list = ls()); source("~/Dropbox/nadiv_parent/makeQuaas1988Fig1Ped.R"); pedigree <- Q1988_2; ggroups <- groups; fuzz <- NULL; det <- FALSE  
#TODO: I THINK THIS ONLY WORKS WHEN NO INBREEDING OCCURS! (hence why DMU has only implemented non-inbreeding scenario) - but see top of p. 4 of Fikse 2009: x=d^-1
makeAinvGGR <- function(pedigree, ggroups = NULL, fuzz = NULL, gOnTop = FALSE, det = FALSE){
  if(!is.null(ggroups)){
     if(length(ggroups) == length(unique(ggroups))){
        ptype <- "D"
        if(any(pedigree[, 2:3] == "0" | pedigree[, 2:3] == "*" | is.na(pedigree[, 2:3]))){
           stop("When specifying the unique genetic groups as a vector in the 'ggroups' argument, all individuals in the pedigree must have a genetic group when a parent is unknown (<NA>, '0' and '*' are considered unknown parents)")
        }
        pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])), dam = c(rep(NA, length(ggroups)), as.character(pedigree[, 2])), sire = c(rep(NA, length(ggroups)), as.character(pedigree[, 3])))
        numped <- suppressWarnings(numPed(pedalt))
        nggroups <- numped[which(numped[, 2] == -998), 1]
     }

     if(length(ggroups) == dim(pedigree)[1]){
        ptype <- "R"
        nonggped <- pedigree[which(is.na(ggroups)), 2:3]
        if(any(nonggped == "0" | nonggped == "*" | is.na(nonggped))){
           stop("All individuals with missing parents (indicated as '0', '*', or <NA>) must have a genetic group specified")
        }
        if(length(which(pedigree[, 2] == 0)) > 0){
           pedigree[which(pedigree[, 2] == 0), 2] <- NA
           warning("Zero in the dam column interpreted as a missing parent")
        }
        if(length(which(pedigree[, 3] == 0)) > 0){
           pedigree[which(pedigree[, 3] == 0), 3] <- NA
           warning("Zero in the sire column interpreted as a missing parent")
        }
        if(length(which(pedigree[, 2] == "*")) > 0){ 
           pedigree[which(pedigree[, 2] == "*"), 2] <- NA
        }
        if(length(which(pedigree[, 3] == "*")) > 0){
           pedigree[which(pedigree[, 3] == "*"), 3] <- NA
        }   
        pedalt <- data.frame(id = I(as.character(pedigree[, 1])), dam = I(as.character(pedigree[, 2])), sire = I(as.character(pedigree[, 3])))
        pedalt[which(is.na(pedigree[, 2])), 2] <- as.character(ggroups[which(is.na(pedigree[, 2]))])
        pedalt[which(is.na(pedigree[, 3])), 3] <- as.character(ggroups[which(is.na(pedigree[, 3]))])
        uggroups <- as.character(unique(ggroups))
        uggroups <- uggroups[!is.na(uggroups)]
        pedalt <- data.frame(id = c(uggroups, as.character(pedalt[, 1])), dam = c(rep(NA, length(uggroups)), as.character(pedalt[, 2])), sire = c(rep(NA, length(uggroups)), as.character(pedalt[, 3])))
        numped <- numPed(pedalt, check = FALSE)
        nggroups <- numped[which(numped[, 2] == -998), 1]       
     }
  }


  if(is.null(ggroups) & is.null(fuzz)){
     ptype <- "A"
     numped <- numPed(pedigree)
     ggroupsD <- as.character(pedigree[numped[, 2] == -998, 1])
     ggroupsS <- as.character(pedigree[numped[, 3] == -998, 1])
     if(!all(ggroupsD == ggroupsS)){
        stop("Only rows identifying genetic groups should have missing parents.  All individuals in the pedigree must have a genetic group when a parent is unknown")
     }
     nggroups <- numped[numped[, 2] == -998, 1]
  }

##################################
#browser()
##### FUZZY BITS!
  if(!is.null(fuzz)){
     if(!is.null(ggroups)) stop("Only one of ggroups or fuzz may be non-NULL")
#     if(nrow(fuzz) != 
  }
###########  END Fuzzy bits  ####
  N <- dim(numped)[1]
  n <- N - length(nggroups)
#  dnmiss <- which(numped[, 2] != -998)
#  snmiss <- which(numped[, 3] != -998)
#  maxcnt <- (length(dnmiss) + length(snmiss) + N)
#Could consider using the number of non-zero elements in A-inverse (and calculating number of non-zero elements in Q-inverse to preallocate Ainv to have the correct number of entries in Ainv@x, etc.) - use sparseMatrix(i, j, x) instead of Matrix(0, N, N)
   Ainv <- Matrix(0, N, N)   # An upper triangle sparse symmetric matrix
   for(i in (max(nggroups)+1):N){
       j <- numped[i, 3]   #sire or sire group (if i is base)
       k <- numped[i, 2]   # dam or dam group (if i is base)
       m <- sum(c(j, k) %in% nggroups)   # number of base parents of individual i
       x <- 4 / (m + 2)
       if(j == k && x == 1){
          Ainv[i, i] <- Ainv[i, i] + 1
          Ainv[j, j] <- Ainv[j, j] + 1
          Ainv[j, i] <- Ainv[j, i] + -1
       } else{
            Ainv[i, i] <- x
            Ainv[j, i] <- Ainv[j, i] + x*-0.5
	    Ainv[k, i] <- Ainv[k, i] + x*-0.5
            Ainv[j, j] <- Ainv[j, j] + x*0.25
            Ainv[k, k] <- Ainv[k, k] + x*0.25
            if(j <= k) Ainv[j, k] <- Ainv[j, k] + x*0.25 else Ainv[k, j] <- Ainv[k, j] + x*0.25
         }
   }
   Ainv <- as(forceSymmetric(Ainv), "dgCMatrix")
   Ainv@Dimnames <- list(pedalt[, 1], NULL)
   if(!gOnTop){ 
      permute <- as(as.integer(c(seq(n+1, N, 1), seq(n))), "pMatrix")
      Ainv <- t(permute) %*% Ainv %*% permute
   }
   if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

# return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedalt[, 1], colnames = c("row", "column", "Ainv")), f = Cout[[3]][-(N+1)], logDet = logDet))
#FIXME: Check rowname assignment for listAinv as well as any other derived output once the 'gOnTop' argument has done its magic
 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedalt[, 1], colnames = c("row", "column", "Ainv")), f = NULL, logDet = logDet))
}




##############################################
makeAinvGG <- function(pedigree, det = FALSE){
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

   Ainv <- Matrix(0, N, N)
   Ainv[1, 2] <- 1
   Ainv@i <- Cout[[9]][1:Cout[[12]]]
   Ainv@p <- Cout[[10]]
   Ainv@x <- Cout[[11]][1:Cout[[12]]]
   Ainv <- as(Matrix::tcrossprod(Ainv), "dgCMatrix")
   Ainv@Dimnames <- list(pedigree[, 1], NULL)
   if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = pedigree[, 1], colnames = c("row", "column", "Ainv")), f = Cout[[3]][-(N+1)], logDet = logDet))
}

