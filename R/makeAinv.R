makeAinv <- function(pedigree, ggroups = NULL, fuzz = NULL, keepPhantoms = TRUE, gOnTop = FALSE, det = FALSE){
  if(is.null(ggroups)){
     ptype <- "O"
     renPed <- order(genAssign(pedigree), pedigree[, 2], pedigree[, 3], na.last = FALSE)
     nPed <- numPed(pedigree[renPed, ])
     groupRows <- NULL
     nggroups <- 0
  } else {
    if(length(ggroups) == dim(pedigree)[1]){
       stop("length(ggroups) should either be:\n  == 1 and a numeric indicating the number of genetic groups (type 'A')\n  == length(unique(ggroups)) and a character vector indicating the names of each unique genetic goup (type 'D')")
    }
    if(!is.null(fuzz)) stop("fuzzy genetic groups not yet implemented: fuzz must be NULL")
    zerNA <- union(which(pedigree[, 2] == "0"), which(pedigree[, 3] == "0"))
    astNA <- union(which(pedigree[, 2] == "*"), which(pedigree[, 3] == "*"))
    naNA <- union(which(is.na(pedigree[, 2])), which(is.na(pedigree[, 3])))
    naPed <- union(union(zerNA, astNA), naNA)

    ####    pedigree type "D"     ####
    if(!is.numeric(ggroups) && length(ggroups) == length(unique(ggroups))){
       ptype <- "D"
       if(is.null(fuzz) && length(naPed) > 0){
          stop("When supplying a vector of unique genetic group names in the 'ggroups' argument, all individuals in the pedigree must have a genetic group when a parent is unknown (<NA>, '0' and '*' are considered unknown parents)")# or be identified as a phantom parent in 'fuzz'")
       }
       if(any(is.na(ggroups))) ggroups <- na.omit(ggroups)
       pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])), dam = c(rep(NA, length(ggroups)), as.character(pedigree[, 2])), sire = c(rep(NA, length(ggroups)), as.character(pedigree[, 3])))
     #TODO: write method for numPed to handle genetic group pedigree
     ## same genetic group for dam and sire won't throw warning about selfing
       nPed <- suppressWarnings(numPed(pedalt))    # FIXME: remove suppressWarnings() throughout file
     ## return command below as attribute
       groupRows <- nPed[which(nPed[, 2] == -998), 1]
    }

  ####    pedigree type "A"     ####
    if(length(ggroups) == 1){
       ptype <- "A"
          if(length(naPed) == ggroups){
             nPed <- suppressWarnings(numPed(pedigree))
             groupRows <- naPed
          } else {
             stop("Only rows identifying genetic groups should have missing parents.\n  All individuals in the pedigree must have a genetic group when a parent is unknown")# or be identified as a phantom parent in 'fuzz'")
         }
    }

  nggroups <- length(groupRows)
  renPed <- order(suppressWarnings(genAssign(nPed)), nPed[, 2], nPed[, 3])
  nPed <- numPed(ronPed(nPed, renPed))
  }
  ####
  N <- nrow(nPed)
  eN <- N - nggroups
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  Tinv.row <- c(nPed[, 1][dnmiss], nPed[, 1][snmiss], 1:N)
  Tinv.col <- c(nPed[, 2][dnmiss], nPed[, 3][snmiss], 1:N)
  el.order <- order(Tinv.col + Tinv.row/(N + 1), decreasing = FALSE)
  if(is.null(ggroups)){
     sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))
  } else {
     sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	p = as.integer(c(match(1:N, Tinv.col[el.order]), length(el.order) + 1) - 1),
	index1 = FALSE, dims = c(N, N), symmetric = FALSE,
	dimnames = list(as.character(nPed[, 1]), NULL))[-groupRows, ]
  }
  Ainv <- t(crossprod(sTinv)) # transpose gives lower triangle
  # 1: Adds Ainv elements in same for loop as calculation of f
  # 2: First checks to see if individual k has same dam and sire as k-1, if so then just assigns k-1's f 
  # 3: simplifies the calculation of the addition to the Ainv element (instead of alphai * 0.25 - defines alphai=alphai*0.25).
  nPed[nPed == -998] <- N + 1
  Cout <- .C("ainvml",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(c(rep(-1, nggroups), rep(0, eN), -1)),	#f
            as.double(rep(0, N)),  				#dii
            as.integer(N),   					#n
            as.integer(nggroups),   				#g
            as.double(rep(0, length(Ainv@i))),  		#xA
	    as.integer(Ainv@i), 				#iA
	    as.integer(Ainv@p), 				#pA
	    as.integer(length(Ainv@i))) 			#nzmaxA
  Ainv <- as(Ainv, "dsCMatrix")
  Ainv@x <- Cout[[7]]
  fsOrd <- as(as.integer(renPed), "pMatrix")
  Ainv <- as(t(fsOrd) %*% Ainv %*% fsOrd, "dgCMatrix")
   if(ptype == "D"){
      Ainv@Dimnames <- list(as.character(pedalt[, 1]), NULL)
   } else {
      Ainv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
     }
  if(!is.null(ggroups) && !gOnTop){ 
     permute <- as(as.integer(c(seq(eN+1, N, 1), seq(eN))), "pMatrix")
     Ainv <- t(permute) %*% Ainv %*% permute
  }
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = Ainv, listAinv = sm2list(Ainv, rownames = rownames(Ainv), colnames = c("row", "column", "Ainv")), f = Cout[[3]][seq(nggroups+1, N, 1)], logDet = logDet))
}

