
#' @aliases makeGGAinv
#' @rdname makeAinv
#' @export
makeGGAinv <- function(pedigree, f = NULL, ggroups = NULL, det = TRUE, ...){
  zerNA <- union(which(pedigree[, 2] == "0"), which(pedigree[, 3] == "0"))
  astNA <- union(which(pedigree[, 2] == "*"), which(pedigree[, 3] == "*"))
  naNA <- union(which(is.na(pedigree[, 2])), which(is.na(pedigree[, 3])))
  naPed <- union(union(zerNA, astNA), naNA)

  ####    pedigree type "D"     ####
  if(!is.numeric(ggroups) && length(ggroups) == length(unique(ggroups))){
    ptype <- "D"
    if(length(naPed) > 0){
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
  N <- nrow(nPed)
  eN <- N - nggroups



  ### Make Q (adapted - with changes - from `ggcontrib()` code) #########
  dnmiss <- which(nPed[, 2] != -998)
  snmiss <- which(nPed[, 3] != -998)
  maxcnt <- (length(dnmiss) + length(snmiss) + N)
  Tinv.row <- Tinv.x <- rep(0, maxcnt)
  Tinv.col <- rep(0, N+1)

  Cout <- .C("reT", PACKAGE = "nadiv",
	as.integer(nPed[, 2] - 1),
	as.integer(nPed[, 3] - 1),
        as.integer(Tinv.row),
	as.integer(Tinv.col),
	as.double(Tinv.x),
	as.integer(maxcnt),
	as.integer(N),
	as.double(c(0.5, 0.5, 1.0, 1.0))) #maternal, paternal, self, diagonal

  Tinv <- t(sparseMatrix(i = Cout[[3]][1:Cout[[6]]],
    p = Cout[[4]], 
    x = Cout[[5]][1:Cout[[6]]],
    dims = c(N, N),
    dimnames = list(as.character(nPed[, 1]), as.character(nPed[, 1])),
    symmetric = FALSE, index1 = FALSE))

  T <- as(solve(Tinv), "CsparseMatrix")
  T@Dimnames <- Tinv@Dimnames
  Q <- as(T[-c(1:nggroups), 1:nggroups], "matrix")

  
  ### Make Ainv_j for each group  #########
  # Make non-group specific objects: create pedigree without genetic groups
  nPedNoGG <- ronPed((nPed - nggroups), -c(1:nggroups))
    nPedNoGG[nPedNoGG[, 2] <= 0, 2] <- -998
    nPedNoGG[nPedNoGG[, 3] <= 0, 3] <- -998
  Tinv2 <- makeTinv(nPedNoGG)
  OneMinDii <- 1 - makeDiiF(nPedNoGG)$D@x  


  ##### Begin group-specific calculations ######
  Ainv_list <- vector("list", length = nggroups)
  if(det){
    logDet_list <- vector("numeric", length = nggroups) 
  } else logDet_list <- rep(NULL, nggroups)
  if(any(renPed[1:nggroups] != 1:nggroups)){
    stop("Somehow, genetic groups were not put at the beginning of the pedigree when it was re-ordered by the function. Contact the package maintainer")
  }
  fsOrd <- as(as.integer(renPed[-c(1:nggroups)] - nggroups), "pMatrix")

  for(j in 1:nggroups){
    # Scaling factor for T by multiplying diagonals by Q[, j]
    ## However, from Muff et al. 2019 (eqn 8) and (eqn 12), pull scaling factor
    ### out for both T matrices in A=TDT and will use calculate DTilde
    #### Operation allows the use of same T (or Tinv) for all groups
    scalingj <- ifelse(Q[, j] > 0, Q[, j]^2, 1e-12)
    ##TODO could make option of calculating group-specific inbreeding to make `Dj`
    #### (would follow Lacy, Alaks, & Walsh 1996 founder-specific F approach) 
    #### this would be implemented with Muff et al. 2019 (eqn 11)
    ## However, currently do Muff et al. 2019 approximation to D_j (eqn 10)
    Dj <- (1 - Q[, j] * OneMinDii)
    DinvTilde <- Diagonal(x = 1 / (Dj * scalingj), n = eN) 
    # Ainv_j created by Muff et al. 2019 (eqn 13)
    Ainv_list[[j]] <- structure(as(crossprod(fsOrd, crossprod(sqrt(DinvTilde) %*% Tinv2)) %*% fsOrd, "CsparseMatrix"),
                                geneticGroups = c(0, 0))
    if(det) logDet_list[j] <- -1*determinant(Ainv_list[[j]], logarithm = TRUE)$modulus[1]

    if(ptype == "D"){
      names(Ainv_list)[j] <- as.character(ggroups[j])
      Ainv_list[[j]]@Dimnames <- list(as.character(pedigree[, 1]), NULL)
      if(det) names(logDet_list)[j] <- as.character(ggroups[j])
    } else {
        names(Ainv_list)[j] <- as.character(pedigree[groupRows[j], 1])
        Ainv_list[[j]]@Dimnames <- list(as.character(pedigree[-groupRows, 1]), NULL)
        if(det) names(logDet_list)[j] <- as.character(pedigree[groupRows[j], 1])
      }
  }  #<-- end for loop

  # re-order Q to original/pedigree order
  Q[] <- Q[invPerm(fsOrd@perm), ]
  if(ptype == "D"){
    attr(Q, "dimnames") <- list(as.character(pedigree[, 1]),
				as.character(ggroups))
  } else{
      attr(Q, "dimnames") <- list(as.character(pedigree[-groupRows, 1]),
				as.character(pedigree[groupRows, 1]))
    }


 return(list(Ainv = Ainv_list,
	listAinv = lapply(Ainv_list, FUN = function(j){ structure(sm2list(j, rownames = rownames(j), colnames = c("row", "column", "Ainv")), geneticGroups = c(0, 0))}),
	logDet = logDet_list,
	Q = Q))
}

