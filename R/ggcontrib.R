ggcontrib <- function(pedigree, ggroups = NULL, fuzz = NULL, output = "matrix"){

  if(is.null(ggroups)){
     ptype <- "A"
     numped <- numPed(pedigree)
     ggroupsD <- as.character(pedigree[which(numped[, 2] == -998), 1])
     ggroupsS <- as.character(pedigree[which(numped[, 3] == -998), 1])
     if(!all(ggroupsD == ggroupsS)){
        stop("Only rows identifying genetic groups should have missing parents.  All individuals in the pedigree must have a genetic group when a parent is unknown")
     }
     nggroups <- numped[which(numped[, 2] == -998), 1]

  } else{
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
               if (length(which(pedigree[, 2] == "*")) > 0){ 
                  pedigree[which(pedigree[, 2] == "*"), 2] <- NA
               }
               if (length(which(pedigree[, 3] == "*")) > 0){
                  pedigree[which(pedigree[, 3] == "*"), 3] <- NA
               }   
               pedalt <- data.frame(id = I(as.character(pedigree[, 1])), dam = I(as.character(pedigree[, 2])), sire = I(as.character(pedigree[, 3])))
               pedalt[which(is.na(pedigree[, 2])), 2] <- as.character(ggroups[which(is.na(pedigree[, 2]))])
               pedalt[which(is.na(pedigree[, 3])), 3] <- as.character(ggroups[which(is.na(pedigree[, 3]))])
               uggroups <- as.character(unique(ggroups))
               uggroups <- uggroups[!is.na(uggroups)]
               pedalt <- data.frame(id = c(uggroups, as.character(pedalt[, 1])), dam = c(rep(NA, length(uggroups)), as.character(pedalt[, 2])), sire = c(rep(NA, length(uggroups)), as.character(pedalt[, 3])))
               numped <- suppressWarnings(numPed(pedalt))
               nggroups <- numped[which(numped[, 2] == -998), 1]       
            }
        }


  N <- dim(numped)[1]
  dnmiss <- which(numped[, 2] != -998)
  snmiss <- which(numped[, 3] != -998)
  maxcnt <- (length(dnmiss) + length(snmiss) + N)
  Tinv.row <- Tinv.x <- rep(0, maxcnt)
  Tinv.col <- rep(0, N+1)

  Cout <- .C("reT",
		as.integer(numped[, 2] - 1),
		as.integer(numped[, 3] - 1),
                as.integer(Tinv.row),
		as.integer(Tinv.col),
		as.double(Tinv.x),
		as.integer(maxcnt),
		as.integer(N),
		as.double(c(0.5, 0.5, 1.0, 1.0))) #maternal, paternal, self, diagonal

# fuzz parameter should alter the entries below to be the input as opposed to -0.5
  Tinv <- Matrix(0, N, N, sparse = TRUE)
  Tinv[1, 2] <- 1
  Tinv@i <- Cout[[3]][1:Cout[[6]]]
  Tinv@p <- Cout[[4]]
  Tinv@x <- Cout[[5]][1:Cout[[6]]]

  T <- as(solve(t(Tinv)), "dgCMatrix")
  if(ptype == "A"){
     T@Dimnames <- list(pedigree[, 1], pedigree[, 1])
  } else{
       T@Dimnames <- list(pedalt[, 1], pedalt[, 1])   
    }
 return(as(T[-c(nggroups), nggroups], output))
}
