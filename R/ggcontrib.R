#' Genetic group contribution
#' 
#' Calculates the genomic contribution each genetic group makes to every
#' individual in a pedigree
#' 
#' The specification of genetic groups is done in one of two approaches, either
#' using fuzzy classification or not.
#' 
#' Fuzzy classification enables phantom parents to be assigned to (potentially)
#' more than one genetic group (Fikse 2009). This method requires unique
#' phantom parent identities to be included in the pedigree for all observed
#' individuals with unknown parents. For 'p' phantom parents, 'p' identities
#' should be listed as individuals in the first 'p' rows of the pedigree and
#' these should be the only individuals in the pedigree with missing values in
#' their Dam and Sire columns (denoted by either 'NA', '0', or '*'). The matrix
#' supplied to the \code{fuzz} argument should have 'p' rows (one for each
#' phantom parent) and 'r' columns, where 'r' is the number of genetic groups.
#' 
#' Non-fuzzy classification can handle the specification of genetic groups in
#' three formats:
#' 
#' (1) similar to ASReml's format for specifying genetic groups, the first 'r'
#' rows of the pedigree (given to the \code{pedigree} argument) contain the
#' label for each genetic group in the ID column and indicate missing values
#' for the Dam and Sire columns (denoted by either 'NA', '0', or '*'). No
#' object is supplied to the \code{ggroups} argument. All individuals in the
#' pedigree must then have one of the 'r' genetic groups as parent(s) for each
#' unknown parent. Note, a warning message indicating \code{In
#' numPed(pedigree): Dams appearing as Sires} is expected, since the dam and
#' sire can be the same for all individuals in the pedigree composing the base
#' population of a genetic group.
#' 
#' (2) similar to Jarrod Hadfield's \code{rbv} function arguments in the
#' \code{MCMCglmm} package, for a pedigree of dimension i x 3 (given to the
#' \code{pedigree} argument), where 'i' is the total number of individuals in
#' the pedigree, a similar vector of length 'i' is given to the \code{ggroups}
#' argument. This vector lists either the genetic group to which each
#' individual's phantom parents belong or NA if the individual is not to be
#' considered part of one of the base populations (genetic groups). NOTE, this
#' approach does not allow phantom dams and phantom sires of a given individual
#' to be from different genetic groups.
#' 
#' (3) similar to DMU's format for specifying genetic groups. For a pedigree of
#' dimension i x 3 (given to the \code{pedigree} argument), where 'i' is the
#' total number of individuals in the pedigree, instead of missing values for a
#' parent, one of the 'r' genetic groups is specified. A character vector of
#' length 'r' with unique genetic group labels is given to the \code{ggroups}
#' argument. Note, that all individuals with a missing parent should have a
#' genetic group substituted instead of the missing value symbol (i.e., either
#' 'NA', '0', or '*').
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param ggroups An optional vector of either: genetic group assignment for
#'   every individual or just the unique genetic groups. \code{fuzz} must be
#'   \code{NULL} if an object is supplied to the \code{ggroups} argument.
#' @param fuzz A matrix containing the fuzzy classification of phantom parents
#'   into genetic groups. \code{ggroups} must be \code{NULL} if an object is
#'   supplied to the \code{fuzz} argument.
#' @param output Format for the output
#'
#' @return Returns i x r genetic group contributions to all 'i' individuals 
#'   from each of the 'r' genetic groups. Default output is an object of class
#'   \code{matrix} (dense), but this format can be changed (e.g., "dgCMatrix"
#'   for a sparse matrix).
#' @author \email{matthewwolak@@gmail.com}
#' @references Fikse, F. 2009. Fuzzy classification of phantom parent groups in
#' an animal model. Genetics, Selection, Evolution. 41:42.
#' 
#' Quaas, R.L. 1988. Additive genetic model with groups and relationships.
#' Journal of Dairy Science. 71:1338-1345.
#' @examples
#' 
#' # Use the pedigree from Quaas 1988 (See `data(Q1988)`)
#' ##########################
#' # Fuzzy classification
#'   ## Fuzzy classification with complete assignment to one group
#'     Q1988fuzz <- Q1988[-c(1:2), c("id", "phantomDam", "phantomSire")]
#'     Qfnull <- matrix(c(1,0,0,1,0, 0,1,1,0,1), nrow = 5, ncol = 2,
#' 	dimnames = list(letters[1:5], c("g1", "g2")))
#'     (Qfuzznull <- ggcontrib(Q1988fuzz, fuzz = Qfnull))
#' 
#'     ## Should be identical to the non-fuzzy classification output
#'     # format (1) from above
#'       (Q <- ggcontrib(Q1988[-c(3:7), c(1,4,5)]))
#'     stopifnot(Qfuzznull == Q)
#' 
#'   ## Fuzzy classification with arbitrary assignments
#'     Qf <- matrix(c(1,0,0.5,0.5,0.5, 0,1,0.5,0.5,0.5), nrow = 5, ncol = 2,
#' 	dimnames = list(letters[1:5], c("g1", "g2")))
#'     (Qfuzz <- ggcontrib(Q1988fuzz, fuzz = Qf))  
#' 
#'   ## Using the pedigree and fuzzy classification in Fikse (2009)
#'     F2009fuzz <- data.frame(id = c(letters[1:7], LETTERS[1:6]),
#' 	dam = c(rep(NA, 7), "a", "c", "e", "A", "C", "D"),
#' 	sire = c(rep(NA, 7), "b", "d", "f", "B", "g", "E"))
#'     Ff <- matrix(c(1,0,1,0,0,0,0.2,
#' 		0,1,0,0.6,0,0.3,0.4,
#' 		0,0,0,0.4,1,0.7,0.4),
#' 		nrow = 7, ncol = 3,
#' 		dimnames = list(letters[1:7], paste0("g", 1:3)))
#'     # Actual Q matrix printed in Fikse (2009)
#'       Fikse2009Q <- matrix(c(0.5,0.5,0,0.5,0.1,0.3,
#' 			0.5,0.3,0.15,0.4,0.275,0.3375, 
#' 			0,0.2,0.85,0.1,0.625,0.3625),
#' 		nrow = 6, ncol = 3,
#' 		dimnames = list(LETTERS[1:6], paste0("g", seq(3))))
#' 
#'     Ffuzz <- ggcontrib(F2009fuzz, fuzz = Ff)
#'       (diffFfuzz <- Ffuzz - Fikse2009Q)
#'       # Encountering some rounding error
#'       stopifnot(length((drop0(diffFfuzz, tol = 1e-12))@x) == 0)
#' 
#' 
#' ##########################
#' # Non-fuzzy classification
#'   # format (1) from above
#'     Q1 <- Q1988[-c(3:7), c(1,4,5)]
#'     (gg1 <- ggcontrib(Q1, ggroups = NULL)) # note the warning message which is typical
#' 
#'   # format (2) from above
#'     Q2 <- Q1988[-c(1:7), 1:3]
#'     # arbitrarily assign individuals genetic groups for unknown parents
#'     ## Means gg2 is NOT comparable to gg1 or gg3!
#'     ggvec.in <- c("g1", "g2", "g1", NA)
#'     (gg2 <- ggcontrib(Q2, ggroups = ggvec.in))
#' 
#'   # format (3) from above
#'     Q3 <- Q1988[-c(1:7), c(1,4,5)]
#'     gg3 <- ggcontrib(Q3, ggroups = c("g1", "g2"))
#' 
#'   stopifnot(gg1 == gg3)
#' 
#' @export
ggcontrib <- function(pedigree, ggroups = NULL, fuzz = NULL, output = "matrix"){
 
  if(!is.null(ggroups) & !is.null(fuzz)){
    stop("Arguments to either 'ggroups' or 'fuzz' can be non-NULL, but not both")
  }
  if(is.null(fuzz)){
  # Without fuzzy classification
    if(is.null(ggroups)){
      ptype <- "A"
      nPed <- numPed(pedigree)
      ggroupsD <- as.character(pedigree[which(nPed[, 2] == -998), 1])
      ggroupsS <- as.character(pedigree[which(nPed[, 3] == -998), 1])
      if(!all(ggroupsD == ggroupsS)){
        stop("Only rows identifying genetic groups should have missing parents.  All individuals in the pedigree must have a genetic group when a parent is unknown")
      }
      nggroups <- nPed[which(nPed[, 2] == -998), 1]

    } else{
        if(length(ggroups) == length(unique(ggroups))){
          ptype <- "D"
          if(any(pedigree[, 2:3] == "0" | pedigree[, 2:3] == "*" | is.na(pedigree[, 2:3]))){
            stop("When specifying the unique genetic groups as a vector in the 'ggroups' argument, all individuals in the pedigree must have a genetic group when a parent is unknown (<NA>, '0' and '*' are considered unknown parents)")
          }
          pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])),
		dam = c(rep(NA, length(ggroups)), as.character(pedigree[, 2])),
		sire = c(rep(NA, length(ggroups)), as.character(pedigree[, 3])))
          nPed <- suppressWarnings(numPed(pedalt))
          nggroups <- nPed[which(nPed[, 2] == -998), 1]
        }

        if(length(ggroups) == dim(pedigree)[1]){
          ptype <- "R"
          nonggped <- pedigree[which(is.na(ggroups)), 2:3]
          if(any(nonggped == "0" | nonggped == "*" | is.na(nonggped))){
            stop("All individuals with missing parents (indicated as '0', '*', or <NA>) must have a genetic group specified")
          }
          if(length(d0 <- which(pedigree[, 2] == 0)) > 0){
            pedigree[d0, 2] <- NA
            warning("Zero in the dam column interpreted as a missing parent")
          }
          if(length(s0 <- which(pedigree[, 3] == 0)) > 0){
            pedigree[s0, 3] <- NA
            warning("Zero in the sire column interpreted as a missing parent")
          }
          if(length(dast <- which(pedigree[, 2] == "*")) > 0){ 
            pedigree[dast, 2] <- NA
          }
          if(length(sast <- which(pedigree[, 3] == "*")) > 0){
            pedigree[sast, 3] <- NA
          }   
          pedalt <- data.frame(id = I(as.character(pedigree[, 1])),
		dam = I(as.character(pedigree[, 2])),
		sire = I(as.character(pedigree[, 3])))
          pedalt[which(is.na(pedigree[, 2])), 2] <- as.character(ggroups[which(is.na(pedigree[, 2]))])
          pedalt[which(is.na(pedigree[, 3])), 3] <- as.character(ggroups[which(is.na(pedigree[, 3]))])
          uggroups <- as.character(unique(ggroups))
          uggroups <- uggroups[!is.na(uggroups)]
          pedalt <- data.frame(id = c(uggroups, as.character(pedalt[, 1])),
		dam = c(rep(NA, length(uggroups)), as.character(pedalt[, 2])),
		sire = c(rep(NA, length(uggroups)), as.character(pedalt[, 3])))
          nPed <- suppressWarnings(numPed(pedalt))
          nggroups <- nPed[which(nPed[, 2] == -998), 1]       
        }
      }
  } else{
  # Below is fuzzy classification
      ptype <- "F"
      nPed <- numPed(pedigree)
      na2 <- nPed[, 2] == -998
      na3 <- nPed[, 3] == -998

      # checks on fuzzy classification matrix and pedigree consistency:
      if(is.null(dimnames(fuzz)[[1]]) | is.null(dimnames(fuzz)[[2]])){
        stop("dimnames(fuzz) must have phantom parent identities matching the pedigree for row names and genetic groups for column names")
      } 
      if(any(colnames(fuzz) %in% pedigree[, 1])){
        stop("genetic groups (column names of `fuzz`) cannot have an entry in the pedigree")
      }
      if(any(rowSums(fuzz) != 1)){
        stop("all rowSums(fuzz) must equal 1")
      }
      if(any(I(na2 + na3) == 1)){
        stop("Observed individuals/phantom parents must have either 2 known parents or 2 missing parent identities (e.g., NA)")
      }
      if(!all(pedigree[na2, 1] %in% rownames(fuzz))){
        stop(paste0(pedigree[na2, 1][which(!pedigree[na2, 1] %in% rownames(fuzz))], ": phantom dams in the pedigree do not have an entry/row name in fuzz"))
      }
      if(!all(pedigree[na3, 1] %in% rownames(fuzz))){
        stop(paste0(pedigree[na3, 1][which(!pedigree[na3, 1] %in% rownames(fuzz))], ": phantom sires in the pedigree do not have an entry/row name in fuzz"))
      }
      if(!all(nPed[match(rownames(fuzz), pedigree[, 1]), 2] == -998)){
        stop(paste0(pedigree[match(rownames(fuzz), pedigree[, 1]), 1][nPed[match(rownames(fuzz), pedigree[, 1]), 2] != -998], ": phantom parents in fuzz do not have unknown/missing dams in the pedigree"))
      }
      if(!all(nPed[match(rownames(fuzz), pedigree[, 1]), 3] == -998)){
        stop(paste0(pedigree[match(rownames(fuzz), pedigree[, 1]), 1][nPed[match(rownames(fuzz), pedigree[, 1]), 3] != -998], ": phantom parents in fuzz do not have unknown/missing sires in the pedigree"))
      }
      # end checks

      Qb <- as(fuzz, "dgCMatrix")
      pp <- which(na2)            #<-- phantom parents
      oid <- which(!na2)          #<-- observed individuals
    }
  # End strictly fuzzy classification section  

    N <- dim(nPed)[1]
    if(is.null(fuzz)){
      dnmiss <- which(nPed[, 2] != -998)
      snmiss <- which(nPed[, 3] != -998)
      maxcnt <- (length(dnmiss) + length(snmiss) + N)
    } else maxcnt <- 2*length(oid) + N
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

  if(ptype == "A" | ptype == "F"){
    Tinv <- t(sparseMatrix(i = Cout[[3]][1:Cout[[6]]], p = Cout[[4]], x = Cout[[5]][1:Cout[[6]]],
	dims = c(N, N),
	dimnames = list(as.character(pedigree[, 1]), as.character(pedigree[, 1])),
	symmetric = FALSE, index1 = FALSE))
  } else{
      Tinv <- t(sparseMatrix(i = Cout[[3]][1:Cout[[6]]], p = Cout[[4]], x = Cout[[5]][1:Cout[[6]]],
	dims = c(N, N),
	dimnames = list(as.character(pedalt[, 1]), as.character(pedalt[, 1])),
	symmetric = FALSE, index1 = FALSE))
    }


  if(is.null(fuzz)){
    T <- as(solve(Tinv), "dgCMatrix")
    T@Dimnames <- Tinv@Dimnames
   return(as(T[-c(nggroups), nggroups], output))
  } else{
      Q <- as(solve(Tinv[oid, oid]), "dgCMatrix") %*% (-1*Tinv[oid, pp]) %*% Qb
      Q@Dimnames <- list(Tinv@Dimnames[[1]][oid], Qb@Dimnames[[2]])
     return(as(Q, output))
    }
}
