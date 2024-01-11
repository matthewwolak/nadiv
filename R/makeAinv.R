# Generic


#' Creates the inverse additive genetic relationship matrix
#' 
#' This returns the inverse of the numerator relationship matrix (inverse
#' additive genetic relatedness matrix). It can also be used to obtain
#' coefficients of inbreeding for the pedigreed population.
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' The functions implement an adaptation of the Meuwissen and Luo (1992)
#' algorithm (particularly, following the description of the algorithm in
#' Mrode 2005) with some code borrowed from the \code{inverseA} function by
#' Jarrod Hadfield in the \code{MCMCglmm} package. Further, providing a
#' non-NULL argument to \code{ggroups} incorporates the Quaas (1988) algorithm
#' for directly obtaining the augmented A-inverse matrix for genetic groups
#' into Meuwissen and Luo's (1992) algorithm, thereby, considering inbreeding
#' during the construction of the A-inverse. Further calculations needed for
#' the algorithm to incorporate inbreeding and genetic groups follow the theory
#' presented in VanRaden (1992). Alternatively, group-specific inverse
#' relatedness matrices can be formed with \code{makeGGAinv}, see below.
#' 
#' At the moment, providing the inbreeding level of individuals or the base
#' population has not been implemented. However, this argument is a placeholder
#' for now.
#' 
#' Genetic groups can be incorporated into a single A-inverse by providing a value
#' to the \code{ggroups} argument in \code{makeAinv}. The value supplied to
#' \code{ggroups} can either be (1) a single integer indicating the number of
#' unique genetic groups or (2) a character vector containing the name for each
#' genetic group. These are referred to as pedigree types "A" and "D",
#' respectively, and further details follow below.
#'
#' (Type="A") the pedigree contains unique IDs for the 'g' genetic groups in the
#' first 'g' lines of the pedigree. The dam and sire of the genetic group rows
#' should contain missing values (e.g., NA, "0", or "*"). All individuals in the
#' pedigree should then have one of the `g' genetic groups instead of an unknown
#' parent.
#' (Type="D") the pedigree contains only individuals in the ID column (no
#' genetic groups have an ID) and there should be no missing values for any dams
#' or sires. Instead, individuals for whom the dam and/or sire is unknown should
#' have one of the genetic groups identified in the vector supplied to
#' \code{ggroups} as the dam or sire.
#' 
#' \sQuote{Fuzzy classification} of genetic groups (Fikse 2009) can be
#' implemented if a \sQuote{matrix} (of class \code{matrix} or \code{Matrix})
#' is supplied to the \code{fuzzy} argument. The fuzzy classification matrix
#' must have row names matching all of the phantom parents in the pedigree and
#' the column names must be present and specify the genetic groups. The fuzzy
#' classification matrix essentially contains probability of group membership
#' for each phantom parent. Therefore, each row should sum to 1. The pedigree
#' must have an identity in a unique row for every phantom parent and cannot
#' have genetic groups as either identities (in the first column) or as dam or
#' sire (second and third columns). Further, if fuzzy classification is
#' desired, the function must specify \code{ggroups = NULL}.
#' 
#' When genetic groups (including the case of fuzzy classification of genetic
#' groups) are included in the A-inverse matrix, the argument to \code{gOnTop}
#' specifies if the genetic group elements in the A-inverse should occupy the
#' top-left (\code{gOnTop = TRUE}) or bottom-right (\code{gOnTop = FALSE}) of
#' the matrix. Depending on how the software implementing an animal model
#' solves the mixed model equations, the equations for the genetic groups (and
#' thus the elements in the augmented A-inverse) should be the first or last
#' set of equations.
#' 
#' @aliases makeAinv makeAinv.default makeAinv.fuzzy makeGGAinv
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param f A numeric indicating the level of inbreeding. See Details
#' @param ggroups Either a vector with the unique name of each genetic group,
#'   or a numeric indicating the number of unique genetic groups. See Details 
#'   for different ways to specify. Note, if NULL then the regular A-inverse
#'   will be constructed. Also, must be NULL if fuzz is non-NULL.
#' @param fuzz A matrix containing the fuzzy classification of phantom parents
#'   into genetic groups. See Details.
#' @param gOnTop A logical indicating if (when including genetic groups) the
#'   A-inverse should be constructed with the \sQuote{g} genetic groups located 
#'   in the first \sQuote{g} rows and columns if \code{TRUE}, else the
#'   \sQuote{g} genetic groups are located in the last \sQuote{g} rows and 
#'   columns of A-inverse
#' @param det Logical, indicating if the (log) determinant of the A matrix
#'   should be returned
#' @param \dots Arguments to be passed to methods
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{Ainv }{the inverse of the additive genetic relationship matrix
#'       in sparse matrix form}
#'     \item{listAinv }{the three column list of the non-zero elements for the 
#'       inverse of the additive genetic relationship matrix with attributes
#'       \code{rowNames} and \code{geneticGroups}. \code{attr(*, "rowNames")}
#'       links the integer for rows/columns to the ID column from the pedigree. 
#'       \code{attr(*, "geneticGroups")} is a two element vector with the first 
#'       integer indicating how many genetic groups are included in the 
#'       pedigree. This last attribute is necessary for some software programs 
#'       to correctly specify the residual degrees of freedom when calculating 
#'       the log-likelihood in a model that implicitly fits fixed genetic group 
#'       effects.}
#'     \item{f }{the individual coefficients of inbreeding for each individual 
#'       in the pedigree (matches the order of the first/ID column of the
#'       pedigree). If the pedigree contains \sQuote{g} genetic groups in the 
#'       first \sQuote{g} rows, then the first \sQuote{g} elements of \code{f} 
#'       are assigned 0. If the pedigree contains \sQuote{p} phantom parents in 
#'       the first \sQuote{p} rows, then the first \sQuote{p} elements of
#'       \code{f} are assigned 0.}
#'     \item{logDet }{the log determinant of the A matrix}
#'     \item{dii }{the (non-zero) elements of the diagonal D matrix of the A=TDT'
#'       decomposition. Contains the variance of Mendelian sampling. Matches
#'       the order of the first/ID column of the pedigree. If the pedigree
#'       contains \sQuote{g} genetic groups in the first \sQuote{g} rows, then
#'       the first \sQuote{g} elements of \code{f} are assigned 0. If the
#'       pedigree contains \sQuote{p} phantom parents in the first \sQuote{p}
#'       rows, then the first \sQuote{p} elements of \code{f} are assigned 0.} 
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeAstarMult}}, \code{\link{makeA}}
#' @references Fikse, F. 2009. Fuzzy classification of phantom parent groups in
#' an animal model. Genetics Selection Evolution 41:42.
#' 
#' Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding coefficients in large
#' populations. Genetics, Selection, Evolution. 24:305-313.
#' 
#' Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
#' Values, 2nd ed.  Cambridge, MA: CABI Publishing.
#' 
#' Quaas, R.L. 1988. Additive genetic model with groups and relationships.
#' Journal of Dairy Science. 71:1338-1345.
#' 
#' VanRaden, P.M. 1992. Accounting for inbreeding and crossbreeding in genetic
#' evaluation of large populations. Journal of Dairy Science. 75:3136-3144.
#' @examples
#' 
#'  ##  Without genetic groups  ##
#'  makeAinv(Mrode2)
#'  
#'  ##  With genetic groups  ##
#'   ## Type A
#'    typeAped <- Q1988[-c(3:7), c("id", "damGG", "sireGG")]
#'    AstarA <- makeAinv(typeAped, ggroups = 2, gOnTop = FALSE)$Ainv
#'   ## Type D
#'    typeDped <- Q1988[-c(1:7), c("id", "damGG", "sireGG")]
#'    AstarD <- makeAinv(typeDped, ggroups = c("g1", "g2"), gOnTop = FALSE)$Ainv
#'   stopifnot(identical(AstarA, AstarD))
#'   
#'   # Show that the augmented A-inverse with genetic groups
#'   # contains the normal A-inverse (i.e., without genetic groups)
#'    ## Augmented A-inverse with genetic groups
#'     ggAinv <- makeAinv(Mrode3[-c(1,2), c("calf", "damGG", "sireGG")],
#' 	ggroups = c("g1", "g2"), gOnTop = FALSE)$Ainv
#'     noggAinv <- makeAinv(Mrode3[-c(1,2), c("calf", "dam", "sire")],
#' 	ggroups = NULL)$Ainv
#'     # First 8 rows & columns of ggAinv are same as A-inverse without 
#'     ## genetic groups
#'     ggAinv[1:8, 1:8]
#'     noggAinv
#'    stopifnot(all.equal(ggAinv[1:8, 1:8], structure(noggAinv, geneticGroups = NULL)))
#'    
#'  ##  With fuzzy classification of genetic groups  ##
#'   ## example in Fikse (2009)
#'   Fped <- F2009[-c(1:3), c("id", "phantomDam", "phantomSire")]
#'     Fped$id <- factor(Fped$id, levels = as.character(unique(Fped$id)))
#'   Ffuzz <- as.matrix(F2009[4:10, c("g1", "g2", "g3")])
#'     dimnames(Ffuzz)[[1]] <- as.character(F2009[4:10, 1])
#'   AstarF <- makeAinv(Fped, fuzz = Ffuzz, gOnTop = FALSE)$Ainv
#' 
#'   ## Show that A-inverse with fuzzy classification of genetic groups
#'   ### can be the same as genetic group A-inverse without fuzzy classification
#'   ### Create a 'null' fuzzy classification matrix for Q1988 pedigree
#'   QfuzzNull <- matrix(c(1,0,0,1,0, 0,1,1,0,1), nrow = 5, ncol = 2,
#' 	dimnames = list(letters[1:5], c("g1", "g2")))
#'   typeFped <- Q1988[-c(1:2), c("id", "phantomDam", "phantomSire")]
#'   AstarNullFuzzy <- makeAinv(typeFped, fuzz = QfuzzNull, gOnTop = FALSE)$Ainv
#'   # Same as above using either pedigree type 'A' or 'D'
#'   stopifnot(identical(AstarNullFuzzy, AstarA),
#' 	identical(AstarNullFuzzy, AstarD))
#' 
#'  ##  With genetic groups  ##
#'   ## Type A
#'    typeAped <- Q1988[-c(3:7), c("id", "damGG", "sireGG")]
#'    (AinvOutA <- makeGGAinv(typeAped, ggroups = 2)$Ainv)
#'   ## Type D
#'    typeDped <- Q1988[-c(1:7), c("id", "damGG", "sireGG")]
#'    (AinvOutD <- makeGGAinv(typeDped, ggroups = c("g1", "g2"))$Ainv)
#'   stopifnot(identical(AinvOutA, AinvOutD))
#' 
#' @export
makeAinv <- function(pedigree, f = NULL, ggroups = NULL, fuzz = NULL, gOnTop = FALSE, det = TRUE, ...){
  if(is(fuzz, "matrix") | is(fuzz, "Matrix")) class(fuzz) <- "fuzzy"
  UseMethod("makeAinv", fuzz)
}


###############################################################################
###############################################################################
# 	Methods: `makeAinv.default()` and `makeAinv.fuzzy()`
###############################################################################
###############################################################################

#' @method makeAinv default
#' @rdname makeAinv
#' @export
makeAinv.default <- function(pedigree, f = NULL,
  ggroups = NULL, fuzz = NULL, gOnTop = FALSE, det = TRUE, ...){
  clnms <- names(match.call())
  if("groups" %in% clnms){
    stop("Found argument 'groups', only use correct argument name 'ggroups'")
  }    
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
  f <- c(rep(-1, nggroups), rep(0, eN), -1) #TODO allow user to specify some f
  Cout <- .C("ainvml", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), 				#dam
	    as.integer(nPed[, 3] - 1),  			#sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N),   					#n
            as.integer(nggroups),   				#g
            as.double(rep(0, length(Ainv@i))),  		#xA
	    as.integer(Ainv@i), 				#iA
	    as.integer(Ainv@p), 				#pA
	    as.integer(length(Ainv@i))) 			#nzmaxA
  Ainv <- as(Ainv, "dMatrix")
  Ainv@x <- Cout[[7]]
  fsOrd <- as(as.integer(renPed), "pMatrix")
  Ainv <- crossprod(fsOrd, Ainv) %*% fsOrd
  if(ptype == "D"){
      Ainv@Dimnames <- list(as.character(pedalt[, 1]), NULL)
      f <- Cout[[3]][t(fsOrd)@perm][-seq(nggroups)]
      dii <- Cout[[4]][t(fsOrd)@perm][-seq(nggroups)]
  } else {
      Ainv@Dimnames <- list(as.character(pedigree[, 1]), NULL)
      f <- c(rep(0, nggroups), Cout[[3]][t(fsOrd)@perm][(nggroups+1):(nggroups + eN)])
      dii <- c(rep(0, nggroups), Cout[[4]][t(fsOrd)@perm][(nggroups+1):(nggroups + eN)])
    }
  if(!is.null(ggroups) && !gOnTop){ 
     permute <- as(as.integer(c(seq(eN+1, N, 1), seq(eN))), "pMatrix")
       dmnms <- Ainv@Dimnames[[1L]][invPerm(permute@perm)]
     Ainv <- crossprod(permute, Ainv) %*% permute
       Ainv@Dimnames[[1L]] <- dmnms
  }
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = structure(Ainv, geneticGroups = c(nggroups, 0)),
	listAinv = structure(sm2list(Ainv, rownames = rownames(Ainv),
	                       colnames = c("row", "column", "Ainv")),
	                     geneticGroups = c(nggroups, 0)),
	f = f,
	logDet = logDet,
	dii = dii))
}





###############################################################################
###############################################################################

#' @method makeAinv fuzzy
#' @rdname makeAinv
#' @export
makeAinv.fuzzy <- function(pedigree, f = NULL, ggroups = NULL, fuzz, gOnTop = FALSE, det = TRUE, ...){

  if(!is.null(ggroups)){
    stop("when 'fuzz' is non-NULL, 'ggroups' should not have any arguments (i.e., 'ggroups==NULL")
  }

  naPed2 <- which(pedigree[, 2] == "0")
  naPed3 <- which(pedigree[, 3] == "0")
  naPed2 <- union(naPed2, which(pedigree[, 2] == "*"))
  naPed3 <- union(naPed3, which(pedigree[, 3] == "*"))
  naPed2 <- union(naPed2, which(is.na(pedigree[, 2])))
  naPed3 <- union(naPed3, which(is.na(pedigree[, 3])))

  # checks on fuzzy classification matrix and pedigree consistency:
  ## 'fuzz' is a type of matrix
  if(!is(fuzz, "matrix") && !is(fuzz, "Matrix")){
    cat("'fuzz' of class", class(fuzz), "\n")
    stop("class of 'fuzz' must be either 'matrix' or 'Matrix'")
  }
  ## rows of 'fuzz' add up to 1
  if(any(rowSums(fuzz) != 1)){
    cat("rows:", which(rowSums(fuzz) != 1), "\ndo not equal 1\n")
    stop("all rowSums(fuzz) must equal 1\n(check for possible rounding errors, e.g., 3*0.33333 != 1)")
  }
  ## fuzz has dimnames
  if(is.null(dimnames(fuzz)[[1]]) | is.null(dimnames(fuzz)[[2]])){
    stop("'fuzz' must have row and column names")
  } 
  ## pedigree does not have genetic groups
  if(any(colnames(fuzz) %in% pedigree[, 1])){
    cat("colnames:", which(colnames(fuzz) %in% pedigree[, 1]), "\nare in 'pedigree'\n")
    stop("colnames of 'fuzz' (genetic groups) must NOT be identities in the first column of 'pedigree'")
  }
  ## pedigree has phantom parents in 'fuzz'
  if(!all(rownames(fuzz) %in% pedigree[, 1])){
    cat("rownames:", which(!rownames(fuzz) %in% pedigree[, 1]), "\nnot in 'pedigree'\n")
    stop("rownames of 'fuzz' (phantom parents) must all be identities in the first column of 'pedigree'. See the `prepPed()` function to help prepare the pedigree")
  }
  ## individuals can only have both parents missing in 'pedigree' or none
  if(length(naPed2) != length(naPed3) | any(!naPed2 %in% naPed3)){
    stop("Individuals must have either two or zero missing parents in 'pedigree'")
  }   
  ## IDs with missing parents (if passed above check, naPed2==naPed3) in 'pedigree' are phantom parents in 'fuzz'
  if(!all(pedigree[naPed2, 1] %in% rownames(fuzz))){
    cat("IDs for 'pedigree' rows:", naPed2[which(!pedigree[naPed2, 1] %in% rownames(fuzz))], "\nare not rownames in 'fuzz'\n")
    stop("Individuals with missing parents (phantom individuals) must have a rowname in 'fuzz'")
  }


  # order of genetic groups in pedalt/A^-1 is same as column order of 'fuzz'
  ggroups <- colnames(fuzz)
  nggroups <- length(ggroups) 			# No. genetic groups
  p <- nrow(fuzz) 				# No. phantom parents
  eN <- nrow(pedigree) - p 			# No. observed IDs
  N <- nggroups + eN 				# No. GGs + observed IDs

  # order of phantom parents in pedalt is same as row order in 'fuzz'
  phantomPars <- seq(p) + nggroups
  groupRows <- seq(nggroups)
  # order pedigree: generations, order phantom parents in fuzz, dam, & sire 
  ## genetic groups first 
  renPed <- c(groupRows, nggroups + order(genAssign(pedigree), match(pedigree[, 1], rownames(fuzz), nomatch = p+1), pedigree[, 2], pedigree[, 3]))
  pedalt <- data.frame(id = c(ggroups, as.character(pedigree[, 1])), dam = c(rep(NA, nggroups), as.character(pedigree[, 2])), sire = c(rep(NA, nggroups), as.character(pedigree[, 3])))[renPed, ]
  nPed <- numPed(pedalt)
  phantomless <- cbind(numPed(nPed[-phantomPars, ], check = FALSE), nPed[-phantomPars, -1]) 

  if(!all(which(phantomless[, 4] == -998) == groupRows)){
    stop("Something wicked happened with the dams in phantomless numeric pedigree:\n  Contact package maintainer: <matthewwolak@gmail.com>\n  or raise an issue on: <https://github.com/matthewwolak/nadiv/issues>")
  }
  if(!all(which(phantomless[, 5] == -998) == groupRows)){
    stop("Something wicked happened with the sires in phantomless numeric pedigree:\n  Contact package maintainer: <matthewwolak@gmail.com>\n  or raise an issue on: <https://github.com/matthewwolak/nadiv/issues>")
  }

 
  groupFuzz <- Diagonal(x = 1, n = nggroups)
  groupFuzz@Dimnames <- list(as.character(ggroups), as.character(ggroups))
  fuzzmat <- as(rbind(groupFuzz, fuzz), "CsparseMatrix")
  # predict non-zero elements of Astar
  ## make H from Quaas 1988:
  ## H = [-Pb Qb : Tinv]
  ## Astar = H' %*% D^-1 %*% H
  ### sTinv: n x n observed IDs (i.e., no GGs or phantom parent IDs)
  dnmiss <- which(phantomless[, 2] != -998)
  snmiss <- which(phantomless[, 3] != -998)
  Tinv.row <- c(c(phantomless[, 1][dnmiss], phantomless[, 1][snmiss]) - nggroups, 1:eN)
  Tinv.col <- c(c(phantomless[, 2][dnmiss], phantomless[, 3][snmiss]) - nggroups, 1:eN)
  el.order <- order(Tinv.col + Tinv.row/(eN + 1), decreasing = FALSE)
  sTinv <- sparseMatrix(i = as.integer(Tinv.row[el.order] - 1),
	  p = as.integer(c(match(1:eN, Tinv.col[el.order]), length(el.order) + 1) - 1),
	  index1 = FALSE, dims = c(eN, eN), symmetric = FALSE,
	  dimnames = list(as.character(phantomless[-groupRows, 1]), NULL))
  ### Pb: n x p version of sTinv
  pdnmiss <- which(phantomless[, 4] %in% phantomPars)
  psnmiss <- which(phantomless[, 5] %in% phantomPars)
  Pb.row <- c(phantomless[, 1][pdnmiss], phantomless[, 1][psnmiss]) - nggroups
  Pb.col <- c(phantomless[, 4][pdnmiss], phantomless[, 5][psnmiss]) - nggroups
  el.order <- order(Pb.col + Pb.row/(p + 1), decreasing = FALSE)
  sPb <- sparseMatrix(i = as.integer(Pb.row[el.order] - 1),
	  p = as.integer(c(match(1:p, Pb.col[el.order]), length(el.order) + 1) - 1),
	  index1 = FALSE, dims = c(eN, p), symmetric = FALSE,
	  dimnames = list(NULL, as.character(pedalt[phantomPars, 1])))
  ### Qb is the fuzzy classification matrix ('fuzz')
  Qb <- fuzzmat[-groupRows, ][match(rownames(fuzzmat)[-groupRows], colnames(sPb)), ]
  sQb <- sparseMatrix(i = Qb@i,
	  p = Qb@p,
	  index1 = FALSE, dims = Qb@Dim, symmetric = FALSE,
	  dimnames = Qb@Dimnames)
  ## sH = [-(sPb %*% sQb) : sTinv]
  sH <- cbind((sPb %*% sQb), sTinv)
  Ainv <- t(crossprod(sH))  # transpose stores lower triangle

  phantomless[phantomless == -998] <- N + 1
  # for now, phantom parents cannot be inbred (just like genetic groups)
  f <- c(rep(-1, nggroups), rep(0, eN), -1) #TODO allow user to specify some f
  Cout <- .C("ainvfuzz", PACKAGE = "nadiv",
	    as.integer(phantomless[, 2] - 1), 			#dam
	    as.integer(phantomless[, 3] - 1),  			#sire
	    as.integer(phantomless[, 4] - 1), 			#phantom dam
	    as.integer(phantomless[, 5] - 1),  			#phantom sire
	    as.double(f),					#f
            as.double(rep(0, N)),  				#dii
            as.integer(N),   					#n
            as.integer(nggroups),   				#g
            as.double(fuzzmat@x),                               #xF
            as.integer(fuzzmat@i),				#iF
            as.integer(fuzzmat@p),				#pF
            as.double(rep(0, length(Ainv@i))),  		#xA
	    as.integer(Ainv@i), 				#iA
	    as.integer(Ainv@p)) 				#pA

  Ainv <- as(Ainv, "dMatrix")
  Ainv@x <- Cout[[12]]
  fsOrd1 <- as(as(as(as.integer(renPed), "pMatrix")[, -c(naPed2 + nggroups)], "nMatrix"), "CsparseMatrix")
  fsOrd <- as(as(fsOrd1 %*% matrix(seq(N), nrow = N), "CsparseMatrix")@x, "pMatrix")
  Ainv <- crossprod(fsOrd, Ainv) %*% fsOrd
  Ainv@Dimnames <- list(as.character(pedalt[crossprod(fsOrd1, matrix(seq(N+p), ncol = 1))@x, 1]), NULL)
  f <- (fsOrd1 %*% Cout[[5]][-c(N+1)])@x[-groupRows]
  dii <- (fsOrd1 %*% Cout[[6]][-c(N+1)])@x[-groupRows]
  if(!gOnTop){ 
    permute <- as(as.integer(c(seq(eN+1, N, 1), seq(eN))), "pMatrix")
      dmnms <- Ainv@Dimnames[[1L]][invPerm(permute@perm)]
    Ainv <- crossprod(permute, Ainv) %*% permute
      Ainv@Dimnames[[1L]] <- dmnms
  }
  if(det) logDet <- -1*determinant(Ainv, logarithm = TRUE)$modulus[1] else logDet <- NULL

 return(list(Ainv = structure(Ainv, geneticGroups = c(nggroups, 0)),
	listAinv = structure(sm2list(Ainv, rownames = rownames(Ainv),
	                               colnames = c("row", "column", "Ainv")),
	                     geneticGroups = c(nggroups, 0)),
	f = f,
	logDet = logDet,
	dii = dii))
}



