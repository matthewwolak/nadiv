# Create Astar by matrix multiplication of A^-1 and Q submatrices


#' Creates the inverse additive genetic relationship matrix with genetic groups
#' 
#' This returns the inverse of the additive genetic relationship matrix with
#' genetic groups (A*). The matrix is set up through matrix multiplication of
#' two sub-matrices instead of directly (as \code{\link{makeAinv}} does).
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' The function implements the matrix multiplication, using sub-matrices
#' \code{Q} and \code{A^-1}, as detailed in Quaas (1988, pp. 1342-1343).
#' 
#' Genetic groups can be incorporated into the A-inverse by providing a value
#' to the \code{ggroups} argument. The value supplied to \code{ggroups} can
#' either be (1) a single integer indicating the number of unique genetic
#' groups or (2) a character vector containing the name for each genetic group.
#' These are referred to as pedigree types "A" and "D", respectively, and
#' further details follow below.  (Type="A") the pedigree contains unique IDs
#' for the 'g' genetic groups in the first 'g' lines of the pedigree. The dam
#' and sire of the genetic group rows should contain missing values (e.g., NA,
#' "0", or "*"). All individuals in the pedigree should then have one of the
#' 'g' genetic groups instead of an unknown parent.  (Type="D") the pedigree
#' contains only individuals in the ID column (no genetic groups have an ID)
#' and there should be no missing values for any dams or sires. Instead,
#' individuals for whom the dam and/or sire is unknown should have one of the
#' genetic groups identified in the vector supplied to \code{ggroups} as the
#' dam or sire.
#' 
#' Fuzzy classification of genetic groups is implemented when \code{fuzz} is
#' non-NULL.
#' 
#' The argument to \code{gOnTop} specifies if the elements in the A-inverse
#' should come at the beginning (\code{gOnTop = TRUE}) or end (\code{gOnTop =
#' FALSE}) of the matrix. Depending on how the software implementing an animal
#' model solves the mixed model equations, the equations for the genetic groups
#' (and thus the elements in the augmented A-inverse) should be the first or
#' last set of equations.
#' 
#' See function \code{\link{makeAinv}} for directly obtaining the inverse of
#' the additive genetic relationship matrix with genetic groups.
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire
#' @param ggroups Either a vector with the unique name of each genetic group,
#'   or a numeric indicating the number of unique genetic groups. See Details 
#'   for different ways to specify. Note, cannot be NULL.
#' @param fuzz A matrix containing the fuzzy classification of individuals into
#'   genetic groups.
#' @param gOnTop A logical indicating if the A-inverse should be constructed
#'   with the \sQuote{g} genetic groups located in the first \sQuote{g} rows 
#'   and columns if TRUE, else the \sQuote{g} genetic groups are located in the 
#'   last \sQuote{g} rows and columns of A-inverse.
#'
#' @return Returns A*, or the inverse of the numerator relationship with
#'   groups, in sparse matrix form.
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{makeAinv}}, \code{\link{ggcontrib}}
#' @references Quaas, R.L. 1988. Additive genetic model with groups and
#' relationships. Journal of Dairy Science. 71:1338-1345.
#' @examples
#' 
#'  # Using the Q1988 dataset in nadiv
#'  ## assign a null fuzzy classification matrix
#'  QfuzzNull <- matrix(c(1,0,0,1,0, 0,1,1,0,1), nrow = 5, ncol = 2,
#' 	dimnames = list(letters[1:5], c("g1", "g2")))
#' 
#'  # Type A
#'  ## no fuzzy classification
#'   Astar_A <- makeAstarMult(Q1988[-c(3:7), c(1,4,5)], ggroups = 2)
#'  ## with fuzzy classification
#'   Astar_Afuzzy <- makeAstarMult(Q1988[, c(1, 6, 7)],
#' 	ggroups = 2, fuzz = QfuzzNull)
#' 
#'  # Type D
#'  ## no fuzzy classification
#'   Astar_D <- makeAstarMult(Q1988[-c(1:7), c(1, 4, 5)], ggroups = c("g1", "g2"))
#'  ## with fuzzy classification
#'   Astar_Dfuzzy <- makeAstarMult(Q1988[-c(1:2), c(1, 6, 7)],
#' 	ggroups = c("g1", "g2"), fuzz = QfuzzNull)
#' 
#' 
#'  # Obtain the matrix directly 
#'  ## no fuzzy classification
#'  Astar_direct <- makeAinv(Q1988[-c(3:7), c(1,4,5)], ggroups = 2)$Ainv
#'  stopifnot(length(drop0(round(Astar_direct
#' 	- (Astar_A - Astar_Afuzzy)
#' 	- (Astar_D - Astar_Dfuzzy)
#' 	- Astar_direct, 10))@x) == 0)
#' 
#'  ## with fuzzy classification
#'  Astar_directF <- makeAinv(Q1988[-c(1:2), c(1, 6, 7)], fuzz = QfuzzNull)$Ainv
#'  stopifnot(length(drop0(round(Astar_directF
#' 	- (Astar_A - Astar_Afuzzy)
#' 	- (Astar_D - Astar_Dfuzzy)
#' 	- Astar_direct, 10))@x) == 0)
#'  
#' 
#' @export
makeAstarMult <- function(pedigree, ggroups, fuzz = NULL, gOnTop = FALSE){
  if(is.null(fuzz)){
    if(inherits(ggroups, "numeric") && length(ggroups) == 1){
      Q <- ggcontrib(pedigree)
      pedAlt <- pedigree[-seq(ggroups), ]
      ggroups <- pedigree[seq(ggroups), 1]
    } else{
        Q <- ggcontrib(pedigree, ggroups)
        pedAlt <- pedigree
      }
    pedAlt[which(pedAlt[, 2] %in% ggroups), 2] <- NA
    pedAlt[which(pedAlt[, 3] %in% ggroups), 3] <- NA
  } else{
      if(inherits(ggroups, "numeric") && length(ggroups) == 1){
        pedAlt <- pedigree[-seq(ggroups), ]
        Q <- ggcontrib(pedAlt, fuzz = fuzz)
        pedAlt <- pedAlt[-match(rownames(fuzz), pedAlt[, 1]), ]
      } else{
          Q <- ggcontrib(pedigree, fuzz = fuzz)
	  pedAlt <- pedigree[-match(rownames(fuzz), pedigree[, 1]), ]
        }
      pedAlt[which(pedAlt[, 2] %in% rownames(fuzz)), 2] <- NA
      pedAlt[which(pedAlt[, 3] %in% rownames(fuzz)), 3] <- NA
    }

  Ainv <- makeAinv(pedAlt)$Ainv
  # U matrix defined pp.1342-1343, Quaas 1988
  if(gOnTop){
    U <- cbind(-Q, Diagonal(x = 1, n = nrow(Ainv)))
  } else{
      U <- cbind(Diagonal(x = 1, n = nrow(Ainv)), -Q)
    }
  Astar <- drop0(zapsmall(crossprod(U, Ainv) %*% U, 12))
    Astar@Dimnames[[1L]] <- if(gOnTop) c(colnames(Q), rownames(Ainv)) else c(rownames(Ainv), colnames(Q))
 Astar
}

