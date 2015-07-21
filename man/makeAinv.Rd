\name{makeAinv}
\alias{makeAinv}
\title{Creates the inverse additive genetic relationship matrix
}
\description{This returns the inverse of the additive relationship matrix. It can also be used to obtain coefficients of inbreeding for the pedigreed population.
}
\usage{
makeAinv(pedigree, f = NULL, ggroups = NULL, fuzz = NULL, keepPhantoms = TRUE,
      gOnTop = FALSE, det = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
 \item{f }{
A numeric indicating the level of inbreeding. See Details
}
 \item{ggroups }{
Either a vector with the unique name of each genetic group, or a numeric indicating the number of unique genetic groups. See Details for different ways to specify. Note, if NULL then the regular A-inverse will be constructed.
 }
  \item{fuzz }{
A matrix containing the fuzzy classification of individuals into genetic groups. NOT yet implemented.
  }
  \item{keepPhantoms }{
A logical value indicating if phantom parents that are added to the pedigree should be kept and returned in the A-inverse
  }
  \item{gOnTop }{
A logical indicating if the A-inverse should be constructed with the `g` genetic groups located in the first `g` rows and columns if TRUE, else the `g` genetic groups are located in the last `g` rows and columns of A-inverse
  }
  \item{det }{
Logical, indicating if the (log) determinant of the A matrix should be returned
}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

The function implements an adaptation of the Meuwissen and Luo (1992) algorithm with some code borrowed from the \code{inverseA} function by Jarrod Hadfield in the \code{MCMCglmm} package. Further, providing a non-NULL argument to \code{ggroups} incorporates the Quaas (1988) algorithm for directly obtaining the augmented A-inverse matrix for genetic groups into Meuwissen and Luo's (1992) algorithm, thereby, considering inbreeding during the construction of the A-inverse.

At the moment, providing the inbreeding level of individuals or the base population has not been implemented. However, this argument has been introduced as a placeholder for future versions to be released (soon hopefully!).

Genetic groups can be incorporated into the A-inverse by providing a value to the \code{ggroups} argument. The value supplied to \code{ggroups} can either be (1) a single integer indicating the number of unique genetic groups or (2) a character vector containing the name for each genetic group. These are referred to as pedigree types "A" and "D", respectively, and further details follow below.
  (Type="A") the pedigree contains unique IDs for the 'g' genetic groups in the first 'g' lines of the pedigree. The dam and sire of the genetic group rows should contain missing values (e.g., NA, "0", or "*"). All individuals in the pedigree should then have one of the 'g' genetic groups instead of an unknown parent.
  (Type="D") the pedigree contains only individuals in the ID column (no genetic groups have an ID) and there should be no missing values for any dams or sires. Instead, individuals for whom the dam and/or sire is unknown should have one of the genetic groups identified in the vector supplied to \code{ggroups} as the dam or sire.
  
At the moment, fuzzy classification of genetic groups has not been implemented.

The argument to \code{gOnTop} specifies if the elements in the A-inverse should come at the beggining (\code{gOnTop = TRUE}) or end (\code{gOnTop = FALSE}) of the matrix. Depending on how the software implementing an animal model solves the mixed model equations, the equations for the genetic groups (and thus the elements in the augmented A-inverse) should be the first or last set of equations.
}

\value{
  \item{Ainv }{the inverse of the additive genetic relationship matrix in sparse matrix form}
  \item{listAinv }{the three column form of the non-zero elements for the inverse of the additive genetic relationship matrix}
  \item{f }{the individual coefficients of inbreeding for each individual in the pedigree (matches the order of the first/ID column of the pedigree). If the pedigree contains 'g' genetic groups in the first 'g' rows, then the first 'g' elements of \code{f} are assigned 0}
  \item{logDet }{the log determinant of the A matrix}
}

\references{
Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.

Quaas, R.L. 1988. Additive genetic model with groups and relationships. Journal of Dairy Science. 71:1338-1345.
}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{makeA}}}

\examples{
 ##  Without genetic groups  ##
 makeAinv(Mrode2)
 
 ##  With genetic groups  ##
  ## Type A
   typeAped <- Q1988[-c(3:7), c("id", "damGG", "sireGG")]
   AstarA <- makeAinv(typeAped, ggroups = 2, gOnTop = FALSE)$Ainv
  ## Type D
   typeDped <- Q1988[-c(1:7), c("id", "damGG", "sireGG")]
   AstarD <- makeAinv(typeDped, ggroups = c("g1", "g2"), gOnTop = FALSE)$Ainv
  stopifnot(identical(AstarA, AstarD))
  
  # Show that the augmented A-inverse with genetic groups
  # contains the normal A-inverse (i.e., without genetic groups)
   ## Augmented A-inverse with genetic groups
    ggAinv <- makeAinv(Mrode3[-c(1,2), c("calf", "damGG", "sireGG")],
	ggroups = c("g1", "g2"), gOnTop = FALSE)$Ainv
    noggAinv <- makeAinv(Mrode3[-c(1,2), c("calf", "dam", "sire")],
	ggroups = NULL)$Ainv
    # First 8 rows & columns of ggAinv are same as the A-inverse without genetic groups
    ggAinv[1:8, 1:8]
    noggAinv
   stopifnot(identical(ggAinv[1:8, 1:8], noggAinv))
   
}

