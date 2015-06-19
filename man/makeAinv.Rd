\name{makeAinv}
\alias{makeAinv}
\title{Creates the Inverse additive genetic relationship matrix
}
\description{This returns the inverse of the additive relationship matrix. It can also be used to obtain coefficients of inbreeding for the pedigreed population.
}
\usage{
makeAinv(pedigree, ggroups = NULL, fuzz = NULL, keepPhantoms = TRUE,
      gOnTop = FALSE, det = FALSE)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
 \item{ggroups }{
Either a vector with the names of the genetic groups, or a numeric indicating the number of genetic groups. See Details for different ways to specify. Note, if NULL then the regular A-inverse will be constructed.
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

The function implements an adaptation of the Meuwissen and Luo (1992) algorithm with some code borrowed from the \code{inverseA} function by Jarrod Hadfield in the \code{MCMCglmm} package. Further, providing a non-NULL argument to \code{ggroups} implements the Quaas (1988) algorithm for directly obtaining the augmented A-inverse matrix for genetic groups. Combination of Meuwissen and Luo's (1992) and Quaas's (1988) algorithm incorporates inbreeding into the construction of the A-inverse.

#TODO: fuzz, gOnTop, keepPhantoms
}

\value{
  \item{Ainv }{the inverse of the additive genetic relationship matrix in sparse matrix form}
  \item{listAinv }{the three column form of the non-zero elements for the inverse of the additive genetic relationship matrix}
  \item{f }{the inbreeding coefficients for each individual in the pedigree}
  \item{logDet }{the log determinant of the A matrix}
}

\references{
Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.

Quaas, R.L. 1988. Additive genetic model with groups and relationships. Journal of Dairy Science. 71:1338-1345.
}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{makeA}}}

\examples{
 makeAinv(Mrode2)
}

