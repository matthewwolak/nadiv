\name{ggcontrib}
\alias{ggcontrib}
\title{
Genetic group contribution}
\description{
Calculates the genomic contribution each genetic group makes to every individual in a pedigree}
\usage{
ggcontrib(pedigree, ggroups = NULL, fuzz = NULL, output = "matrix")
}

\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire
}
  \item{ggroups }{
An optional vector of either: genetic group assignment for every individual or just the unique genetic groups. \code{fuzz} must be \code{NULL} if an object is supplied to the \code{ggroups} argument.
}
  \item{fuzz }{
A matrix containing the fuzzy classification of phantom parents into genetic groups. \code{ggroups} must be \code{NULL} if an object is supplied to the \code{fuzz} argument.
}
  \item{output }{
Format for the output
}
}

\details{
The specification of genetic groups is done in one of two approaches, either using fuzzy classification or not.

Fuzzy classification enables phantom parents to be assigned to (potentially) more than one genetic group (Fikse 2009). This method requires unique phantom parent identities to be included in the pedigree for all observed individuals with unknown parents. For 'p' phantom parents, 'p' identities should be listed as individuals in the first 'p' rows of the pedigree and these should be the only individuals in the pedigree with missing values in their Dam and Sire columns (denoted by either 'NA', '0', or '*'). The matrix supplied to the \code{fuzz} argument should have 'p' rows (one for each phantom parent) and 'r' columns, where 'r' is the number of genetic groups.

Non-fuzzy classification can handle the specification of genetic groups in three formats:

(1) similar to ASReml's format for specifying genetic groups, the first 'r' rows of the pedigree contain the label for each genetic group in the ID column and indicate missing values for the Dam and Sire columns (denoted by either 'NA', '0', or '*').  All individuals in the pedigree must then have one of the 'r' genetic groups as parent(s) for each unknown parent.  Note, a warning message indicating \code{In numPed(pedigree) : Dams appearing as Sires} is expected, since the dam and sire can be the same for all individuals in the pedigree composing the base population of a genetic group.

(2) similar to Jarrod Hadfield's \code{rbv} function arguments in the \code{MCMCglmm} package, for a pedigree of dimension i x 3 (specified by the \code{pedigree} argument), where 'i' is the total number of individuals in the pedigree, a similar vector of length 'i' can be specified by the \code{ggroups} argument.  This vector lists either the genetic group to which each individual belongs or NA if the individual is not to be considered part of one of the base populations (genetic groups).

(3) similar to DMU's format for specifying genetic groups. For a pedigree of dimension i x 3, where 'i' is the total number of individuals in the pedigree, instead of missing values for the parents the genetic groups are specified.  Note, that all individuals with a missing parent should have a genetic group substituted instead of the missing value symbol (i.e., either 'NA', '0', or '*').
}

\value{
Returns i x r genetic group contributions to all 'i' individuals from each of the 'r' genetic groups. Default output is a matrix (dense), but this format can be changed (e.g., "dgCMatrix" for a sparse matrix).
}
\references{
Fikse, F. 2009. Fuzzy classification of phantom parent groups in an animal model. Genetics, Selection, Evolution. 41:42.

Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding Values, 2nd ed.  Cambridge, MA: CABI Publishing.

Quaas, R.L. 1988. Additive genetic model with groups and relationships. Journal of Dairy Science. 71:1338-1345.
}

\author{\email{matthewwolak@gmail.com}}

\examples{
##########################
# Fuzzy classification
  ## Uses the pedigree from Quaas 1988 (See `data(Q1988)`)
  ## Fuzzy classification with complete assignment to one group
    Q1988fuzz <- Q1988[-c(1:2), c("id", "phantomDam", "phantomSire")]
    Qfnull <- matrix(c(1,0,0,1,0, 0,1,1,0,1), nrow = 5, ncol = 2,
	dimnames = list(letters[1:5], c("g1", "g2")))
    Qfuzznull <- ggcontrib(Q1988fuzz, fuzz = Qfnull)

    ## Should be identical to the non-fuzzy classification output
    # format (1) from above
      (Q <- ggcontrib(Q1988[-c(3:7), c(1,4,5)]))
    stopifnot(Qfuzznull == Q)

  ## Fuzzy classification with arbitrary assignments
    Qf <- matrix(c(1,0,0.5,0.5,0.5, 0,1,0.5,0.5,0.5), nrow = 5, ncol = 2,
	dimnames = list(letters[1:5], c("g1", "g2")))
    (Qfuzz <- ggcontrib(Q1988fuzz, fuzz = Qf))  

  ## Using the pedigree and fuzzy classification in Fikse (2009)
    F2009fuzz <- data.frame(id = c(letters[1:7], LETTERS[1:6]),
	dam = c(rep(NA, 7), "a", "c", "e", "A", "C", "D"),
	sire = c(rep(NA, 7), "b", "d", "f", "B", "g", "E"))
    Ff <- matrix(c(1,0,1,0,0,0,0.2,
		0,1,0,0.6,0,0.3,0.4,
		0,0,0,0.4,1,0.7,0.4),
		nrow = 7, ncol = 3,
		dimnames = list(letters[1:7], paste0("g", 1:3)))
    # Actual Q matrix printed in Fikse (2009)
      Fikse2009Q <- matrix(c(0.5,0.5,0,0.5,0.1,0.3,
			0.5,0.3,0.15,0.4,0.275,0.3375, 
			0,0.2,0.85,0.1,0.625,0.3625),
		nrow = 6, ncol = 3,
		dimnames = list(LETTERS[1:6], paste0("g", seq(3))))

    Ffuzz <- ggcontrib(F2009fuzz, fuzz = Ff)
    # Encountering some rounding error
      stopifnot(length((drop0(Ffuzz - Fikse2009Q, tol = 1e-12))@x) == 0)


##########################
# Non-fuzzy classification
## Below arbitrarily assigns genetic groups to the pedigree
## from Mrode (2005), chapter 2 (See `data(Mrode2)`).

  # format (1) from above
    ExPed1.tmp <- Mrode2
    ExPed1.tmp[c(1,2,4), 2] <- c("a", "b", "b")
    ExPed1.tmp[1:2, 3] <- c("a", "b")
    ExPed1 <- data.frame(id = c(letters[1:2], as.character(ExPed1.tmp$id)),
	dam = c(NA, NA, as.character(ExPed1.tmp$dam)),
	sire = c(NA, NA, as.character(ExPed1.tmp$sire)))
    gg1 <- ggcontrib(ExPed1) # note the warning message which is typical

  # format (2) from above
    ExPed2 <- Mrode2
    ggvec.in <- c("a", "b", NA, "b", NA, NA)
    (gg2 <- ggcontrib(ExPed2, ggroups = ggvec.in))
    str(as(gg2, "matrix"))

  # format (3) from above
    ExPed3 <- ExPed1.tmp
    gg3 <- ggcontrib(ExPed3, ggroups = c("a", "b"))

}


