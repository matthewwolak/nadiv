\name{ggcontrib}
\alias{ggcontrib}
\title{
Genetic group contribution}
\description{
Calculates the genomic contribution each genetic group makes to every individual in a pedigree}
\usage{
ggcontrib(pedigree, ggroups = NULL, fuzz = NULL)
}

\arguments{
  \item{pedigree}{
A pedigree where the columns are ordered ID, Dam, Sire
}
  \item{ggroups}{
An optional vector of either: genetic group assignment for every individual or just the unique genetic groups}
  \item{fuzz}{
A matrix containing the fuzzy classification of individuals into genetic groups.  Not yet implemented}
}

\details{
This function can handle the specification of genetic groups in three formats:

(1) similar to ASReml's format for specifying genetic groups, the first 'n' rows of the pedigree contain the label for each genetic group in the ID column and indicate missing values for the Dam and Sire columns (denoted by either 'NA', '0', or '*').  All individuals in the pedigree must then have one of the 'n' genetic groups as parent(s) for each unknown parent.  Note, a warning message indicating \code{In numPed(pedigree) : Dams appearing as Sires} is expected, since the dam and sire should be the same for all individuals in the pedigree composing the base population of a genetic group.

(2) similar to Jarrod Hadfield's \code{rbv} function arguments in the \code{MCMCglmm} package, for a pedigree of dimension i x 3 (specified by the \code{pedigree} argument), where 'i' is the total number of individuals in the pedigree, a similar vector of length 'i' can be specified by the \code{ggroups} argument.  This vector lists either the genetic group to which each individual belongs or NA if the individual is not to be considered part of one of the base populations (genetic groups).

(3) similar to DMU's format for specifying genetic groups, for a pedigree of dimension i x 3 (specified by the \code{pedigree} argument), where 'i' is the total number of individuals in the pedigree, instead of missing values for the parents, the genetic groups are specified.  Note, that all individuals with a missing parent should have a genetic group substituted instead of the missing value symbol (i.e., either 'NA', '0', or '*').
}

\value{
Returns i x n sparse matrix ("dgCMatrix") of genetic contributions to all 'i' individuals from each of the 'n' genetic groups.
}
\references{
Fikse, F. 2009. Fuzzy classification of phantom parent groups in an animal model. Genetics, Selection, Evolution. 41:42.
}

\author{\email{matthewwolak@gmail.com}}

\examples{
#### Below arbitrarily assigns genetic groups
### to the pedigree from Mrode (2005), chapter 2
### (See data(Mrode2).

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
    gg2 <- ggcontrib(ExPed2, ggroups = ggvec.in)
    as(gg2, "matrix")  # changes object to a dense matrix 
       #for easier incorporation into dataframes etc.
    str(as(gg2, "matrix"))

  # format (3) from above
    ExPed3 <- ExPed1.tmp
    gg3 <- ggcontrib(ExPed3, ggroups = c("a", "b"))

}


