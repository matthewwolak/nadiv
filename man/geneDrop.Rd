\name{geneDrop}
\alias{geneDrop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to conduct gene dropping through a pedigree
}

\description{
Functions that perform and summarize gene dropping conducted on supplied pedigrees
}
\usage{
geneDrop(pedigree, N, parallel = FALSE, ncores = getOption("mc.cores", 2L), ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{pedigree }{
A pedigree with columns organized: ID, Dam, Sire}
  \item{N }{
The number of times to iteratively trace alleles through the pedigree}
  \item{parallel }{
A logical indicating whether or not to use parallel processing.  Note, this may only be available for Mac and Linux operating systems.}
  \item{ncores }{
The number of cpus to use when constructing the dominance relatedness matrix.  Default is all available.}
  \item{\dots }{
Other arguments that can be supplied to alter what summaries are reported.}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

\code{parallel} = TRUE should only be used on Linux or Mac operating systems (i.e., not Windows).

Founder allelic values (the alleles assigned to an individual's maternal, paternal, or both haplotypes when the maternal, paternal, or both parents are missing) are equivalent positive and negative integer values corresponding to the maternal and paternal haplotypes, respectively. For example, if the first individual in the pedigree has two unknown parents it will have the following two allelic values: 1=maternal haplotype and -1=paternal haplotype.  
}

\value{
  \item{IDs }{Original identities in the pedigree}
  \item{maternal }{Simulated maternal haplotypes}
  \item{paternal }{Simulated paternal haplotypes}
  \item{numericPedigree }{Pedigree in class \code{numPed} for convenient post-processing of haplotypes}
}

\author{\email{matthewwolak@gmail.com}
}

\seealso{
 \code{\link{makeDsim}} 
}

\examples{
  geneDrop(Mrode2, N = 10)
}
