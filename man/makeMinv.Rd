% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeMinv.R
\name{makeMinv}
\alias{makeMinv}
\alias{makeMinvML}
\title{Create the inverse (additive) mutational effects relationship matrix}
\usage{
makeMinv(pedigree, ...)

makeMinvML(pedigree, ...)
}
\arguments{
\item{pedigree}{A pedigree where the columns are ordered ID, Dam, Sire}

\item{\dots}{Arguments to be passed to methods}
}
\value{
a \code{list}:
  \describe{
    \item{Minv }{the inverse of the (additive) mutational effects
      relationship matrix in sparse matrix form}
    \item{listMinv }{the three column list of the non-zero elements for the 
      inverse of the (additive) mutational effects relationship matrix.
      \code{attr(*, "rowNames")} links the integer for rows/columns to the ID
      column from the pedigree.}
    \item{h }{the amount by which segregation variance is reduced by
      inbreeding. Similar to the individual coefficients of inbreeding (f)
      derived during the construction of the inverse numerator relatedness matrix.
      in the pedigree (matches the order of the first/ID column of the
      pedigree).}
    \item{logDet }{the log determinant of the M matrix}
    \item{dii }{the (non-zero) elements of the diagonal D matrix of the M=TDT'
      decomposition. Contains the variance of Mendelian sampling. Matches
      the order of the first/ID column of the pedigree. Note Wray (1990) and
      Casellas and Medrano (2008) algorithms use \code{v=sqrt(dii)}.} 
  }
}
\description{
Returns the inverse of the (additive) mutational effects relationship matrix.
It can also be used to obtain components needed for the calculations in the
underlying algorithm.
}
\details{
Missing parents (e.g., base population) should be denoted by either 'NA',
'0', or '*'.

Note the assumption under the infinitesimal model, that mutation has essentially
zero probability of affecting an inbred locus (hence removing inbred
identity-by-descent), however, mutations may themselves be subject to
inbreeding (Wray 1990).

By default, the algorithm described in Casellas and Medrano (2008) is
implemented here, in which the inverse-M is separate from the typical inverse
relatedness matrix (inverse-A). Casellas and Medrano's algorithm allows
separate partitioning of additive genetic variance attributed to inheritance
of allelic variation present in the base population (inverse-A) from
additive genetic variance arising from mutation and subsequent sharing of
mutant alleles identical-by-descent. Alternatively, Wray (1990) formulates
an algorithm which combines both of these processes (i.e., the A-inverse with
the M-inverse matrices). If the Wray algorithm is desired, this can be
implemented by specifying a numeric value to an argument named \code{theta}.
The value used for \code{theta} should be as described in Wray (1990). See
examples below for use of this argument.
}
\examples{

 ##  Example pedigree from Wray 1990
 #### Implement Casellas & Medrano (2008) algorithm
   Mout <- makeMinv(Wray90[, 1:3])
 #### Wray (1990) algorithm with extra argument `theta`
   Mwray <- makeMinv(Wray90[, 1:3], theta = 10.0)$Minv # compare to Wray p.184
}
\references{
Casellas, J. and J.F. Medrano. 2008. Within-generation mutation
variance for litter size in inbred mice. Genetics. 179:2147-2155. 

Meuwissen, T.H.E & Luo, Z. 1992. Computing inbreeding
coefficients in large populations. Genetics, Selection, Evolution. 24:305-313.

Mrode, R.A. 2005. Linear Models for the Prediction of Animal Breeding
Values, 2nd ed.  Cambridge, MA: CABI Publishing.

Wray, N.A. 1990. Accounting for mutation effects in the additive genetic
variance-covariance matrix and its inverse. Biometrics. 46:177-186.
}
\author{
\email{matthewwolak@gmail.com}
}
