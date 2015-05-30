\name{makeA}
\alias{makeA}
\title{Creates the additive genetic relationship matrix
}
\description{This returns the additive relationship matrix in sparse matrix format.
}
\usage{
makeA(pedigree)
}
\arguments{
  \item{pedigree }{
A pedigree where the columns are ordered ID, Dam, Sire}
}

\details{
Missing parents (e.g., base population) should be denoted by either 'NA', '0', or '*'.

Used as a support function to \code{\link{makeD}}.

See function \code{\link{makeAinv}} for directly obtaining the inverse of the additive genetic relationship matrix.
}

\value{Returns A, or the numerator relationship matrix, in sparse matrix form.}

\author{\email{matthewwolak@gmail.com}}

\seealso{\code{\link{makeD}}, \code{\link{makeS}}}

\examples{
 makeA(Mrode2)
}

