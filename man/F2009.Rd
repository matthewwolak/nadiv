\name{F2009}
\alias{F2009}
\docType{data}
\title{
Pedigree, adapted from Fikse (2009)
}
\description{
An example pedigree with genetic groups and fuzzy classification of genetic groups.
}
\usage{data("F2009")}
\format{
  A data frame with 16 observations on the following 11 variables.
  \describe{
    \item{\code{id} }{a factor with levels indicating the unique individuals (including phantom parents) and genetic groups}
    \item{\code{dam} }{a factor of observed maternal identities}
    \item{\code{sire} }{a factor vector of observed paternal identities}
    \item{\code{damGG} }{a factor of maternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{sireGG} }{a factor of paternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{phantomDam} }{a factor of maternal identities with phantom parents inserted instead of \code{NA}}
    \item{\code{phantomSire} }{a factor of paternal identities with phantom parents inserted instead of \code{NA}}
    \item{\code{group} }{a factor of genetic groups to which each phantom parent belongs}
    \item{\code{g1} }{a numeric vector with probabilities of group \code{g1} membership for each phantom parent}
    \item{\code{g2} }{a numeric vector with probabilities of group \code{g2} membership for each phantom parent}
    \item{\code{g3} }{a numeric vector with probabilities of group \code{g3} membership for each phantom parent}
  }
}


\source{
Fikse, F. 2009. Fuzzy classification of phantom parent groups in an animal model. Genetics Selection Evolution 41:42.
}

\keyword{datasets}
