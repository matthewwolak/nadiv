\name{Q1988}
\alias{Q1988}
\docType{data}
\title{
Pedigree, adapted from Quaas (1988) equation [5]}
\description{
An example pedigree with genetic groups}
\usage{data("Q1988")}
\format{
  A data frame with 11 observations on the following 8 variables.
  \describe{
    \item{\code{id} }{a factor with levels indicating the unique individuals (including phantom parents) and genetic groups}
    \item{\code{dam} }{a factor of observed maternal identities}
    \item{\code{sire} }{a factor vector of observed paternal identities}
    \item{\code{damGG} }{a factor of maternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{sireGG} }{a factor of paternal identities with genetic groups inserted instead of \code{NA}}
    \item{\code{phantomDam} }{a factor of maternal identities with phantom parents inserted instead of \code{NA}}
    \item{\code{phantomSire} }{a factor of paternal identities with phantom parents inserted instead of \code{NA}}
    \item{\code{group} }{a factor of genetic groups to which each phantom parent belongs}
  }
}
\source{
Quaas, R.L. 1988. Additive genetic model with groups and relationships. Journal of Dairy Science 71:1338-1345.}
\keyword{datasets}
