\name{aiCI}
\alias{aiCI}

\title{
Confidence Intervals for Variance Components
}

\description{
Produces the 1-alpha Upper and Lower Confidence Limits for the variance components in an ASReml-R model.}

\usage{
aiCI(asr.model, Dimnames = NULL, alpha = 0.05)
}

\arguments{
  \item{asr.model }{
Object from a call to \code{asreml}
}
  \item{Dimnames }{
A vector of characters if names are desired for the output.  If not specified, the default labels from the \code{asreml} object will be used.
}
  \item{alpha }{
A numeric value indicating the level of Type I error for constructing the Confidence Intervals.
}}

\details{
Variances from the inverse of the Average Information matrix of an ASReml model are translated according to the \code{\link{varTrans}} function and used in constructing the 1-alpha Confidence Interval.}

\value{
A matrix is returned with a row for each variance component.  The three columns correspond to the Lower Confidence Limit, estimate from the \code{asreml} model, and Upper Confidence Limit for each variance component.}

\author{\email{matthewwolak@gmail.com}}

\note{
The vector of \code{Dimnames} should match the same order of variance components specified in the model.}


\seealso{
\code{\link{aiFun}}
\code{\link{proLik}}
}

\examples{
  \dontrun{
    library(asreml)
    ginvA <- asreml.Ainverse(warcolak)$ginv
    ginvD <- makeD(warcolak[,1:3])$listDinv
    warcolak$idd <- warcolak$id
    warcolak.mod <- asreml(trait1 ~ sex, random = ~ped(id) + giv(idd), 
	ginverse = list(id = ginvA, idd = ginvD), data = warcolak) 
    summary(warcolak.mod)$varcomp
    aiCI(warcolak.mod)    
   }
}

