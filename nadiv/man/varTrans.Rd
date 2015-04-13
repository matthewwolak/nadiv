\name{varTrans}
\alias{varTrans}
\title{
Transforms ASReml-R gamma variances to component scale
}
\description{
The inverse of the Average Information matrix in an ASReml-R object produces the sampling variances of the random effects on the gamma scale.  This function scales these variances to the original component scale.  This allows for Confidence Intervals to be constructed about the variance component estimates.}

\usage{
varTrans(asr.object)
}
\arguments{
  \item{asr.object }{
Object from a call to \code{asreml}}
}

\value{
Returns a numeric vector of variances for each variance component in an ASReml-R model.}

\author{\email{matthewwolak@gmail.com}}



\examples{
  \dontrun{
    library(asreml)
    ginvA <- asreml.Ainverse(warcolak)$ginv
    ginvD <- makeD(warcolak[,1:3])$listDinv
    warcolak$IDD <- warcolak$ID
    warcolak.mod <- asreml(trait1 ~ 1, random = ~ped(ID) + giv(IDD), 
	ginverse = list(ID = ginvA, IDD = ginvD), data = warcolak) 
    summary(warcolak.mod)$varcomp
    varTrans(warcolak.mod)    
   }
}

