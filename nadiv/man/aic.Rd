\name{aic}
\alias{aic}
\title{
Akaike Information Criterion}
\description{
Calculates AIC/AICc values, AIC differences, Likelihood of models, and model probabilities.
}
\usage{
aic(logLik, fp, n = NULL)
}
\arguments{
  \item{logLik }{
A vector of model log-Likelihoods
}
  \item{fp }{
A vector containing the numbers of free parameters of each model included in the logLik vector
}
  \item{n }{
An optional vector of sample sizes for each model.  Used to calculate AICc (small sample un-biased AIC).
}
}

\details{
Calculations and notation follows chapter 2 of Burnham and Anderson (2002).
}

\value{
  \item{AIC }{vector containing AIC/AICc (depending on value of \code{n}}
  \item{delta_AIC }{vector containing AIC differences from the minimum AIC(c)}
  \item{AIClik }{vector containing likelihoods for each model, given the data.  Represents the relative strength of evidence for each model.}
  \item{w }{Akaike weights.}

}
\references{
Burnham, K.P. and D.R. Anderson. 2002. Model Selection and Multimodel Inference. A Practical Information-Theoretic Approach, 2nd edn. Springer, New York.
}

\author{\email{matthewwolak@gmail.com}
}

\examples{
   aic(c(-3139.076, -3136.784, -3140.879, -3152.432), c(8, 7, 8, 5)) 
}

