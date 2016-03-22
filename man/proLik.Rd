\name{proLik}
  \alias{proLik}
  \alias{is.proLik}
  \alias{plot.proLik}

\title{
Profile Likelihoods
}

\description{
Estimation, checking, and plotting of profile likelihoods and objects of class \code{proLik} from a mixed model in ASReml-R.  
}

\usage{
  proLik(full.model, component, G = TRUE, negative = FALSE, 
       nsample.units = 3, nse = 3, alpha = 0.05, tolerance = 0.001,
       parallel = FALSE, ncores = getOption("mc.cores", 2L))

  \method{is }{proLik}(x)

  \method{plot }{proLik}(x, CL = TRUE, alpha = 0.05, type = "l", \dots)
}

\arguments{
  \item{full.model }{An \code{asreml} model object}
  \item{component }{A character indicating for which variance component the profile likelihood will be constructed. Must be an object in \code{full.model$gammas}.}
  \item{G }{Logical indicating whether component is part of the G structure. If the component is part of the R structure, G = FALSE.}
  \item{negative }{Logical indicating whether or not the \code{component} can take on a negative value (i.e., a covariance)}
  \item{nsample.units }{Number of sample units to be used in constructing the area between the confidence limits for the profile likelihood}
  \item{nse }{Number of standard errors on either side of the estimate, over which the confidence limits should be evaluated.}
  \item{alpha }{The critical value for determining the Confidence Interval}
  \item{tolerance }{Acceptable distance, between actual sample values and interpolated values, for determining the upper and lower limits of the Confidence Interval. Actual sample points will be no more than this distance from the true value of the estimate.}
  \item{parallel }{A logical indicating whether or not parallel processing will be used. Note, may only be available for Mac and Linux operating systems.}
  \item{ncores }{Argument indicating number of cpu units to use. Default is all available.}
  \item{x }{Object of class \code{proLik}.}
  \item{CL }{A logical indicating whether a line representing the Confidence Limit is to be drawn when plotting.}
  \item{type }{The type of plot to be generated, see arguments to \code{\link[graphics]{plot}}.}
  \item{\dots }{other arguments to \code{\link[graphics]{plot}}.}

}


\details{
For the \code{negative} argument, this should be used if the profile likelihood of a covariance component is to be constructed.

\code{parallel} = TRUE should only be used on Linux or Mac OSes (i.e., not Windows).

The function uses the \code{optimize} function to obtain the approximate confidence limits. Therefore, \code{nse} should be carefully thought about beforehand when running the function. Increasing this value will ensure the confidence limits are contained by the search space, but at a cost to time. 

If \code{negative} is FALSE, and the lower bound of the sampling interval extends beyond zero, this will instead be set to effectively zero.

Obtaining the profile likelihood for residual variances may necessitate explicitly specifying the R structure of the ASReml model. See example below.    
}


\value{
An S3 object of class proLik
  \item{lambdas }{negative log-Likelihood ratio test statistic. Estimated from the log-Likelihood of the \code{full.model} and the log-Likelihood of the model with the \code{component} constrained to a value in the sampling interval}
  \item{var.estimates }{value along the sampling interval for which the \code{component} was constrained}
  \item{UCL }{approximate Upper Confidence Limit}
  \item{LCL }{approximate Lower Confidence Limit}
  \item{component }{the component for which the profile likelihood surface has been constructed}
}

\author{\email{matthewwolak@gmail.com}}

\section{Warning }{
May be unfeasible to estimate profile likelihoods for complex models with many variance components
}

\seealso{
\code{\link{aiFun}}
}

\examples{
  \dontrun{
    library(asreml)
    ginvA <- asreml.Ainverse(warcolak[, c(1,3,2)])$ginv
    ginvD <- makeD(warcolak[,1:3])$listDinv
    warcolak$IDD <- warcolak$ID
    warcolak.mod <- asreml(trait1 ~ sex, random = ~ ped(ID) + giv(IDD), 
	ginverse = list(ID = ginvA, IDD = ginvD), rcov = ~ idv(units), data = warcolak) 
    summary(warcolak.mod)$varcomp
    profileA <- proLik(full.model = warcolak.mod, component = "ped(ID)!ped", 
        G = TRUE, negative = FALSE, nsample.units = 3, nse = 3)
    profileA
    profileD <- proLik(warcolak.mod, component = "giv(IDD).giv", 
	G = TRUE, negative = FALSE, nsample.units = 3, nse = 3)
    profileE <- proLik(warcolak.mod, component = "R!units.var", G = FALSE, negative = FALSE)

    x11(w = 6, h = 8)
    par(mfrow = c(3,1))
      plot(profileA) 
      plot(profileD)
      plot(profileE)
   }

}


