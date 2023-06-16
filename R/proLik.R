#' Profile Likelihoods
#' 
#' Estimation, checking, and plotting of profile likelihoods and objects of
#' class \code{proLik} from a mixed model in ASReml-R.
#' 
#' For the \code{negative} argument, this should be used if the profile
#' likelihood of a covariance component is to be constructed.
#' 
#' \code{parallel} = TRUE should only be used on Linux or Mac OSes (i.e., not
#' Windows).
#' 
#' The function uses the \code{optimize} function to obtain the approximate
#' confidence limits. Therefore, \code{nse} should be carefully thought about
#' beforehand when running the function. Increasing this value will ensure the
#' confidence limits are contained by the search space, but at a cost to time.
#' 
#' \code{tolerance} is expressed in terms of the \code{alpha} value desired.
#' Convergence to a confidence limit will only be achieved when the constrained
#' value proposed for the confidence limit yields a likelihood ratio test
#' statistic that is no worse than \code{alpha + tolerance}.
#' 
#' If \code{negative} is FALSE, and the lower bound of the sampling interval
#' extends beyond zero, this will instead be set to effectively zero.
#' 
#' Obtaining the profile likelihood for residual variances may necessitate
#' explicitly specifying the R structure of the ASReml model. See example
#' below.
#' 
#' @aliases proLik proLik4 is.proLik plot.proLik
#' @param full.model An \code{asreml} model object
#' @param component A character (alternatively for \code{proLik4} this could also
#'   be a \code{formula}) indicating for which variance component the
#'   profile likelihood will be constructed. For \code{proLik}, must be an object
#'   in \code{full.model$gammas}. For \code{proLik4}, must be an object in 
#'   \code{full.model$vparameters}. To specify as a \code{formula}, components are
#'   written as \code{Vx}, where \dQuote{x} is a number between 1 and
#'   \code{length(full.model$vparameters)} (e.g., \code{component = ~ V1}).
#' @param G Logical indicating whether component is part of the G structure. If
#'   the component is part of the R structure, G = FALSE.
#' @param negative Logical indicating whether or not the \code{component} can
#'   take on a negative value (i.e., a covariance)
#' @param nsample.units Number of sample units to be used in constructing the
#'   area between the confidence limits for the profile likelihood
#' @param nse Number of standard errors on either side of the estimate, over
#'   which the confidence limits should be evaluated.
#' @param alpha The critical value for determining the Confidence Interval
#' @param tolerance Acceptable distance between actual and desired minimization
#'   criteria for determining the upper and lower limits of the confidence
#'   interval. See Details for more
#' @param parallel A logical indicating whether or not parallel processing will
#'   be used. Note, may only be available for Mac and Linux operating systems.
#' @param ncores Argument indicating number of cpu units to use. Default is all
#'   available.
#' @param x Object of class \code{proLik}.
#' @param CL A logical indicating whether a lines at the confidence limits are
#'   to be drawn when plotting.
#' @param type,main,xlab,ylab See arguments to \code{\link[graphics]{plot}}.
#' @param \dots other arguments to \code{\link[graphics]{plot}}.
#'
#' @return An S3 object of class \dQuote{proLik} containing:
#'   \describe{
#'     \item{lambdas }{negative log-Likelihood ratio test statistic. Estimated 
#'       from the log-Likelihood of the \code{full.model} and the 
#'       log-Likelihood of the model with the \code{component} constrained to 
#'       a value in the sampling interval}
#'     \item{var.estimates }{value along the sampling interval for which the 
#'       \code{component} was constrained}
#'     \item{UCL }{approximate Upper Confidence Limit}
#'     \item{LCL }{approximate Lower Confidence Limit}
#'     \item{component }{the component for which the profile likelihood surface 
#'       has been constructed}
#'     \item{alpha }{the alpha level at which the confidence limits have been 
#'       calculated}
#'   }
#' @section Warning : May be unfeasible to estimate profile likelihoods for
#'   complex models with many variance components
#' @author \email{matthewwolak@@gmail.com}
#' @seealso \code{\link{aiFun}}
#' @examples
#' 
#'   \dontrun{
#'     library(asreml)
#'     ginvA <- asreml.Ainverse(warcolak[, c(1,3,2)])$ginv
#'     ginvD <- makeD(warcolak[,1:3])$listDinv
#'     warcolak$IDD <- warcolak$ID
#'     warcolak.mod <- asreml(trait1 ~ sex, random = ~ ped(ID) + giv(IDD), 
#' 	ginverse = list(ID = ginvA, IDD = ginvD), rcov = ~ idv(units), data = warcolak) 
#'     summary(warcolak.mod)$varcomp
#'     profileA <- proLik(full.model = warcolak.mod, component = "ped(ID)!ped", 
#'         G = TRUE, negative = FALSE, nsample.units = 3, nse = 3)
#'     profileA
#'     profileD <- proLik(warcolak.mod, component = "giv(IDD).giv", 
#' 	G = TRUE, negative = FALSE, nsample.units = 3, nse = 3)
#'     profileE <- proLik(warcolak.mod, component = "R!units.var", G = FALSE, negative = FALSE)
#' 
#'     x11(w = 6, h = 8)
#'     par(mfrow = c(3,1))
#'       plot(profileA) 
#'       plot(profileD)
#'       plot(profileE)
#'    }
#' 
#' 
#' @export
proLik <- function(full.model, component,
	G = TRUE, negative = FALSE,
	nsample.units = 3, nse = 3,
	alpha = 0.05, tolerance = 0.001,
	parallel = FALSE, ncores = getOption("mc.cores", 2L)){
  s2 <- full.model$sigma2
  #TODO check if `component` is in  `paste0("V", seq(length(full.model$gammas)))`
  ## then use `V1` to `Vn` notation `else` code below
  gamma.ind <- which(names(full.model$gammas) == component)
  bound.check <- names(full.model$gammas.con)[gamma.ind]
    warned <- FALSE
    if(bound.check == "Boundary"){
      warned <- TRUE
      warning("Boundary parameter: confidence interval estimation may produce strange behavior - proceed with caution)")
      }
  gamma.est <- full.model$gammas[gamma.ind][[1]]
  if(full.model$gammas[gamma.ind] == summary(full.model)$varcomp[gamma.ind, 2]) s2 <- 1
  std.err <- sqrt(diag(aiFun(full.model))[gamma.ind])
  if(is.na(std.err) | std.err == 0) std.err <- 0.1 
  full.mod2 <- asreml::update.asreml(object = full.model, start.values = TRUE)$gammas.table
  chi.val <- 0.5 * qchisq(alpha, df = 1, lower.tail = FALSE)
  chi.tol <- chi.val - 0.5 * qchisq(alpha + tolerance, df = 1, lower.tail = FALSE)

  proLik_keep_uniQUe_UCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTU <- function(st, chi){
      proLik_keep_uniQUe_UCL[[1]] <<- c(proLik_keep_uniQUe_UCL[[1]], st)
      lrt <- constrainFun(st, full.model, full.mod2, component, G)
      proLik_keep_uniQUe_UCL[[2]] <<- c(proLik_keep_uniQUe_UCL[[2]], lrt)
      abs(chi - lrt)
  }
  proLik_keep_uniQUe_LCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTL <- function(st, chi){
      proLik_keep_uniQUe_LCL[[1]] <<- c(proLik_keep_uniQUe_LCL[[1]], st)
      lrt <- constrainFun(st, full.model, full.mod2, component, G)
      proLik_keep_uniQUe_LCL[[2]] <<- c(proLik_keep_uniQUe_LCL[[2]], lrt)
      abs(chi - lrt)
  }
  UCL <- LCL <- list(minimum = NA, objective = chi.val)
  ltol <- utol <- .Machine$double.eps^0.25
  cnt <- 0
  # Are 3 iterations of `while()` sufficient?
  ## with `optim(tol=...)` changed by 2 decimal places each successive try?
  while(((UCL[[2L]] > chi.tol) + (LCL[[2L]] > chi.tol)) > 0 & cnt < 3){
    if(cnt == 0){
      Uint <- c(gamma.est, gamma.est + (nse * std.err))
      Lint <- c(gamma.est - (nse * std.err), gamma.est)
    } else{
        if(UCL[[2L]] > chi.tol){
	  if(with(proLik_keep_uniQUe_UCL, any(lambdas > chi.val) && any(lambdas < chi.val))){
	    Uint <- with(proLik_keep_uniQUe_UCL, c(gam[which.min(ifelse((chi.val - lambdas) > 0, (chi.val - lambdas), Inf))], gam[which.min(ifelse((lambdas - chi.val) > 0, (lambdas - chi.val), Inf))]))
	  } else{
              Uext <- chi.val / max(proLik_keep_uniQUe_UCL$lambdas)
              Uint <- c(gamma.est + (nse * std.err), gamma.est + (Uext*(nse*std.err) + std.err))
            }
          # Adjust `optimize()` tolerance: see details in `?optimize`
#TODO could make `optimize(tol)` argument correspond more directly with `proLik(tolerance)`
## Use default `tol` in first pass
## then see, given approx. quadratic relationship between `gam` and `lambdas`
## what `optimize(tol)` needed to find `gam` value giving `lambda` within `tolerance`/`chi.tol` of `chi.val`
          x_0 <- with(proLik_keep_uniQUe_UCL, gam[which.min(abs(chi.val - lambdas))])
          d <- sqrt(.Machine$double.eps) * abs(x_0) + ((utol)/3)
          if(d > 1e-6) ltol <- 3*(d*0.01 - sqrt(.Machine$double.eps * x_0^2))
        }
        if(LCL[[2L]] > chi.tol){
	  if(with(proLik_keep_uniQUe_LCL, any(lambdas > chi.val) && any(lambdas < chi.val))){
	    Lint <- with(proLik_keep_uniQUe_LCL, c(gam[which.min(ifelse((lambdas - chi.val) > 0, (lambdas - chi.val), Inf))], gam[which.min(ifelse((chi.val - lambdas) > 0, (chi.val - lambdas), Inf))]))
	  } else{
              Lext <- chi.val / max(proLik_keep_uniQUe_LCL$lambdas)
              Lint <- c(gamma.est - (Lext*(nse*std.err) + std.err), gamma.est - (nse*std.err))
            }
          # Adjust `optimize()` tolerance: see details in `?optimize`
          x_0 <- with(proLik_keep_uniQUe_LCL, gam[which.min(abs(chi.val - lambdas))])
          d <- sqrt(.Machine$double.eps) * abs(x_0) + ((ltol)/3)
          if(d > 1e-6) ltol <- 3*(d*0.01 - sqrt(.Machine$double.eps * x_0^2))
        }          
      }
    if(!negative & Lint[1L] < 0) Lint[1L] <- 1e-8
      if(!negative & Lint[2L] < 0) Lint[2L] <- 1e-7
    #FIXME next two lines assume correlations and won't work for covariance
    if(negative == TRUE & Uint[2L] > 1) Uint[2L] <- 1.0 - 1e-8
    if(negative == TRUE & Lint[1L] < -1) Lint[1L] <- -1.0 + 1e-8

    if(UCL[[2L]] > chi.tol){
      if(parallel){
        tmpUCL <- parallel::mcparallel(expr = expression(c(optimize(f=tmpLRTU, interval = Uint, chi = chi.val, tol = utol), proLik_keep_uniQUe_UCL)))
      } else{
          UCL <- optimize(f = tmpLRTU, interval = Uint, chi = chi.val, tol = utol)
        }
    }
    if(LCL[[2L]] > chi.tol){
      LCL <- optimize(f = tmpLRTL, interval = Lint, chi = chi.val, tol = ltol)
    }
    if(parallel & UCL[[2L]] > chi.tol){
      tmpUCL.out <- parallel::mccollect(tmpUCL, wait = TRUE)[[1]]
      UCL <- list(minimum = tmpUCL.out$minimum[[1]], objective = tmpUCL.out$objective[[1]])
      proLik_keep_uniQUe_UCL <- list(gam = tmpUCL.out$gam, lambdas = tmpUCL.out$lambdas)
    }
    cnt <- cnt + 1
  }  #<-- end `while()`

  gamma.vec <- c(gamma.est, seq(gamma.est - (gamma.est - LCL$minimum)/2, gamma.est + (UCL$minimum - gamma.est)/2, length.out = nsample.units))

  if(parallel){
    if(length(gamma.vec) < ncores) ncores <- length(gamma.vec)
    prof <- list(lambdas = parallel::pvec(v = seq(1,length(gamma.vec),1), FUN = parConstrainFun, parameters = gamma.vec, full = full.model, fm2 = full.mod2, comp = component, G = G, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE), var.estimates = gamma.vec)
    } else{
    prof <- list(lambdas = vapply(gamma.vec, FUN = constrainFun, FUN.VALUE = vector("numeric", length = 1), full = full.model, fm2 = full.mod2, comp = component, G = G), var.estimates = gamma.vec)
      }

  prof$var.estimates <- c(prof$var.estimates, proLik_keep_uniQUe_LCL$gam, proLik_keep_uniQUe_UCL$gam)
  prof$lambdas <- c(prof$lambdas, proLik_keep_uniQUe_LCL$lambdas, proLik_keep_uniQUe_UCL$lambdas)
  ord.index <- order(prof$var.estimates)

    if(!negative & LCL$minimum < 0.01 & !warned){
      warning("Boundary parameter: confidence interval estimation may produce strange behavior - proceed with caution)")
      }

  # Remove CI limit estimates if they haven't been reached
  if(UCL[[2L]] > chi.tol) UCL[[1L]] <- NA
  if(LCL[[2L]] > chi.tol) LCL[[1L]] <- NA


 return(structure(list(lambdas = prof$lambdas[ord.index], 
	var.estimates = prof$var.estimates[ord.index] * s2, 
	UCL = UCL$minimum * s2, 
	LCL = LCL$minimum * s2, 
	component = component,
	alpha = alpha),
  class = c("proLik", class(prof))))
}










################################################################################
#' @rdname proLik
#' @export
proLik4 <- function(full.model, component,
	G = TRUE, negative = FALSE,
	nsample.units = 3, nse = 3,
	alpha = 0.05, tolerance = 0.001,
	parallel = FALSE, ncores = getOption("mc.cores", 2L)){

  if(inherits(component, "formula")){
    vin <- as.character(component[[length(component)]])
    v <- substr(vin, start = 1, stop = 1)
    vInd <- as.integer(substr(vin, start = 2, stop = nchar(vin)))
    if(!v %in% c("V", "v")){
      warning(cat(" Assuming formula position", length(component),
          "indicates component to profile.\n",
        "Also assuming right hand side of formula, character string position(s)",
          seq(nchar(vin))[-1],
          "\n represent an integer denoting component position in 'full.model$vparameters'\n",
          "(interpreted as position", vInd, ")\n"))
    }

  } else{
      if(!inherits(component, "character")){
        stop("'component' must be either a formula (e.g., ~ V1) or a character")
      }
      vInd <- match(component, names(full.model$vparameters))
    }  
   bndChk <- asreml::vpc.char(full.model)[[vInd]]
    warned <- FALSE
    if(bndChk == "B"){
      warned <- TRUE
      warning(cat("Boundary parameter: CI estimation may produce strange behavior\n
        proceed with caution\n"))
    }

  s2 <- full.model$sigma2
  vEst <- full.model$vparameters[[vInd]]
  std.err <- sqrt(diag(full.model$ai)[vInd])
  if(is.na(std.err) | std.err == 0) std.err <- 0.1 
  full.mod2 <- asreml::update.asreml(object = full.model,
    start.values = TRUE)$vparameters.table
  # constrain component
  full.mod2[vInd, "Constraint"] <- "F"
  chi.val <- 0.5 * qchisq(1 - alpha, df = 1)
  chi.tol <- chi.val - 0.5 * qchisq(alpha + tolerance, df = 1, lower.tail = FALSE)

  # Set asreml option to continue even if singularity in AI
  oldAIsingFailOpt <- asreml::asreml.options("ai.sing", "fail")
  asreml::asreml.options(ai.sing = TRUE, fail = "soft")
  on.exit(asreml::asreml.options(oldAIsingFailOpt))
  
  ###########################################################
  # Internal constrain function for asreml version 4
  ## So not confused with version 3 `nadiv::constrainFun()`
  conFun <- function(parameterVal, maxit = 600){ 
    full.mod2[vInd, "Value"] <- parameterVal
    ## rerun model with constraint
    if(G){
      conMod <- asreml::update.asreml(full.model,
        random = ~ ., G.param = full.mod2,
        maxiter = maxit, trace = FALSE)
    } else{
        conMod <- asreml::update.asreml(full.model,
          residual = ~ ., R.param = full.mod2,
          maxiter = maxit, trace = FALSE)
      }  #<-- end if/else


    cnt <- 0
    while(!conMod$converge & cnt <= 5){
      conMod <- asreml::update.asreml(conMod)
      cnt <- cnt + 1
    }  #<-- end while


    cnt <- 0
    if(conMod$converge){
      pcntChng <- sum(!conMod$vparameters.pc < 0.015, na.rm = TRUE) == 0
      while(!pcntChng & cnt <= 5){
        conMod <- asreml::update.asreml(conMod, maxiter = maxit)
        if(conMod$converge){
          pcntChng <- sum(!conMod$vparameters.pc < 0.015, na.rm = TRUE) == 0
        } else pcntChng <- FALSE
        cnt <- cnt + 1
      }  #<-- end while
      conMod$converge <- pcntChng
    }  #<-- end if converge

    if(conMod$converge) return(asreml::lrt(conMod, full.model)$'LR-statistic')
      else return(NA)
  }  #<-- end conFun()
  
  # Internal parallel friendly wrapper of constrain function for asreml version 4
  ## So not confused with version 3 `nadiv::parConstrainFun()`
  parConFun <- function(x, parameterVals){
    vapply(parameterVals[min(x):max(x)],
      FUN = conFun, FUN.VALUE = vector("numeric", length = 1))
  }  #<-- end parConFun() 
  #####################

  
  proLik_keep_uniQUe_UCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTU <- function(st, chi){
      proLik_keep_uniQUe_UCL[[1]] <<- c(proLik_keep_uniQUe_UCL[[1]], st)
      lrt <- conFun(st)
      proLik_keep_uniQUe_UCL[[2]] <<- c(proLik_keep_uniQUe_UCL[[2]], lrt)
      abs(chi - lrt)
  }
  proLik_keep_uniQUe_LCL <- list(gam = NULL, lambdas = NULL)
  tmpLRTL <- function(st, chi){
      proLik_keep_uniQUe_LCL[[1]] <<- c(proLik_keep_uniQUe_LCL[[1]], st)
      lrt <- conFun(st)
      proLik_keep_uniQUe_LCL[[2]] <<- c(proLik_keep_uniQUe_LCL[[2]], lrt)
      abs(chi - lrt)
  }
  UCL <- LCL <- list(minimum = NA, objective = chi.val)
  ltol <- utol <- .Machine$double.eps^0.25


  cnt <- 0
  # Are 3 iterations of `while()` sufficient?
  ## with `optim(tol=...)` changed by 2 decimal places each successive try?
  while(((UCL[[2L]] > chi.tol) + (LCL[[2L]] > chi.tol)) > 0 & cnt < 3){
    if(cnt == 0){
      Uint <- c(vEst, vEst + (nse * std.err))
      Lint <- c(vEst - (nse * std.err), vEst)
    } else{
        if(UCL[[2L]] > chi.tol){
	  if(with(proLik_keep_uniQUe_UCL,
	    any(lambdas > chi.val) && any(lambdas < chi.val))){
	      Uint <- with(proLik_keep_uniQUe_UCL,
	        c(gam[which.min(ifelse((chi.val - lambdas) > 0,
	            (chi.val - lambdas), Inf))],
	          gam[which.min(ifelse((lambdas - chi.val) > 0,
	            (lambdas - chi.val), Inf))]))
	  } else{
              Uext <- chi.val / max(proLik_keep_uniQUe_UCL$lambdas)
              Uint <- c(vEst + (nse * std.err),
                vEst + (Uext*(nse*std.err) + std.err))
            }
          # Adjust `optimize()` tolerance: see details in `?optimize`
#TODO could make `optimize(tol)` argument correspond more directly with `proLik(tolerance)`
## Use default `tol` in first pass
## then see, given approx. quadratic relationship between `gam` and `lambdas`
## what `optimize(tol)` needed to find `gam` value giving `lambda` within `tolerance`/`chi.tol` of `chi.val`
          x_0 <- with(proLik_keep_uniQUe_UCL,
            gam[which.min(abs(chi.val - lambdas))])
          d <- sqrt(.Machine$double.eps) * abs(x_0) + ((utol)/3)
          if(d > 1e-6) ltol <- 3*(d*0.01 - sqrt(.Machine$double.eps * x_0^2))
        }
        if(LCL[[2L]] > chi.tol){
	  if(with(proLik_keep_uniQUe_LCL,
	    any(lambdas > chi.val) && any(lambdas < chi.val))){
	      Lint <- with(proLik_keep_uniQUe_LCL,
	        c(gam[which.min(ifelse((lambdas - chi.val) > 0,
	            (lambdas - chi.val), Inf))],
	          gam[which.min(ifelse((chi.val - lambdas) > 0,
	            (chi.val - lambdas), Inf))]))
	  } else{
              Lext <- chi.val / max(proLik_keep_uniQUe_LCL$lambdas)
              Lint <- c(vEst - (Lext*(nse*std.err) + std.err),
                vEst - (nse*std.err))
            }
          # Adjust `optimize()` tolerance: see details in `?optimize`
          x_0 <- with(proLik_keep_uniQUe_LCL,
            gam[which.min(abs(chi.val - lambdas))])
          d <- sqrt(.Machine$double.eps) * abs(x_0) + ((ltol)/3)
          if(d > 1e-6) ltol <- 3*(d*0.01 - sqrt(.Machine$double.eps * x_0^2))
        }          
      }
    if(!negative & Lint[1L] < 0) Lint[1L] <- 1e-8
      if(!negative & Lint[2L] < 0) Lint[2L] <- 1e-7
    #FIXME next two lines assume correlations and won't work for covariance
    if(negative == TRUE & Uint[2L] > 1) Uint[2L] <- 1.0 - 1e-8
    if(negative == TRUE & Lint[1L] < -1) Lint[1L] <- -1.0 + 1e-8


    if(UCL[[2L]] > chi.tol){
      if(parallel){
        tmpUCL <- parallel::mcparallel(expr = expression(c(optimize(f=tmpLRTU,
            interval = Uint,
            chi = chi.val,
            tol = utol),
          proLik_keep_uniQUe_UCL)))
      } else{
          UCL <- optimize(f = tmpLRTU,
            interval = Uint, chi = chi.val, tol = utol)
        }
    }


    if(LCL[[2L]] > chi.tol){
      LCL <- optimize(f = tmpLRTL, interval = Lint, chi = chi.val, tol = ltol)
    }
    if(parallel & UCL[[2L]] > chi.tol){
      tmpUCL.out <- parallel::mccollect(tmpUCL, wait = TRUE)[[1]]
      UCL <- list(minimum = tmpUCL.out$minimum[[1]],
        objective = tmpUCL.out$objective[[1]])
      proLik_keep_uniQUe_UCL <- list(gam = tmpUCL.out$gam,
        lambdas = tmpUCL.out$lambdas)
    }
    cnt <- cnt + 1
  }  #<-- end `while()`

  vEstVec <- c(vEst,
      seq(vEst - (vEst - LCL$minimum)/2, vEst + (UCL$minimum - vEst)/2,
        length.out = nsample.units))


  if(parallel){
    if(length(vEstVec) < ncores) ncores <- length(vEstVec)
    prof <- list(lambdas = parallel::pvec(v = seq(1,length(vEstVec),1),
          FUN = parConFun, parameterVals = vEstVec,
          mc.set.seed = FALSE, mc.silent = FALSE,
          mc.cores = ncores, mc.cleanup = TRUE),
      var.estimates = vEstVec)
    } else{
    prof <- list(lambdas = vapply(vEstVec,
        FUN = conFun, FUN.VALUE = vector("numeric", length = 1)),
      var.estimates = vEstVec)
      }

  prof$var.estimates <- c(prof$var.estimates, proLik_keep_uniQUe_LCL$gam,
    proLik_keep_uniQUe_UCL$gam)
  prof$lambdas <- c(prof$lambdas, proLik_keep_uniQUe_LCL$lambdas,
    proLik_keep_uniQUe_UCL$lambdas)
  ord.index <- order(prof$var.estimates)

    if(!negative & LCL$minimum < 0.01 & !warned){
      warning(cat("Boundary parameter: CI estimation may produce strange behavior\n
        proceed with caution\n"))
      }

  # Remove CI limit estimates if they haven't been reached
  if(UCL[[2L]] > chi.tol) UCL[[1L]] <- NA
  if(LCL[[2L]] > chi.tol) LCL[[1L]] <- NA


 return(structure(list(lambdas = prof$lambdas[ord.index], 
	var.estimates = prof$var.estimates[ord.index] * s2, 
	UCL = UCL$minimum * s2, 
	LCL = LCL$minimum * s2, 
	component = component,
	alpha = alpha),
  class = c("proLik", class(prof))))
}

################################################################################



#' @method is proLik
#' @rdname proLik
#' @export
is.proLik <- function(x) inherits(x, "proLik")





#' @method plot proLik
#' @rdname proLik
#' @export
plot.proLik <- function(x, CL = TRUE, alpha = NULL, type = "l",
	main = NULL, xlab = NULL, ylab = NULL, ...)
{
  if(is.null(alpha)) alpha <- x$alpha
  if(is.null(main)) main <- deparse(substitute(x))
  if(is.null(xlab)) xlab <- x$component
  if(is.null(ylab)) ylab <- "LRT statistic"

  plot(x$lambdas ~ x$var.estimates,
     main = main, 
     xlab = xlab, ylab = ylab, 
     type = type, ...)
     if(CL){  
        chi <- (0.5 * qchisq(alpha, df = 1, lower.tail = FALSE))
        abline(h = chi, lty = "dotted", col = "red", lwd = 2)
        abline(v = unlist(x[c("LCL", "UCL")]), lty = "dashed", col = "blue", lwd = 2)
    }  
}

