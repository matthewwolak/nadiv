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

