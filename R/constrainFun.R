#' Function used in conjunction with others to produce a profile likelihood for
#' a variance component
#' 
#' Given a model object from \code{asreml} and a range of estimates of the
#' parameter, the function will supply the likelihood ratio test statistic for
#' the comparison of the full model to one where the parameter of interest is
#' constrained.
#' 
#' Used internally in the \code{proLik} function
#' 
#' @param parameter.val a value for which the log-Likelihood of a model is to
#'   be calculated
#' @param full the full model \code{asreml} object
#' @param fm2 starting values for the full model
#' @param comp which variance component to constrain
#' @param G logical, indicating if the component is part of the G structure
#' @param mit numeric, indicating maximum number of iterations for the
#' constrained asreml model
#'
#' @author \email{matthewwolak@@gmail.com}
#' @seealso See also \code{\link{proLik}}
#' @export
constrainFun <- function(parameter.val, full, fm2, comp, G, mit = 600){
  row <- which(fm2$Gamma == comp)
  fm2[row, 2:3] <- c(parameter.val, "F")
  if(G) full$G.param <- fm2 else full$R.param <- fm2
  con.mod <- asreml::update.asreml(object = full, maxiter = mit, trace = FALSE)
  cnt <- 0
  while(!con.mod$converge & cnt <= 5){
     con.mod <- asreml::update.asreml(con.mod)
     cnt <- cnt + 1
  }
  cnt <- 0
  if(con.mod$converge){
     pcc.out <- pcc(con.mod, silent = TRUE)
     while(!pcc.out & cnt <= 5){
        con.mod <- asreml::update.asreml(con.mod, maxiter = mit)
        if(con.mod$converge) pcc.out <- pcc(con.mod, silent = TRUE)
        cnt <- cnt + 1
     }
     con.mod$converge <- pcc.out
  }
  if(con.mod$converge) return(LRTest(full$loglik, con.mod$loglik)$lambda) else return(NA)
}

