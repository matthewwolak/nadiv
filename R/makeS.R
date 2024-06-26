#' Creates the additive genetic relationship matrix for the shared sex
#' chromosomes
#' 
#' The function returns the inverse of the additive relationship matrix in
#' sparse matrix format for the sex chromosomes (e.g., either X or Z).
#' 
#' Missing parents (e.g., base population) should be denoted by either 'NA',
#' '0', or '*'.
#' 
#' The inverse of the sex-chromosome additive genetic relationship matrix
#' (S-matrix) is constructed implementing the Meuwissen and Luo (1992)
#' algorithm to directly construct inverse additive relationship matrices
#' (borrowing code from Jarrod Hadfield's MCMCglmm function, \code{inverseA})
#' and using equations presented in Fernando & Grossman (1990; see Wolak et al.
#' 2013).  Additionally, the S-matrix itself can be constructed (although this
#' takes much longer than computing S-inverse directly).
#' 
#' The choices of dosage compensation models are: no global dosage compensation
#' ("ngdc"), random inactivation in the homogametic sex ("hori"), doubling of
#' the single shared sex chromosome in the heterogametic sex ("hedo"), halving
#' expression of both sex chromosomes in the homogametic sex ("hoha"), or
#' inactivation of the paternal sex chromosome in the homogametic sex ("hopi").
#' 
#' @param pedigree A pedigree where the columns are ordered ID, Dam, Sire, Sex
#' @param heterogametic Character indicating the label corresponding to the
#'   heterogametic sex used in the \dQuote{Sex} column of the pedigree
#' @param DosageComp A character indicating which model of dosage compensation.
#'   If \code{NULL} then the \dQuote{ngdc} model is assumed.
#' @param returnS Logical statement, indicating if the relationship matrix
#'   should be constructed in addition to the inverse
#'
#' @return a \code{list}:
#'   \describe{
#'     \item{model }{the model of sex-chromosome dosage compensation assumed.}
#'     \item{S }{the sex-chromosome relationship matrix in sparse matrix
#'       form or NULL if \code{returnS} = FALSE}
#'     \item{logDet }{the log determinant of the S matrix}
#'     \item{Sinv }{the inverse of the S matrix in sparse matrix form}
#'     \item{listSinv }{the three column form of the non-zero elements for the 
#'       inverse of the S matrix}
#'     \item{inbreeding }{the sex-linked inbreeding coefficients for all 
#'       individuals in the pedigree}
#'     \item{vii }{a vector of the (non-zero) elements of the diagonal V matrix
#'       of the S=TVT' decomposition. Contains the variance of Mendelian
#'       sampling for a sex-linked locus}
#'   }
#' @author \email{matthewwolak@@gmail.com}
#' @references Wolak, M.E., D.A. Roff, and D.J. Fairbairn. in prep. The
#' contribution of sex chromosomal additive genetic (co)variation to the
#' phenotypic resemblance between relatives under alternative models of dosage
#' compensation.
#' 
#' Fernando, R.L. & Grossman, M. 1990. Genetic evaluation with autosomal and
#' X-chromosomal inheritance. Theoretical and Applied Genetics, 80:75-80.
#' 
#' Meuwissen, T.H.E. and Z. Luo. 1992. Computing inbreeding coefficients in
#' large populations. Genetics, Selection, Evolution, 24:305-313.
#' @examples
#' 
#'  makeS(FG90, heterogametic = "0", returnS = TRUE)
#' 
#' @export
makeS <- function(pedigree, heterogametic,
	DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"),
	returnS = FALSE){

    if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")
    if(!heterogametic %in% pedigree[, 4]){
      stop("Error: value given to 'heterogametic' argument (", heterogametic,
        ") not a value in the 'Sex' column of the pedigree\n")
    }
    nPed <- numPed(pedigree[, 1:3])

    damsex <- pedigree[unique(nPed[, 2])[-1], 4]
    if(any(damsex == heterogametic)){
       pedname <- names(pedigree)
       pedigree <- pedigree[, c(1,3,2,4)]
       names(pedigree) <- pedname
       nPed <- numPed(pedigree[, 1:3])
      warning("Assuming female heterogametic (e.g., ZZ/ZW) sex chromosome system\n")
    } else warning("Assuming male heterogametic (e.g., XX/XY) sex chromosome system\n")

    sex <- rep(-998, dim(pedigree)[1])
    sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
    sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
    N <- dim(nPed)[1]
    N2 <- N + 1

    dc.model <- match.arg(DosageComp)
    if(is.null(dc.model)){
      warning("Assuming 'ngdc' dosage compensation model")
      dc.model <- "ngdc"
    }


    if(dc.model == "ngdc"){
          dnmiss <- which(nPed[,2] != -998)
          fsnmiss <- which(nPed[,3] != -998 & sex == 1)
          bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
          nA <- N + 2 * length(dnmiss) + 2 * length(fsnmiss)
          nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
          Q.col <- c(nPed[,1][dnmiss], nPed[,1][fsnmiss], 1:N) 
          Q.row <- c(nPed[,2][dnmiss], nPed[,3][fsnmiss], 1:N)
          Q.x <- c(rep(-0.5, length(dnmiss)), rep(-1, length(fsnmiss)), rep(1, N))
          ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
          Q <- sparseMatrix(i = as.integer(Q.row[ord] - 1),
	    p = as.integer(c(match(1:N, Q.col[ord]), length(ord) + 1) - 1),
            x = as.double(Q.x[ord]),
            index1 = FALSE, dims = c(N, N), symmetric = FALSE,
            dimnames = list(NULL, NULL))

          nPed[nPed == -998] <- N2
          Vii <- (sex + 1)/2
          f <- c(rep(0, N), -1)

          Cout <- .C("sinv", PACKAGE = "nadiv",
	    as.integer(nPed[, 2] - 1), #dam
	    as.integer(nPed[, 3] - 1),  #sire
	    as.double(f),  #f
	    as.double(Vii),  #vii
            as.integer(Q@i),  #iQP
	    as.integer(c(Q@p, length(Q@x))),  #pQP
            as.double(Q@x),  #xQP
            as.integer(N),   #nQP
	    as.integer(length(Q@x)), #nzmaxQP	
	    as.integer(rep(0, nA)), #iSP
	    as.integer(rep(0, N + 1)), #pSP
	    as.double(rep(0, nA)), #xSP
	    as.integer(nA), #nzmaxSP
            as.double(0.25), #DC
            as.double(Vii))  #sex

          Sinv <- sparseMatrix(i = Cout[[10]][1:Cout[[13]]],
	    p = Cout[[11]],
            x = Cout[[12]][1:Cout[[13]]],
            index1 = FALSE, dims = c(N, N),
            dimnames = list(as.character(pedigree[, 1]),
                           as.character(pedigree[, 1])))
          Vii <- Cout[[4]]
          f <- Cout[[3]][1:N]



    } else{
         if(dc.model != "hopi"){
             fdnmiss <- which(nPed[,2] != -998 & sex == 1)
             mdnmiss <- which(nPed[,2] != -998 & sex == 0)
             fsnmiss <- which(nPed[,3] != -998 & sex == 1)
             bnmiss <- which(nPed[, 2] != -998 & nPed[, 3] != -998)
             nA <- N + 2 * (length(fdnmiss) + length(mdnmiss)) + 2 * length(fsnmiss)
             nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[bnmiss]) == FALSE)
             Q.col <- c(nPed[,1][fdnmiss], nPed[,1][mdnmiss], nPed[,1][fsnmiss], 1:N) 
             Q.row <- c(nPed[,2][fdnmiss], nPed[,2][mdnmiss], nPed[,3][fsnmiss], 1:N)
             Q.x <- c(rep(-0.5, length(fdnmiss)), rep(-1, length(mdnmiss)), rep(-0.5, length(fsnmiss)), rep(1, N))
             ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
             Q <- sparseMatrix(i = as.integer(Q.row[ord] - 1),
	       p = as.integer(c(match(1:N, Q.col[ord]), length(ord)+1) - 1),
               x = as.double(Q.x[ord]),
               index1 = FALSE, dims = c(N, N), symmetric = FALSE,
               dimnames = list(NULL, NULL))

             nPed[nPed == -998] <- N2
             Vii <- (2 - sex)
             f <- c(rep(0, N), -1)

             Cout <- .C("sinv", PACKAGE = "nadiv",
	       as.integer(nPed[, 2] - 1), #dam
	       as.integer(nPed[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(1), #DC
               as.double(Vii))  #sex

             Sinv <- sparseMatrix(i = Cout[[10]][1:Cout[[13]]],
	       p = Cout[[11]],
               x = Cout[[12]][1:Cout[[13]]],
               index1 = FALSE, dims = c(N, N),
               dimnames = list(as.character(pedigree[, 1]),
                           as.character(pedigree[, 1])))
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]


         } else{
                dnmiss <- which(nPed[,2] != -998)
                nA <- N + 2 * length(dnmiss)
                nA <- nA + 2 * sum(duplicated(paste(nPed[, 2], nPed[, 3])[dnmiss]) == FALSE)
                Q.col <- c(nPed[,1][dnmiss], 1:N) 
                Q.row <- c(nPed[,2][dnmiss], 1:N)
                Q.x <- c(rep(-0.5, length(dnmiss)), rep(1, N))
                ord <- order(Q.col + Q.row/(N+1), decreasing = FALSE)
                Q <- sparseMatrix(i = as.integer(Q.row[ord] - 1),
	          p = as.integer(c(match(1:N, Q.col[ord]), length(ord)+1) - 1),
                  x = as.double(Q.x[ord]),
                  index1 = FALSE, dims = c(N, N), symmetric = FALSE,
                  dimnames = list(NULL, NULL))

                nPed[nPed == -998] <- N2
                Vii <- rep(1, N) 
                f <- rep(0, N)

             Cout <- .C("sinv", PACKAGE = "nadiv",
	       as.integer(nPed[, 2] - 1), #dam
	       as.integer(nPed[, 3] - 1),  #sire
	       as.double(f),  #f
	       as.double(Vii),  #vii
               as.integer(Q@i),  #iQP
	       as.integer(c(Q@p, length(Q@x))),  #pQP
               as.double(Q@x),  #xQP
               as.integer(N),   #nQP
	       as.integer(length(Q@x)), #nzmaxQP	
	       as.integer(rep(0, nA)), #iSP
	       as.integer(rep(0, N + 1)), #pSP
	       as.double(rep(0, nA)), #xSP
	       as.integer(nA), #nzmaxSP
               as.double(0), #DC
               as.double(Vii))  #sex

             Sinv <- sparseMatrix(i = Cout[[10]][1:Cout[[13]]],
	       p = Cout[[11]],
               x = Cout[[12]][1:Cout[[13]]],
               index1 = FALSE, dims = c(N, N),
               dimnames = list(as.character(pedigree[, 1]),
                           as.character(pedigree[, 1])))
             Vii <- Cout[[4]]
             f <- Cout[[3]][1:N]

           }
      }              

    logDet <- -1*determinant(Sinv, logarithm = TRUE)$modulus[1]
    listSinv <- sm2list(Sinv, rownames = as.character(pedigree[, 1]),
                              colnames = c("Row", "Column", "Sinverse"))
    if(returnS){
          T <- as(solve(Q), "CsparseMatrix")
          S <- as(crossprod(T, Diagonal(N, Vii)) %*% T, "dgCMatrix")
            S@Dimnames <- list(as.character(pedigree[, 1]), NULL)
    } else{
         S <- NULL
      }

return(list(model = dc.model, S = S, logDet = logDet,
	Sinv = Sinv, listSinv = listSinv,
	inbreeding = f, vii = Vii))
}

