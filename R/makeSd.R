#' @aliases makeSd makeD
#' @rdname makeD
#' @export
makeSd <- function(pedigree, heterogametic,
	DosageComp = c(NULL, "ngdc", "hori", "hedo", "hoha", "hopi"),
	parallel = FALSE, ncores = getOption("mc.cores", 2L),
	invertSd = TRUE, returnS = FALSE, det = TRUE, verbose = TRUE){


  if(length(unique(pedigree[,4])) > 2) stop("Error: more than 2 sexes specified")

  dc.model <- match.arg(DosageComp)
  if(is.null(dc.model)){
    warning("Assuming 'ngdc' dosage compensation model")
    dc.model <- "ngdc"
  }
  if(dc.model == "hopi" | dc.model == "hori"){
    warning("Assume sex chromosomal dominance allelic interactions do not occur under 'hopi' or 'hori'\n")
    return(NULL)
  }


  Sout <- makeS(pedigree, heterogametic = heterogametic,
	DosageComp = dc.model, returnS = TRUE)
  # makeA() returns `dsCMatrix`, but S is `dgCMatrix` from above
  ## makeD()-like code below expects symmetric matrix ('dsCMatrix')
  S <- forceSymmetric(Sout$S)

  nPed <- numPed(pedigree[, 1:3])
  damsex <- pedigree[unique(nPed[, 2])[-1], 4]
  if(any(damsex == heterogametic)){
    pedname <- names(pedigree)
    pedigree <- pedigree[, c(1,3,2,4)]
    names(pedigree) <- pedname
    nPed <- numPed(pedigree[, 1:3])
  }

  sex <- rep(-998, dim(pedigree)[1L])
  sex[homs <- which(pedigree[,4] != heterogametic)] <- 1
  sex[hets <- which(pedigree[,4] == heterogametic)] <- 0
  nhom <- sum(sex)  # Number of individuals with homogametic sex chromosomes
  N <- dim(nPed)[1L]

#FIXME turned off next check so can test parallel=TRUE on small pedigrees
#  if(parallel){
#    if(length(S@x)/ncores < 10){
#      warning("pedigree too small - 'parallel' set to FALSE instead")
#      parallel <- FALSE
#    }
#  }

  if(!parallel){
     if(verbose) cat("starting to make Sd...")

     Cout <- .C("sdij", PACKAGE = "nadiv",
                as.integer(nPed[, 2] - 1), 		# [[1]] dam ID/No.
		as.integer(nPed[, 3] - 1), 		# [[2]] sire ID/No.
		as.integer(S@i), 			# [[3]] S@i
		as.integer(S@p),                        # [[4]] S@p
		as.double(S@x),                         # [[5]] S@x
		as.integer(N),                          # [[6]] No. in pedigree
		as.double(rep(0, length(S@x))),         # [[7]] Sd@x
		as.integer(rep(0, length(S@i))),        # [[8]] Sd@i
                as.integer(rep(0, N)),                  # [[9]] Sd@p
		as.integer(0),				# [[10]] cnt/count
		as.integer(sex))	      		# [[11]] sex                  

     Sd <- sparseMatrix(i = Cout[[8]][1:Cout[[10]]],
     		p = Cout[[9]][1:(nhom+1)],
     		x = Cout[[7]][1:Cout[[10]]],
     		dims = c(nhom, nhom),
     		dimnames = list(as.character(pedigree[homs, 1]), NULL),
     		symmetric = TRUE, index1 = FALSE)
     diag(Sd) <- 1 - Sout$inbreeding[homs]

     if(!returnS) S <- NULL
     rm("Cout")

   } else{
#TODO
stop("code not yet written to parallelize function") #FIXME
#        listA <- data.frame(Row = as.integer(rep(1:length(A@p[-1]), diff(A@p))), Column = as.integer(A@i + 1))
#        wrap_dij <- function(x){
#           sub_lA <- listA[min(x):max(x), 1:2]
#           lA_r <- dim(sub_lA)[1]
#           Cout <- .C("dijp", PACKAGE = "nadiv",
#		as.integer(numeric.pedigree[, 2] - 1),
#		as.integer(numeric.pedigree[, 3] - 1), 
#                as.integer(lA_r),  
#		as.integer(sub_lA[, 1] - 1),  
#		as.integer(sub_lA[, 2] - 1),  
#		as.integer(A@i),  
#		as.integer(A@p), 
#		as.double(A@x/2),  
#		as.double(rep(0, lA_r)))
#         Cout[[9]]
#        }

#        if(verbose) cat("starting to make D...")
#        Dijs <- parallel::pvec(seq(1, dim(listA)[1], 1), FUN = wrap_dij, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
  
#        D <- sparseMatrix(i = A@i, p = A@p, x = Dijs, ...)
#        if(!returnA) A <- NULL
#        D <- drop0(D)
#        diag(D) <- 2 - dA

     }

  if(verbose) cat(".done", "\n")
  
  if(det) logDet <- determinant(Sd, logarithm = TRUE)$modulus[1] else logDet <- NULL
  if(invertSd){
    if(verbose) cat("starting to invert Sd...")
    Sdinv <- as(solve(Sd), "dgCMatrix")
      Sdinv@Dimnames <- Sd@Dimnames
    if(verbose) cat(".done", "\n")
    listSdinv <- sm2list(Sdinv, rownames = Sd@Dimnames[[1L]],
	colnames = c("row", "column", "Sdinverse"))
 return(list(S = S, Sd = Sd, logDet = logDet,
		Sdinv = Sdinv, listSdinv = listSdinv))
  } else{
    return(list(S = S, Sd = Sd, logDet = logDet))
    } 
}


