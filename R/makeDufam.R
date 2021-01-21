makeDufam <- function(pedigree, parallel = FALSE,
	ncores = getOption("mc.cores", 2L), invertD = TRUE,
	returnA = FALSE, det = TRUE){

  N <- nrow(pedigree)
  pedigree <- cbind(pedigree, gen = genAssign(pedigree), oseq = seq.int(N))
  pedigree <- pedigree[order(pedigree$gen, pedigree[, 2], pedigree[, 3]), ]
  numeric.pedigree <- numPed(pedigree[, 1:3]) 
  A <- makeA(pedigree[, 1:3])
  dA <- diag(A)

  if(parallel){
     if(length(A@x)/ncores < 10){
        warning("pedigree too small - 'parallel' set to FALSE instead")
        parallel <- FALSE
     }
  }

  if(!parallel){
     cat(paste("starting to make D..."))
     Cout <- .C("dijjskip", PACKAGE = "nadiv",
                as.integer(numeric.pedigree[, 2] - 1), 
		as.integer(numeric.pedigree[, 3] - 1), 
		as.integer(A@i), 			
		as.integer(A@p),                        
		as.double(A@x/2),                       
		as.integer(N),                           
		as.double(rep(0, length(A@x))),         
		as.integer(rep(0, length(A@i))),        
                as.integer(rep(0, N)),                  
		as.integer(0))	                        

     D <- sparseMatrix(i = Cout[[8]][1:Cout[[10]]],
     		p = c(Cout[[9]], Cout[[10]]),
     		x = Cout[[7]][1:Cout[[10]]],
     		dims = c(N, N),	dimnames = list(as.character(pedigree[, 1]), NULL),
     		symmetric = TRUE, index1 = FALSE)
    diag(D) <- 2 - dA

     if(!returnA) A <- NULL
     rm("Cout")

   } else{
        listA <- data.frame(Row = as.integer(rep(1:length(A@p[-1]), diff(A@p))), Column = as.integer(A@i + 1))
        wrap_dij <- function(x){
           sub_lA <- listA[min(x):max(x), 1:2]
           lA_r <- dim(sub_lA)[1]
           Cout <- .C("dijp", PACKAGE = "nadiv",
		as.integer(numeric.pedigree[, 2] - 1),
		as.integer(numeric.pedigree[, 3] - 1), 
                as.integer(lA_r),  
		as.integer(sub_lA[, 1] - 1),  
		as.integer(sub_lA[, 2] - 1),  
		as.integer(A@i),  
		as.integer(A@p), 
		as.double(A@x/2),  
		as.double(rep(0, lA_r)))
         Cout[[9]]
        }

        cat(paste("starting to make D..."))
        Dijs <- parallel::pvec(seq(1, dim(listA)[1], 1), FUN = wrap_dij, mc.set.seed = FALSE, mc.silent = FALSE, mc.cores = ncores, mc.cleanup = TRUE)
  
        D <- sparseMatrix(i = A@i,
     		p = A@p,
     		x = Dijs,
     		dims = c(N, N),	dimnames = list(as.character(pedigree[, 1]), NULL),
     		symmetric = TRUE, index1 = FALSE)
        if(!returnA) A <- NULL
        D <- drop0(D)
        diag(D) <- 2 - dA

     }

  cat(paste(".done", "\n"))
  D <- D[pedigree$oseq, pedigree$oseq]
  if(returnA) A <- A[pedigree$oseq, pedigree$oseq]
  if(det) logDet <- determinant(D, logarithm = TRUE)$modulus[1] else logDet <- NULL
  if(invertD){
    Dinv <- as(solve(D), "dgCMatrix")
    Dinv@Dimnames <- list(as.character(pedigree[pedigree$oseq, 1]), NULL)
    listDinv <- sm2list(Dinv, rownames = pedigree[pedigree$oseq, 1], colnames=c("row", "column", "Dinverse"))
 return(list(A = A, D = D, logDet = logDet, Dinv=Dinv, listDinv=listDinv))
  } else{
    return(list(A = A, D = D, logDet = logDet))
    } 
}

