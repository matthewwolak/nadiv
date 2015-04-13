varTrans <- function(asr.object){
    if(all(summary(asr.object)$varcomp$gamma == summary(asr.object)$varcomp$component)){
       vars <- asr.object$ai
    } else{
       Rcomp <- which(asr.object$gammas == 1.00)
       AI <- aiFun(asr.object)
       comps <- summary(asr.object)$varcomp$component

       vars <- c(((asr.object$gammas[-Rcomp]^2) * diag(AI)[Rcomp] + comps[Rcomp]^2 * diag(AI)[-Rcomp] + 2*asr.object$gammas[-Rcomp]*comps[Rcomp]*AI[Rcomp, -Rcomp]), diag(AI)[Rcomp])
      }
vars
}

