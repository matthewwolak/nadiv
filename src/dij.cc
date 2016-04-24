#include "nadivcc.h"

extern "C"{  

void dij(
	int *dam,
	int *sire,
        int *iAP,     
	int *pAP,	         
	double *xAP,
	int *nAP,
	double *dij,
	int *Di,
        int *Dp,
	int *cnt
){         

  int     k, j, m, kDam, kSire, jDam, jSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;

  for(k = 0; k < nAP[0]; k++){  // iterate through each column of "A"
    Dp[k] = cnt[0];
    kDam = dam[k];
    kSire = sire[k];
    if((kDam != -999) && (kSire != -999)){
      for(j = pAP[k]; j < pAP[k+1]; j++){ //iterate through all rows of k column
         if(k != iAP[j]){ 
           jDam = dam[iAP[j]];
           jSire = sire[iAP[j]];
           if((jDam != -999) && (jSire != -999)){
             rmfp = 0.0;
             rmmp = 0.0;
             rfmp = 0.0;
             rffp = 0.0;

             for(m = pAP[max(kDam, jSire)]; m < pAP[max(kDam, jSire)+1]; m++){
                if(iAP[m] == min(kDam, jSire)){
                  rmfp += xAP[m];
               }
                if(rmfp != 0.0){
                  break;
                }
             }

             for(m = pAP[max(kDam, jDam)]; m < pAP[max(kDam, jDam)+1]; m++){
                if(iAP[m] == min(kDam, jDam)){
                   rmmp += xAP[m];
                }
                if(rmmp != 0.0){
                   break;
                }
             }

             for(m = pAP[max(kSire, jDam)]; m < pAP[max(kSire, jDam)+1]; m++){  
                if(iAP[m] == min(kSire, jDam)){
                   rfmp += xAP[m];
                }
                if(rfmp != 0.0){
	           break;
                }
             }

             for(m = pAP[max(kSire, jSire)]; m < pAP[max(kSire, jSire)+1]; m++){
                if(iAP[m] == min(kSire, jSire)){
                   rffp += xAP[m];
                }
                if(rffp != 0.0){
                   break;
                }
             }


             dij_tmp = (rmmp*rffp) + (rmfp*rfmp);
             if(dij_tmp != 0.0){
                dij[cnt[0]] = dij_tmp;   
                Di[cnt[0]] = iAP[j];
                cnt[0]++;
             }
           }   
         }

      }   
    }   
  }   
         

            
}
}




/////////////////////////////////////////////


extern "C"{  

void dijp(
	int *dam,
	int *sire,
	int *lAr,
	int *indk,
        int *indj,
	int *iAP,     
	int *pAP,	         
	double *xAP,
	double *dij
){         

  int     k, m, kDam, kSire, jDam, jSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;

  for(k = 0; k < lAr[0]; k++){ 
    kDam = dam[indk[k]];
    kSire = sire[indk[k]];
    if((kDam != -999) && (kSire != -999)){
       if(indk[k] != indj[k]){ 
         jDam = dam[indj[k]];
         jSire = sire[indj[k]];
         if((jDam != -999) && (jSire != -999)){
           rmfp = 0.0;
           rmmp = 0.0;
           rfmp = 0.0;
           rffp = 0.0;
           for(m = pAP[max(kDam, jSire)]; m < pAP[max(kDam, jSire)+1]; m++){
             if(iAP[m] == min(kDam, jSire)){
               rmfp += xAP[m];
             }
             if(rmfp != 0.0){
                break;
             }
           }

           for(m = pAP[max(kDam, jDam)]; m < pAP[max(kDam, jDam)+1]; m++){
             if(iAP[m] == min(kDam, jDam)){
               rmmp += xAP[m];
             }
             if(rmmp != 0.0){
                break;
             }
           }

           for(m = pAP[max(kSire, jDam)]; m < pAP[max(kSire, jDam)+1]; m++){  
             if(iAP[m] == min(kSire, jDam)){
               rfmp += xAP[m];
             }
             if(rfmp != 0.0){
	        break;
             }
           }

           for(m = pAP[max(kSire, jSire)]; m < pAP[max(kSire, jSire)+1]; m++){
             if(iAP[m] == min(kSire, jSire)){
               rffp += xAP[m];
             }
             if(rffp != 0.0){
                break;
             }
           }


           dij_tmp = (rmmp*rffp) + (rmfp*rfmp);
           if(dij_tmp != 0.0){
              dij[k] = dij_tmp;   
           }   
         }
       }   
    }   
  }   
                    
}
}








//////////////////////////////////////////

extern "C"{  

void dijjskip(
	int *dam,
	int *sire,
        int *iAP,     
	int *pAP,	         
	double *xAP,
	int *nAP,
	double *dij,
	int *Di,
        int *Dp,
	int *cnt
){         

  int     k, j, m, kDam, kSire, jDam, jSire, jmoDam, jmoSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;
 
  dij_tmp = 0.0;
  for(k = 0; k < nAP[0]; k++){  // iterate through each column of "A"
    Dp[k] = cnt[0];
    kDam = dam[k];
    kSire = sire[k];
    if((kDam != -999) && (kSire != -999)){
      jmoDam = -999;
      jmoSire = -999;
      for(j = pAP[k]; j < pAP[k+1]; j++){ //iterate through all rows of k column
         if(k != iAP[j]){ 
           jDam = dam[iAP[j]];
           jSire = sire[iAP[j]];
           if((jDam != -999) && (jSire != -999)){
             if((jDam == jmoDam) && (jSire == jmoSire)){
                if(dij_tmp != 0.0){
                   dij[cnt[0]] = dij_tmp;   
                   Di[cnt[0]] = iAP[j];
                   cnt[0]++;
                }
             }
                else{
                  rmfp = 0.0;
                  rmmp = 0.0;
                  rfmp = 0.0;
                  rffp = 0.0;

                  for(m = pAP[max(kDam, jSire)]; m < pAP[max(kDam, jSire)+1]; m++){
                     if(iAP[m] == min(kDam, jSire)){
                       rmfp += xAP[m];
                    }
                     if(rmfp != 0.0){
                       break;
                     }
                  }

                  for(m = pAP[max(kDam, jDam)]; m < pAP[max(kDam, jDam)+1]; m++){
                     if(iAP[m] == min(kDam, jDam)){
                        rmmp += xAP[m];
                     }
                     if(rmmp != 0.0){
                        break;
                     }
                  }

                  for(m = pAP[max(kSire, jDam)]; m < pAP[max(kSire, jDam)+1]; m++){  
                     if(iAP[m] == min(kSire, jDam)){
                        rfmp += xAP[m];
                     }
                     if(rfmp != 0.0){
	                break;
                     }
                  }

                  for(m = pAP[max(kSire, jSire)]; m < pAP[max(kSire, jSire)+1]; m++){
                     if(iAP[m] == min(kSire, jSire)){
                        rffp += xAP[m];
                     }
                     if(rffp != 0.0){
                        break;
                     }
                  }


                  jmoDam = jDam;
                  jmoSire = jSire;
                  dij_tmp = (rmmp*rffp) + (rmfp*rfmp);
                  if(dij_tmp != 0.0){
                     dij[cnt[0]] = dij_tmp;   
                     Di[cnt[0]] = iAP[j];
                     cnt[0]++;
                  }
                }
           }   
         }

      }   
    }   
  }   
         

            
}
}
