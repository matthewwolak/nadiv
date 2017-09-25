#include "nadivcc.h"
//////////////////////////////////////
// since nadiv >v2.15.0 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 

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

  int     lb, step, it, k, j, m, kDam, kSire, jDam, jSire;
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

             m = pAP[max(kDam, jSire)];
             lb = pAP[max(kDam, jSire)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kDam, jSire)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


             m = pAP[max(kDam, jDam)];
             lb = pAP[max(kDam, jDam)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kDam, jDam)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


             m = pAP[max(kSire, jDam)];
             lb = pAP[max(kSire, jDam)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kSire, jDam)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


             m = pAP[max(kSire, jSire)];
             lb = pAP[max(kSire, jSire)+1] - 1 - m;
             while(lb > 0){
               step = lb/2;
               it = m + step;
               if(iAP[it] < min(kSire, jSire)){
                 m=++it;
                 lb-=step+1;
               }
               else lb = step;
             }
             if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


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
// since nadiv >v2.14.3 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 

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

  int     lb, step, it, k, m, kDam, kSire, jDam, jSire;
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


           m = pAP[max(kDam, jSire)];
           lb = pAP[max(kDam, jSire)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kDam, jSire)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


           m = pAP[max(kDam, jDam)];
           lb = pAP[max(kDam, jDam)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kDam, jDam)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


           m = pAP[max(kSire, jDam)];
           lb = pAP[max(kSire, jDam)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kSire, jDam)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


           m = pAP[max(kSire, jSire)];
           lb = pAP[max(kSire, jSire)+1] - 1 - m;
           while(lb > 0){
             step = lb/2;
             it = m + step;
             if(iAP[it] < min(kSire, jSire)){
               m=++it;
               lb-=step+1;
             }
             else lb = step;
           }
           if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


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















//////////////////////////////////////
//    SEX-LINKED DOMINANCE
// since nadiv >v2.15.0 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 

extern "C"{  

void sdij(
	int *dam,
	int *sire,
        int *iAP,     
	int *pAP,	         
	double *xAP,
	int *nAP,
	double *dij,
	int *Di,
        int *Dp,
	int *cnt,
	int *sex
){         

  int     lb, step, it, k, j, m, kDam, kSire, jDam, jSire;
  double rmmp, rffp, rmfp, rfmp, dij_tmp;

  for(k = 0; k < nAP[0]; k++){  // iterate through each column of "S" called "A"
    if(sex[k] == 1){
      Dp[k] = cnt[0];
      kDam = dam[k];
      kSire = sire[k];
      if((kDam != -999) && (kSire != -999)){
        for(j = pAP[k]; j < pAP[k+1]; j++){ //iterate through all rows of k column
           if((k != iAP[j]) && (sex[iAP[j]] == 1)){ 
             jDam = dam[iAP[j]];
             jSire = sire[iAP[j]];
             if((jDam != -999) && (jSire != -999)){
               rmfp = 0.0;
               rmmp = 0.0;
               rfmp = 0.0;
               rffp = 0.0;

               m = pAP[max(kDam, jSire)];
               lb = pAP[max(kDam, jSire)+1] - 1 - m;
               while(lb > 0){
                 step = lb/2;
                 it = m + step;
                 if(iAP[it] < min(kDam, jSire)){
                   m=++it;
                   lb-=step+1;
                 }
                 else lb = step;
               }
               if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


               m = pAP[max(kDam, jDam)];
               lb = pAP[max(kDam, jDam)+1] - 1 - m;
               while(lb > 0){
                 step = lb/2;
                 it = m + step;
                 if(iAP[it] < min(kDam, jDam)){
                   m=++it;
                   lb-=step+1;
                 }
                 else lb = step;
               }
               if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


               m = pAP[max(kSire, jDam)];
               lb = pAP[max(kSire, jDam)+1] - 1 - m;
               while(lb > 0){
                 step = lb/2;
                 it = m + step;
                 if(iAP[it] < min(kSire, jDam)){
                   m=++it;
                   lb-=step+1;
                 }
                 else lb = step;
               }
               if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


               m = pAP[max(kSire, jSire)];
               lb = pAP[max(kSire, jSire)+1] - 1 - m;
               while(lb > 0){
                 step = lb/2;
                 it = m + step;
                 if(iAP[it] < min(kSire, jSire)){
                   m=++it;
                   lb-=step+1;
                 }
                 else lb = step;
               }
               if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


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
    }  // if k has homogametic sex chromosomes
  }  // end for k loop   
         

            
}
}















//////////////////////////////////////////
// since nadiv >v2.14.3 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 

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

  int     lb, step, it, k, j, m, kDam, kSire, jDam, jSire, jmoDam, jmoSire;
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

                  m = pAP[max(kDam, jSire)];
                  lb = pAP[max(kDam, jSire)+1] - 1 - m;
                  while(lb > 0){
                    step = lb/2;
                    it = m + step;
                    if(iAP[it] < min(kDam, jSire)){
                      m=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iAP[m] == min(kDam, jSire)) rmfp += xAP[m];


                  m = pAP[max(kDam, jDam)];
                  lb = pAP[max(kDam, jDam)+1] - 1 - m;
                  while(lb > 0){
                    step = lb/2;
                    it = m + step;
                    if(iAP[it] < min(kDam, jDam)){
                      m=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iAP[m] == min(kDam, jDam)) rmmp += xAP[m];


                  m = pAP[max(kSire, jDam)];
                  lb = pAP[max(kSire, jDam)+1] - 1 - m;
                  while(lb > 0){
                    step = lb/2;
                    it = m + step;
                    if(iAP[it] < min(kSire, jDam)){
                      m=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iAP[m] == min(kSire, jDam)) rfmp += xAP[m];


                  m = pAP[max(kSire, jSire)];
                  lb = pAP[max(kSire, jSire)+1] - 1 - m;
                  while(lb > 0){
                    step = lb/2;
                    it = m + step;
                    if(iAP[it] < min(kSire, jSire)){
                      m=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iAP[m] == min(kSire, jSire)) rffp += xAP[m];


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

