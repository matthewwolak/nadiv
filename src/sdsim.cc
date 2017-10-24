#include "nadivcc.h"

extern "C"{  

void sdsim(int *da,
	int *sa,
        int *eN,
	int *en,
	int *dam,
	int *sire,
	int *Di,
	int *Dp,
	int *sdij
){         

  int i, j, k, m, p, l;
  int mi, si, cdama, csirea, rdama, rsirea;
  int adij = 0;
  GetRNGstate();

  for(i = 0; i < en[0]; i++){
     mi = dam[i];
     si = sire[i];
     if(mi != -999){
     k = i*eN[0];
        for(j = 0; j < eN[0]; j++){
           if(runif(0.0, 2.0) > 1.0){
              da[k] = da[(mi*eN[0]) + j];
           } 
              else {
                da[k] = sa[(mi*eN[0]) + j];
	      }
           k++;
        }
     }

     if(si != -999){
     k = i*eN[0];
        for(j = 0; j < eN[0]; j++){
           if(runif(0.0, 2.0) > 1.0){
	      sa[k] = da[(si*eN[0]) + j];
	   } 
              else {
	        sa[k] = sa[(si*eN[0]) + j];
              }
           k++;
        }
     }
  }  

  PutRNGstate();


  for(m = 0; m < en[0]; m++){
     for(p = Dp[m]; p < Dp[m+1]; p++){
        sdij[p] = 0;
        adij = 0;
        for(l = 0; l < eN[0]; l++){
           cdama = da[m*eN[0]+l];
           csirea = sa[m*eN[0]+l];
           rdama = da[Di[p]*eN[0]+l];
	   rsirea = sa[Di[p]*eN[0]+l];

           if(cdama == rdama){
              if(csirea == rsirea){
                 adij += 1;
              }
           } 
           else{
              if(cdama == rsirea){
                 if(csirea == rdama){
                    adij += 1;
                 }
              }
           }


        }

        sdij[p] += adij;

     }
  }               

}
}
