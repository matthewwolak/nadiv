#include "slsims.h"

extern "C"{  

void slsims(
	int *dam,
	int *sire,
	int *sex,
	int *n,
	int *nloci,
	int *seg,
	double *hapa,
	double *hapb
){         

  int     k, l, m, kdam, ksire, cnt;

  cnt = 0;
  for(k = 0; k < n[0]; k++){
     kdam = dam[k];
     ksire = sire[k];
     if((kdam != -999)){
        for(l = 0; l < nloci[0]; l++){
           m = k*nloci[0] + l;
           if(seg[cnt] == 1){
              hapa[m] += hapa[kdam*nloci[0] + l];
           }
             else{
                hapa[m] += hapb[kdam*nloci[0] + l];
             }

           if(sex[k] == 2){
              hapb[m] += hapa[ksire*nloci[0] + l];
           }
        cnt++;
        }
     }
  }
        
}
}
