#include "nadivcc.h"

extern "C"{  

void ga(
	int *dam,
	int *sire,
	int *generation,
	int *n
){         

  int     k, kdam, ksire;

  for(k = 0; k < n[0]; k++){
     kdam = dam[k];
     ksire = sire[k];
     if((kdam != -999) && (ksire != -999)){
        generation[k] = max(generation[kdam], generation[ksire]) + 1;
     }
     else{
        if((kdam != -999)){
           generation[k] = generation[kdam] + 1;
        }
        if((ksire != -999)){
           generation[k] = generation[ksire] + 1;
        }
     }
  }
        
}
}


//////////////////////////////////////

extern "C"{  

void gaUnsort(
	int *dam,
	int *sire,
	int *dgen,
	int *sgen,
	int *n
){         

  int     en, k, dk, sk;

  en = n[0];  // initialize flag value for infinite loops
  for(k = 0; k < n[0]; k++){
     dk = k;
     sk = k;
     while(dam[dk] != -999){
        dk = dam[dk];
        dgen[k] += 1;
        if(dgen[k] > n[0]){  // catch any infinite loops
          en = k;  // stick ID involved in ped loop into this spot to return to R
          dgen[0] = -999;  // flag that ped loop occurred through dam ancestors
          break;  // end while early
        }
     }
     if(en < n[0])  break;  // should end for loop early

     while(sire[sk] != -999){
        sk = sire[sk];
        sgen[k] += 1;
        if(sgen[k] > n[0]){  // catch any infinite loops
          en = k;  // stick ID involved in ped loop into this spot to return to R
          sgen[0] = -999;  // flag that ped loop occurred through sire ancestors
          break;  // end while early
        }
     }
     if(en < n[0]) break;  // should end for loop early
       
  }  // end for k

  n[0] = en;

}
}
