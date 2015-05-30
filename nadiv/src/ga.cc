#include "ga.h"

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
     if((kdam != -999) & (ksire != -999)){
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

  int     k, dk, sk;

  for(k = 0; k < n[0]; k++){
     dk = k;
     sk = k;
     while(dam[dk] != -999){
        dk = dam[dk];
        dgen[k] += 1;
     }
     while(sire[sk] != -999){
        sk = sire[sk];
        sgen[k] += 1;
     }
  }
          
}
}
