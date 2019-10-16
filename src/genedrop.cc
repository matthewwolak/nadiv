#include "nadivcc.h"

extern "C"{  

void genedrop(int *da,     // Ndalleles
	int *sa,           // Nsalleles
        int *eN,           // N
	int *en,	   // n; pedigree size
	int *dam,      
	int *sire
){         

  int i, j, k, l, mi, si;

  GetRNGstate();
  for(i = 0; i < en[0]; i++){
    mi = dam[i];
    si = sire[i];
    if(mi != -999){
      k = i*eN[0];
      l = mi*eN[0];
      for(j = 0; j < eN[0]; j++){
        if(runif(0.0, 2.0) > 1.0){
          da[k] += da[l];
        } 
          else {
            da[k] += sa[l];
        }
        k++;
        l++;
      }
    }

    if(si != -999){
      k = i*eN[0];
      l = si*eN[0];
      for(j = 0; j < eN[0]; j++){
        if(runif(0.0, 2.0) > 1.0){
	  sa[k] += da[l];
        } 
          else {
            sa[k] += sa[l];
          }
        k++;
        l++;
      }
    }
  } 
  PutRNGstate();

}
}

