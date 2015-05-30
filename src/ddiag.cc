#include "ddiag.h"

extern "C"{  

void ddiag(
        int *dam,       
        int *sire,         
        double *f,     
        double *dii,    
	int *iTP,              
	int *pTP,	         
	double *xTP,
	int *nTP,
 	int *nzmaxTP
){         

  int     k, j, l;
  double  sumTcol;
  cs *T;

  T = cs_spalloc(nTP[0], nTP[0], nzmaxTP[0], true, false);  

         for (k = 0 ; k < nzmaxTP[0] ; k++){
           T->i[k] = iTP[k];
           T->x[k] = xTP[k];
         }
         for (k = 0 ; k <= nTP[0]; k++){
           T->p[k] = pTP[k];
         }

  for(k=0; k < nTP[0]; k++){  // iterate through each row of T 
      dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);
      sumTcol=0.0;      // set "sumTcol" (sum of a column of T) to zero

      for(j = T->p[k]; j < T->p[k+1]; j++){
          l = T->i[j];
          if(l <= k){
              sumTcol += T->x[j] * T->x[j] * dii[l];
          }
     }

     f[k] = sumTcol-1.0;
   }
            
  cs_spfree(T);
}
}
