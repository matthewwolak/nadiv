#include "ddiag.h"

extern "C"{  

void sinv(
        int *dam,       
        int *sire,         
        double *f,     
        double *vii,    
	int *iQP,              
	int *pQP,	         
	double *xQP,
        int *nQP,
	int *nzmaxQP,
	int *iSP,
	int *pSP,
	double *xSP,
	int *nzmaxSP,
	double *DC,
	double *sex
){         

  int     j, k, h, ntwo, cnt, sj, dj;
  double  si;
  double  *AN = new double[2*nQP[0]];
  double  *li = new double[nQP[0]];
  cs *Q, *V, *tQ, *QV, *tS, *S;

  for(k=0; k<nQP[0]; k++){
     li[k]=0.0;               
  }
  for(k=0; k<nQP[0]; k++){
     AN[k]=-1;               
  }

  ntwo = nQP[0];

  Q = cs_spalloc(nQP[0], nQP[0], nzmaxQP[0], true, false);  
         for (k = 0 ; k < nzmaxQP[0] ; k++){
           Q->i[k] = iQP[k];
           Q->x[k] = xQP[k];
         }
         for (k = 0 ; k <= nQP[0]; k++){
           Q->p[k] = pQP[k];
         }
   
  tQ = cs_transpose(Q, true);

  V = cs_spalloc(nQP[0], nQP[0], nzmaxQP[0], true, false);
	 for (k = 0; k < nQP[0]; k++){
	    V->i[k] = k;
	    V->x[k] = 1.0;
	    V->p[k] = k;
	 }
	 V->p[nQP[0]] = nQP[0];


if(DC[0] == 0.0){
  for(k=0; k<nQP[0]; k++){   

    li[k] = 1.0;                   
    si=0.0;                        

    if(dam[k] == ntwo){
      vii[k] = 1.0;
    } 
      else{
        vii[k] = 0.25*(3 - f[dam[k]]);
      }


    j=k;
    cnt=0;
    while(j>=0){

      dj=dam[j];

      if(dj != ntwo){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      si += li[j]*li[j]*vii[j];

      j=-1;

      for(h=0; h<cnt; h++){   // find eldest individual
       if(AN[h]>j){
         j = AN[h];
       }
      }
      for(h=0; h<cnt; h++){   
        if(AN[h]==j){
          AN[h] -= ntwo; 
        }
      }
    }  // end of while
    f[k] = si - sex[k];

    for(h=0; h<nQP[0]; h++){
      li[h]  = 0.0;            
    }
  } // end of for


} // end hopi
 else{    
  for(k=0; k<nQP[0]; k++){   

    li[k] = 1.0;                   
    si=0.0;                        

    if(sex[k] != 1.0){
      vii[k] = DC[0] * (1.0 - f[dam[k]]);
    }
      else{
        if(sire[k] != ntwo){
          vii[k] = 0.25 * (1.0 - f[dam[k]]);
        }
          else{
            vii[k] = 0.25 * (3.0 - f[dam[k]]);
          }
      }


    j=k;
    cnt=0;
    while(j>=0){

      sj=sire[j];
      dj=dam[j];

      if(sex[j] == 1.0){
        if(sj != ntwo){
          AN[cnt] = sj;
          if(DC[0] == 0.25){
            li[sj] += 1.0*li[j];
          }
            else{
              li[sj] += 0.5*li[j];
            }
          cnt++;
        }
      }

      if(dj != ntwo){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      si += li[j]*li[j]*vii[j];

      j=-1;

      for(h=0; h<cnt; h++){   // find eldest individual
       if(AN[h]>j){
         j = AN[h];
       }
      }
      for(h=0; h<cnt; h++){   
        if(AN[h]==j){
          AN[h] -= ntwo; 
        }
      }
    }  // end of while
    f[k] = si - sex[k];

    for(h=0; h<nQP[0]; h++){
      li[h]  = 0.0;            
    }
  } // end of for


}



  for(k=0; k<nQP[0]; k++){   
      V->x[k] = 1.0 / vii[k];
  }

  QV = cs_multiply(Q, V);
  tS = cs_multiply(QV, tQ);
  S = cs_transpose(tS, TRUE);
  
  for (k = 0 ; k < S->nzmax; k++){
    iSP[k] = S->i[k];
    xSP[k] = S->x[k];
  }
  for (k = 0 ; k <= S->n; k++){
    pSP[k] = S->p[k];
  }
  nzmaxSP[0] = S->nzmax;

  cs_spfree(Q);
  cs_spfree(tQ);
  cs_spfree(V);
  cs_spfree(QV);
  cs_spfree(tS);
  cs_spfree(S);
  delete[] AN;
  delete[] li;
}
}
