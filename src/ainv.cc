#include "ddiag.h"

////////////////////////////////

////   H algorithm (nadiv < v2.13.4)

extern "C"{  

void acinv(
        int *dam,       
        int *sire,         
        double *f,     
	int *iTinvP,              
	int *pTinvP,	         
	double *xTinvP,
        int *nTinvP,
	int *nzmaxTinvP,
	int *iAP,
	int *pAP,
	double *xAP,
	int *nzmaxAP
){         

  int     j, k, h, ntwo, cnt, sj, dj;
  double  ai;
  double  *AN = new double[2*nTinvP[0]];
  double  *li = new double[nTinvP[0]];
  cs *Tinv, *D, *tTinv, *tTD;

  for(k=0; k<nTinvP[0]; k++){
     li[k]=0.0;               // set l to zero
  }
  for(k=0; k<nTinvP[0]; k++){
     AN[k]=-1;               // set AN to zero
  }

  ntwo = nTinvP[0];

  Tinv = cs_spalloc(nTinvP[0], nTinvP[0], nzmaxTinvP[0], true, false);  
         for (k = 0 ; k < nzmaxTinvP[0] ; k++){
           Tinv->i[k] = iTinvP[k];
           Tinv->x[k] = xTinvP[k];
         }
         for (k = 0 ; k <= nTinvP[0]; k++){
           Tinv->p[k] = pTinvP[k];
         }
   
  tTinv = cs_transpose(Tinv, true);

  D = cs_spalloc(nTinvP[0], nTinvP[0], nzmaxTinvP[0], true, false);
	 for (k = 0; k < nTinvP[0]; k++){
	    D->i[k] = k;
	    D->x[k] = 1.0;
	    D->p[k] = k;
	 }
	 D->p[nTinvP[0]] = nTinvP[0];
    
  for(k=0; k<nTinvP[0]; k++){  // iterate through each row of l 

    li[k] = 1.0;                   // set l_ii to one
    ai=0.0;                        // set a_ii to zero

    D->x[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);

    j=k;
    cnt=0;

    while(j>=0){

      sj=sire[j];
      dj=dam[j];

      if(sj!= ntwo){
        AN[cnt] = sj;
        li[sj] += 0.5*li[j];
        cnt++;
      }

      if(dj!= ntwo){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      ai += li[j]*li[j]*D->x[j];

      j=-1;

      for(h=0; h<cnt; h++){   // find eldest individual
       if(AN[h]>j){
         j = AN[h];
       }
      }
      for(h=0; h<cnt; h++){   // delete duplicates
        if(AN[h]==j){
          AN[h] -= ntwo; 
        }
      }
    }  // end of while
    f[k] = ai-1.0;
    for(h=0; h<nTinvP[0]; h++){
      li[h]  = 0.0;            // reset l to zero except l_ii =1
    }
  } // end of for


  for(k=0; k<nTinvP[0]; k++){  // iterate through each row of l 
      D->x[k] = sqrt(1.0/(D->x[k]));
  }

  tTD = cs_multiply(tTinv, D);
  
  for (k = 0 ; k < tTD->nzmax; k++){
    iAP[k] = tTD->i[k];
    xAP[k] = tTD->x[k];
  }
  for (k = 0 ; k <= tTD->n; k++){
    pAP[k] = tTD->p[k];
  }
  nzmaxAP[0] = tTD->nzmax;

  cs_spfree(Tinv);
  cs_spfree(tTinv);
  cs_spfree(D);
  cs_spfree(tTD);
  delete[] AN;
  delete[] li;
}
}










/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//   M&L algorithm below
extern "C"{  

void acinvML(
        int *dam,       
        int *sire,         
        double *f, 
        double *dii,    
//	int *iTinvP,              
//	int *pTinvP,	         
//	double *xTinvP,
        int *n,
//	int *nzmaxTinvP,

        double *xA
	int *iA,
	int *pA,
	int *nzmaxA
){         

  int     j, k, h, cnt, sj, dj, tstart;
  double  ai, alphai;
  double  *AN = new double[2*n[0]];
  double  *li = new double[n[0]];
//  cs *Tinv, *D, *tTinv, *tTD;
//  cs *Ainv;

  for(k=0; k<n[0]; k++){
     li[k]=0.0;               // set l to zero
  }
  for(k=0; k<n[0]; k++){
     AN[k]=-1;               // set AN to zero
  }

    
  for(k=0; k<n[0]; k++){  // iterate through each row of l 

    li[k] = 1.0;                   // set l_ii to one
    ai=0.0;                        // set a_ii to zero
    dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);

    j=k;
    cnt=0;

    while(j>=0){

      sj=sire[j];
      dj=dam[j];

      if(sj!= n[0]){
        AN[cnt] = sj;
        li[sj] += 0.5*li[j];
        cnt++;
      }

      if(dj!= n[0]){ 
        AN[cnt] = dj;
        li[dj] += 0.5*li[j];
        cnt++;
      }

      ai += li[j]*li[j]*dii[j];

      j=-1;

      for(h=0; h<cnt; h++){   // find eldest individual
       if(AN[h]>j){
         j = AN[h];
       }
      }
      for(h=0; h<cnt; h++){   // delete duplicates
        if(AN[h]==j){
          AN[h] -= n[0]; 
        }
      }
    }  // end of while
    f[k] = ai-1.0;
// ADDED bit
//FIXME: try sticking her the part with the separate for loop through k<n[0] (N individs)
// END added bit
    for(h=0; h<n[0]; h++){
      li[h]  = 0.0;            // reset l to zero except l_ii =1
    }
  } // end of for



//  Ainv = cs_spalloc(n[0], n[0], nzmaxA[0], true, false);
  cnt = 0;
  for(k=0; k<n[0]; k++){  // iterate through each individual 

// LEFT OFF HERE
      alphai = sqrt(1.0/(dii[k]));
      tstart = k*6;
      // i,i
      xA[tstart] += alphai;
      if(sire[j] != n[0]){
         // i,sire
         xA[tstart+1] += alphai * -0.5;
         // sire,sire
         xA[tstart+3] += alphai * 0.25;
      }
      if(dam[j] != n[0]){
         // i,dam
         xA[tstart+2] += alphai * -0.5;
         // dam,dam
         xA[tstart+5] += alphai * 0.25;
         if(sire[j] != n[0]){
            //sire,dam/dam,sire
            xA[tstart+4] += alphai * 0.25;
         }
      }
  }

  
//  for (k = 0 ; k < tTD->nzmax; k++){
//    iAP[k] = tTD->i[k];
//    xAP[k] = tTD->x[k];
//  }
//  for (k = 0 ; k <= tTD->n; k++){
//    pAP[k] = tTD->p[k];
//  }
//  nzmaxAP[0] = tTD->nzmax;

  delete[] AN;
  delete[] li;
}
}
