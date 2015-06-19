#include "ainv.h"

/////////////////////////
//   M&L 1992 algorithm 
extern "C"{  

void ainvml(
        int *dam,       
        int *sire,         
        double *f, 
        double *dii,    
        int *n,
	int *g,
        double *xA,
	int *iA,
	int *pA,
	int *nzmaxA
){         

  int     j, k, h, cnt, sj, dj, istart;
  double  ai, alphai;
  double  *AN = new double[2*n[0]];
  double  *li = new double[n[0]];

  for(k=g[0]; k<n[0]; k++){
     li[k]=0.0;               // set l to zero
  }
  for(k=g[0]; k<n[0]; k++){
     AN[k]=-1;               // set AN to zero
  }

  for(k=g[0]; k<n[0]; k++){  // iterate through each row of l 
    dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);
    if((k > 0) && (dam[k] == dam[k-1]) && (sire[k] == sire[k-1])){
      f[k] += f[k-1];
    } 
    else {
      li[k] = 1.0;                   // set l_ii to one
      ai=0.0;                        // set a_ii to zero
      j=k;
      cnt=0;
      while(j>=0){
        sj=sire[j];
        dj=dam[j];

        if((sj >= g[0]) && (sj!= n[0])){
          AN[cnt] = sj;
          li[sj] += 0.5*li[j];
          cnt++;
        }

        if((dj >= g[0]) && (dj!= n[0])){ 
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
      for(h=0; h<=k; h++){
        li[h]  = 0.0;            // reset l to zero except l_ii =1
      }

    } // end else for checking if k has same parents as k-1

    // check to see if k has 2 phantom parents from same genetic group
    sj = sire[k];
    dj = dam[k];
    if((sj == dj) && (sj < g[0])){
      xA[pA[k]] += 1.0;       // k,k
      istart = pA[sj];
      xA[istart] += 1.0;      // sire,sire (same as dam,dam)
      for(j=istart; j<pA[sj+1]; j++){
        if(iA[j] == k){
          xA[j] += -1.0;      //sire/dam,k
          break;
        }
      }
    }
    // if k doesn't have two phantom parents from same genetic group
    else {
      alphai = 1.0/(dii[k] * 4.0);
      istart = pA[k];
      // k,k
      xA[istart] += alphai * 4.0;
      if(sj != n[0]){
         istart = pA[sj];
         // sire,sire
         xA[istart] += alphai;
         // sire,dam
         if(sj <= dj){
            if(dj != n[0]){
               for(j=istart; j<pA[sj+1]; j++){
                 if(iA[j] == dj){
                   xA[j] += alphai; 
                   break;
                 }
               }
            }
         }
         // sire,k
         for(j=istart; j<pA[sj+1]; j++){
            if(iA[j] == k){
              xA[j] += alphai * -2.0;
              break;
            }
         }
      }
      if(dj != n[0]){
         istart = pA[dj];
         // dam,dam
         xA[istart] += alphai;
         for(j=istart; j<pA[dj+1]; j++){
         // sire
         // dam,k
            if(iA[j] == k){
              xA[j] += alphai * -2.0;
            }
         }
         // dam,sire
         if(dj <= sj){
            if(sj != n[0]){
               for(j=istart; j<pA[dj+1]; j++){
                 if(iA[j] == sj){
                   xA[j] += alphai;
                   break;
                 }
               }
            }
         }
       }
    }
  } // end of for
  delete[] AN;
  delete[] li;
}
}



///////////////////////////////////////////////////////////////////////
////   H implementation of Meuwissen & Luo algorithm (nadiv < v2.13.4)

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

