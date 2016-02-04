#include "nadivcc.h"

extern "C"{  

void reT(
	int *dam,
	int *sire,
	int *i,
	int *p,
	double *x,
	int *maxcnt,
	int *n,
	double *tx
	
){         

  int     k, kdam, ksire, cnt;

  cnt = 0;
  for(k = 0; k < n[0]; k++){
     p[k] = cnt;
     kdam = dam[k];
     ksire = sire[k];
     if(kdam == ksire){
        if(kdam != -999){
           i[cnt] += kdam;
           x[cnt] -= tx[2];
           cnt++;
        }
     }
        else{
           if(kdam < ksire){
              if(kdam != -999){
                 i[cnt] += kdam;
                 x[cnt] -= tx[0];
                 cnt++;
              }
              if(ksire != -999){
                 i[cnt] += ksire;
                 x[cnt] -= tx[1];
                 cnt++;
              }
           }
              else{
                 if(ksire != -999){
                    i[cnt] += ksire;
                    x[cnt] -= tx[1];
                    cnt++;
                 }
                 if(kdam != -999){
                    i[cnt] += kdam;
                    x[cnt] -= tx[0];
                    cnt++;
                 }
              }
        }
     i[cnt] += k;
     x[cnt] += tx[3];
     cnt++;
  }
  p[n[0]] += cnt;
  maxcnt[0] = cnt;
        
}
}
