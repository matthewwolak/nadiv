#include "nadivcc.h"

///////////////////////////////////////////////////////////////////////
//     Mutational Effects (additive) inverse relatedness matrix
////   algorithm in Casellas & Medrano 2008
////   Uses the Quaas 1976 algorithm for A-inverse as a basis
//     Note difference from Wray 1990 algorithm (coded below as `minvw`)
////   Casellas & Medrano make separate invA_0 and invM matrices; Wray combines

extern "C"{  

void minvq(
        int *dam,       
        int *sire,         
        double *h,     // essentially coeff. of inbreeding (f)     
	double *v,     // ultimately `dii`, but starts as sqrt(dii)=v
        int *n,
        double *xMinv,
	int *iMinv,
	int *pMinv,
	double *logDet,
	double *theta	
){         

  int     j, k, m, p, q, sk, dk, lb, step, istart, it;
  double  vi2, detM;
  double  *u = new double[n[0]];
  
  for(k = 0; k < n[0]; k++) u[k] = 0.0;      // initialize u to zero

  detM = 1.0;  // determinant of M=TDT'=prod(diag(D))
    
  for(k = 0; k < n[0]; k++){  // iterate through each individual k 

    // intialize/guess and change if different
    //// if sire < dam OR sire=dam=UNKNOWN
    p = sire[k];
    q = dam[k];
    //// otherwise, change so p is always < q (unless p=q)
    if(sire[k] > dam[k]){
      p = dam[k];
      q = sire[k];
    }
  
    if(p != n[0] && q != n[0]){  
      v[k] = sqrt(0.25 * (u[p] + u[q]) - 0.5 * (h[p] + h[q]) + theta[0]);
    }
    if(p < n[0] && q == n[0]){
      if(theta[0] == 1.0){
        v[k] = sqrt(0.25*u[p] - 0.5*h[p] + 0.5);
      } else{
          v[k] = sqrt(0.5 + 0.25*u[p] - 0.5*h[p] + theta[0]);
        }
    }
    if(p == n[0]) v[k] = 1.0; // because p <= q then if p=n[0]=missing ID THEN so will q


    for(j = k; j < n[0]; j++){
      if(j > k){
        // intialize/guess and change if different
        //// if sire < dam OR sire=dam=UNKNOWN
        p = sire[j];
        q = dam[j];
        //// otherwise, change so p is always < q (unless p=q)
        if(sire[j] > dam[j]){
          p = dam[j];
          q = sire[j];
        }

        if(p != n[0] && q != n[0]){
          if(p >= k){
            v[j] = 0.5 * (v[p] + v[q]);
            h[j] += 0.5 * v[p] * v[q];
          }
          if(p < k && k <= q) v[j] = 0.5 * v[q];
          if(q < k) v[j] = 0.0;
        }  // end if p and q BOTH known    
        if(p != n[0] && q == n[0]){
          if(k <= p) v[j] = 0.5 * v[p];
          if(p < k) v[j] = 0.0;
        }
        if(p == n[0]) v[j] = 0.0;  // assumes if q>p then if p==n[0] so does q
        
      }  // end if j != k    
      u[j] += v[j] * v[j];

    }  // end for j


    ////////////////////////////////////////
    // Now add contributions to M-inverse
    vi2 = v[k] * v[k];
    detM *= vi2;
    
    sk = sire[k];
    dk = dam[k];
    istart = pMinv[k];
    // k,k
    xMinv[istart] += 1.0 / vi2;    
    if(sk != n[0]){
     istart = pMinv[sk];
     // sire,sire
     xMinv[istart] += 1.0 / (4.0 * vi2);
     // sire,dam
     if(sk <= dk){
        if(dk != n[0]){
          m = istart;
          lb = pMinv[sk+1] - 1 - m;
              while(lb > 0){
                step = lb/2;
                it = m + step;
                if(iMinv[it] < dk){
                  m = ++it;
                  lb -= step+1;
                }
                else lb = step;
              }
              if(iMinv[m] == dk) xMinv[m] += 1.0 / (4.0 * vi2);
            }
         }
         // sire,k
         m = istart;
         lb = pMinv[sk+1] - 1 - m;
         while(lb > 0){
           step = lb/2;
           it = m + step;
           if(iMinv[it] < k){
             m = ++it;
             lb -= step+1;
           }
           else lb = step;
         }
         if(iMinv[m] == k) xMinv[m] += -1.0 / (2.0 * vi2);
      }  // end if sire KNOWN
      
      if(dk != n[0]){
         istart = pMinv[dk];
         // dam,dam
         xMinv[istart] += 1.0 / (4.0 * vi2);
         // dam,k
         m = istart;
         lb = pMinv[dk+1] - 1 - m;
         while(lb > 0){
           step = lb/2;
           it = m + step;
           if(iMinv[it] < k){
             m = ++it;
             lb -= step+1;
           }
           else lb = step;
         }
         if(iMinv[m] == k) xMinv[m] += -1.0 / (2.0 * vi2);
         // dam,sire
         if(dk <= sk){
            if(sk != n[0]){
              m = istart;
              lb = pMinv[dk+1] - 1 - m;
              while(lb > 0){
                step = lb/2;
                it = m + step;
                if(iMinv[it] < sk){
                  m = ++it;
                  lb -= step+1;
                }
                else lb = step;
              }
              if(iMinv[m] == sk) xMinv[m] += 1.0 / (4.0 * vi2);
            }
         }
       }  // end if dam KNOWN

   
  }  // end for k


  // now replace v with v-squared to be `dii` to return to R
  for(k = 0; k < n[0]; k++) v[k] *= v[k];
  // Calculate log determinant of M (needed for REML)    
  logDet[0] += log(detM);

  delete[] u;
}
}



///////////////////////////////////////////////////////////////////////









///////////////////////////////////////////////////////////////////////
//     Mutational Effects (additive) inverse relatedness matrix
////   algorithm in Casellas & Medrano 2008
////   Uses the Quaas 1976 algorithm for A-inverse as a basis
//     Note difference from Wray 1990 algorithm (coded below as `minvw`)
////   Casellas & Medrano make separate invA_0 and invM matrices; Wray combines

extern "C"{  

void minvml(
        int *dam,       
        int *sire,         
        double *h,     // essentially coeff. of inbreeding (f)     
	double *dii,
        int *n,
        double *xMinv,
	int *iMinv,
	int *pMinv,
	double *logDet
){         

  int     j, k, m, p, q, cnt, sk, dk, sj, dj, lb, step, istart, it;
  double  g, alphai, detM;
  double  *AN = new double[2*n[0]];
  double  *li = new double[n[0]];
  double  *u = new double[n[0]];
  
  for(k = 0; k < n[0]; k++){
     li[k]=0.0;               // set l to zero
     AN[k]=-1;                // set AN to zero
     u[k] = 0.0;              // set u to zero
  }

  detM = 1.0;  // determinant of M=TDT'=prod(diag(D))
    
  for(k = 0; k < n[0]; k++){  // iterate through each row of L 

    // intialize/guess and change if different
    //// if sire < dam OR sire=dam=UNKNOWN
    p = sire[k];
    q = dam[k];
    //// otherwise, change so p is always < q (unless p=q)
    if(sire[k] > dam[k]){
      p = dam[k];
      q = sire[k];
    }
  
    if(p != n[0] && q != n[0]){
      dii[k] = 0.25 * (u[p] + u[q]) - 0.5 * (h[p] + h[q]) + 1.0;
    }
    if(p < n[0] && q == n[0]) dii[k] = 0.25*u[p] - 0.5*h[p] + 0.5;  
    if(p == n[0]) dii[k] = 1.0; // because p <= q then if p=n[0]=missing ID THEN so will q

    alphai = 1 / dii[k];
    detM *= dii[k];    // add contributions to determinant of M
    
    li[k] = 1.0;     // set L_ii to one
    j = k;
    cnt = 0;
    g = 0.0;
    while(j >= 0){
      sj = sire[j];
      dj = dam[j];
      
      if(sj != n[0]){
        AN[cnt] = sj;
        li[sj] += 0.5 * li[j];
        cnt++;
      }

      if(dj != n[0]){
        AN[cnt] = dj;
        li[dj] += 0.5 * li[j];
        cnt++;
      }

      u[k] += li[j] * li[j] * dii[j];
      g += li[j];
      
      j =- 1;
      
      for(m = 0; m < cnt; m++){    // find the eldest individual
        if(AN[m] > j){
          j = AN[m];
        }
      }
      for(m = 0; m < cnt; m++){    // delete duplicates
        if(AN[m] == j){
          AN[m] -= n[0];  // set to negative value so never `AN[m]>j` in above
        }
      }
      
    }  // end of while

    h[k] = u[k] - g;

    for(m = 0; m <= k; m++) li[m] = 0.0;    // reset li to zero
            

    ////////////////////////////////////////
    // Now add contributions to M-inverse
    sk = sire[k];
    dk = dam[k];
    istart = pMinv[k];
    // k,k
    xMinv[istart] += alphai;    
    if(sk != n[0]){
     istart = pMinv[sk];
     // sire,sire
     xMinv[istart] += 0.25 * alphai;
     // sire,dam
     if(sk <= dk){
        if(dk != n[0]){
          m = istart;
          lb = pMinv[sk+1] - 1 - m;
              while(lb > 0){
                step = lb/2;
                it = m + step;
                if(iMinv[it] < dk){
                  m = ++it;
                  lb -= step+1;
                }
                else lb = step;
              }
              if(iMinv[m] == dk) xMinv[m] += 0.25 * alphai;
            }
         }
         // sire,k
         m = istart;
         lb = pMinv[sk+1] - 1 - m;
         while(lb > 0){
           step = lb/2;
           it = m + step;
           if(iMinv[it] < k){
             m = ++it;
             lb -= step+1;
           }
           else lb = step;
         }
         if(iMinv[m] == k) xMinv[m] += -0.5 * alphai;
      }  // end if sire KNOWN
      
      if(dk != n[0]){
         istart = pMinv[dk];
         // dam,dam
         xMinv[istart] += 0.25 * alphai;
         // dam,k
         m = istart;
         lb = pMinv[dk+1] - 1 - m;
         while(lb > 0){
           step = lb/2;
           it = m + step;
           if(iMinv[it] < k){
             m = ++it;
             lb -= step+1;
           }
           else lb = step;
         }
         if(iMinv[m] == k) xMinv[m] += -0.5 * alphai;
         // dam,sire
         if(dk <= sk){
            if(sk != n[0]){
              m = istart;
              lb = pMinv[dk+1] - 1 - m;
              while(lb > 0){
                step = lb/2;
                it = m + step;
                if(iMinv[it] < sk){
                  m = ++it;
                  lb -= step+1;
                }
                else lb = step;
              }
              if(iMinv[m] == sk) xMinv[m] += 0.25 * alphai;
            }
         }
       }  // end if dam KNOWN

   
  }  // end for k


  // Calculate log determinant of M (needed for REML)    
  logDet[0] += log(detM);

  delete[] u;
  delete[] li;
  delete[] AN;
  
}
}




