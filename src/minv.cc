#include "nadivcc.h"

///////////////////////////////////////////////////////////////////////
//     Mutational Effects (additive) inverse relatedness matrix
////   algorithm in Casellas & Medrano 2008
////   Uses the Quaas 1976 algorithm for A-inverse as a basis
//     Note difference from Wray 1990 algorithm (coded below as `minvw`)
////   Casellas & Medrano make separate invA_0 and invM matrices; Wray combines

extern "C"{  

void minv(
        int *dam,       
        int *sire,         
        double *f,     
	double *dii,
        int *n,
        double *xMinv,
	int *iMinv,
	int *pMinv,
	double *logDet	
){         

  int     i, j, k, m, p, q, sk, dk, lb, step, istart, it;
  double  vi2, detL;
  double  *u = new double[n[0]];
  double  *v = new double[n[0]];
  double  *h = new double[n[0]];
  
  for(k = 0; k < n[0]; k++){    // set u, v, and to zero
    u[k] = 0.0;      
    v[k] = 0.0;      
    h[k] = 0.0;      
  }

  detL = 1.0;  // determinant of LL'=M (log(detL) + log(detL) = log(det(M))
    
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
  
  // TODO default values for u and h so when parent missing yields correct
  // probably make u and h n[0]+1 and assign position n[0]=0 for all but last case
    if(p != n[0] && q != n[0]){  
      v[k] = sqrt(0.25 * (u[p] + u[q]) - 0.5 * (h[p] + h[q]) + 1);
    }
    if(p < n[0] && q == n[0]){
      v[k] = sqrt(0.25*u[p] - 0.5*h[p] + 0.5);
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
          //TODO / FIXME check to see if v below adds to element or replaces?!
          if(p >= k){
            v[j] = 0.5 * (v[p] + v[q]);
            h[j] += 0.5 * v[p] * v[q];
          }
          if(p < k && k <= q) v[j] = 0.5 * v[q];
          if( q < k) v[j] = 0.0;
        }  // end if p and q BOTH known    
        if(p != n[0] && q == n[0]){
          if(k <= p) v[j] = 0.5 * v[p];
          if(p < k) v[j] = 0.0;
        }
        if(p == n[0]) v[j] = 0.0;  // assumes if q>p then if p==n[0] so does q
        
      }  // end if j != k    
      u[j] += v[j] * v[j];

    }  // end for j
    
Rprintf("\n\t k=%i", k);    
    for(i = 0; i < n[0]; i++){
Rprintf("\n u[%i]=%6.3f | v[%i]=%6.3f | h[%i]=%6.3f", i, u[i], i, v[i], i, h[i]);
    }

    ////////////////////////////////////////
    // Now add contributions to M-inverse
    vi2 = v[k] * v[k];
    detL *= v[k];
    
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

  logDet[0] += log(detL) + log(detL);

  delete[] h;
  delete[] v;
  delete[] u;
}
}





///////////////////////////////////////////////////////////////////////
//     Mutational Effects (additive) inverse relatedness matrix
////   algorithm in Wray 1990
////   Uses the Quaas 1976 algorithm for A-inverse as a basis

extern "C"{  

void minvw(
        int *dam,       
        int *sire,         
        double *f,     
	double *dii,
        int *n,
        double *xMinv,
	int *iMinv,
	int *pMinv,
	double *logDet	
){         

  int     i, j, k, m, p, q, sk, dk, lb, step, istart, it;
  double  vi2, detL;
  double  *u = new double[n[0]];
  double  *v = new double[n[0]];
  double  *h = new double[n[0]];
  
  for(k = 0; k < n[0]; k++){    // set u, v, and to zero
    u[k] = 0.0;      
    v[k] = 0.0;      
    h[k] = 0.0;      
  }

  detL = 1.0;  // determinant of LL'=M (log(detL) + log(detL) = log(det(M))
    
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
  
  // TODO default values for u and h so when parent missing yields correct
  // probably make u and h n[0]+1 and assign position n[0]=0 for all but last case
    if(p != n[0] && q != n[0]){  
      v[k] = sqrt(0.25 * (u[p] + u[q]) - 0.5 * (h[p] + h[q]) + 10);
    }
    if(p < n[0] && q == n[0]){
      v[k] = sqrt(0.5 + 0.25*u[p] - 0.5*h[p] + 10);
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
          //TODO / FIXME check to see if v below adds to element or replaces?!
          if(p >= k){
            v[j] = 0.5 * (v[p] + v[q]);
            h[j] += 0.5 * v[p] * v[q];
          }
          if(p < k && k <= q) v[j] = 0.5 * v[q];
          if( q < k) v[j] = 0.0;
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
    detL *= v[k];
    
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

  logDet[0] += log(detL) + log(detL);

  delete[] h;
  delete[] v;
  delete[] u;
}
}






