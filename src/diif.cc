#include "nadivcc.h"

// R-interface/wrapper for M&L routine
extern "C"{  

void diif(
	int *dam,
	int *sire,
	double *f,
	double *dii,
	int *n,
	int *g,
	int *fmiss
){

  // Meuwissen and Luo 1992 algorithm to obtain f and dii values
  ml(dam, sire, f, dii, n[0], g[0], fmiss[0]);
}
}  






// R-interface/wrapper for mutational effects M&L routine
extern "C"{  

void mdiif(
	int *dam,
	int *sire,
	double *h,
	double *dii,
	int *n
){

  // Meuwissen and Luo 1992 algorithm to obtain f and dii values
  //// Extends Wray 1990; Casellas and Medrano 2008
  mml(dam, sire, h, dii, n[0]);
}
}  






////////////////////////////////////////////////
//   based on M&L 1992 algorithm (for `ainvml`)
//   as presented in Mrode 2005
//// replaces elements of f and dii with calculated values in place
void ml(int *dam, int *sire,
	double *f, double *dii,
	int n, int g, int fmiss
){

  int     j, k, h, cnt, sj, dj;
  double  ai;
  double  *AN = new double[2*n];
  double  *li = new double[n];

  for(k = g; k < n; ++k){
     li[k] = 0.0;               // set l to zero
  }
  for(k = g; k < n; ++k){
     AN[k] = -1;      // set AN to "empty" 
                      //// (since lowest ID is 0, make empty with 1 less than it)
  }

  for(k = g; k < n; ++k){  // iterate through each row of l 
    dii[k] = 0.5 - 0.25*(f[dam[k]] + f[sire[k]]);
    
    if(fmiss == 1){  // only do below if f coefficients NOT supplied by user
      if((k > 0) && (dam[k] == dam[k-1]) && (sire[k] == sire[k-1])){
        f[k] += f[k-1];
      } 
      else {
        li[k] = 1.0;                   // set l_ii to one
        ai = 0.0;                      // set a_ii to zero
        j = k;
        cnt = 0;
        while(j >= 0){
          sj = sire[j];
          dj = dam[j];

          if((sj >= g) && (sj != n)){
            AN[cnt] = sj;
            li[sj] += 0.5*li[j];
            ++cnt;
          }

          if((dj >= g) && (dj != n)){
            AN[cnt] = dj;
            li[dj] += 0.5*li[j];
            ++cnt;
          }

          ai += li[j]*li[j]*dii[j];
          j -= n;          // set to empty-value lower than all known identities

          for(h = 0; h < cnt; ++h){   // find eldest individual
            if(AN[h] > j){
              j = AN[h];
            }
          }
          for(h = 0; h < cnt; ++h){   // delete duplicates
            if(AN[h] == j){
              AN[h] -= n;    // set to empty-value lower than all known identities
            }
          }
        }  // end of while
        
        f[k] = ai - 1.0;

        for(h = 0; h <= k; ++h){
          li[h] = 0.0;            // reset l to zero
        }

      } // end else for checking if k has same parents as k-1
    }  // end if f missing
  } // end of for
  delete[] AN;
  delete[] li;

}
  














// Mutational effects inbreeding and dii
//// based on Meuwissen and Luo 1992 algorithm to obtain f and dii values
//// Extends Wray 1990; Casellas and Medrano 2008

void mml(int *dam, int *sire,
	double *h, double *dii,
	int n
){

  int     j, k, m, p, q, cnt, sj, dj;
  double  g;
  double  *AN = new double[2*n];
  double  *li = new double[n];
  double  *u = new double[n];
  
  for(k = 0; k < n; ++k){
     li[k] = 0.0;     // set l to zero
     AN[k] = -1;      // set AN to "empty" 
                      //// (since lowest ID is 0, make empty with 1 less than it)
     u[k] = 0.0;      // set u to zero
  }


  for(k = 0; k < n; ++k){  // iterate through each row of L 

    // intialize/guess and change if different
    //// if sire < dam OR sire=dam=UNKNOWN
    p = sire[k];
    q = dam[k];
    //// otherwise, change so p is always < q (unless p=q)
    if(sire[k] > dam[k]){
      p = dam[k];
      q = sire[k];
    }
  
    if(p != n && q != n){
      dii[k] = 0.25 * (u[p] + u[q]) - 0.5 * (h[p] + h[q]) + 1.0;
    }
    if(p < n && q == n) dii[k] = 0.25*u[p] - 0.5*h[p] + 0.5;  
    // because p <= q then if p=n=missing ID, THEN so will q
    if(p == n) dii[k] = 1.0; 

    
    li[k] = 1.0;     // set L_ii to one
    j = k;
    cnt = 0;
    g = 0.0;
    while(j >= 0){
      sj = sire[j];
      dj = dam[j];
      
      if(sj != n){
        AN[cnt] = sj;
        li[sj] += 0.5 * li[j];
        ++cnt;
      }

      if(dj != n){
        AN[cnt] = dj;
        li[dj] += 0.5 * li[j];
        ++cnt;
      }

      u[k] += li[j] * li[j] * dii[j];
      g += li[j];
      
      j -= n;          // set to empty, value lower than all known identities
      
      for(m = 0; m < cnt; ++m){    // find the eldest individual
        if(AN[m] > j){
          j = AN[m];
        }
      }
      for(m = 0; m < cnt; ++m){    // delete duplicates
        if(AN[m] == j){
          AN[m] -= n;  // set to negative value so never `AN[m]>j` in above
        }
      }
      
    }  // end of while

    h[k] = u[k] - g;

    for(m = 0; m <= k; ++m) li[m] = 0.0;    // reset li to zero

  } // end of for k
  
  delete[] AN;
  delete[] li;
  delete[] u;

}

