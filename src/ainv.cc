#include "nadivcc.h"


//////////////////////////////////////
//   M&L 1992 algorithm 
//   as presented in Mrode 2005
// since nadiv >v2.15.0 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 
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

  int     lb, step, it, k, h, sj, dj, istart;
  double  alphai;


  // Meuwissen and Luo 1992 algorithm to obtain f and dii values
  //TODO FIXME: last argument assumes f NOT supplied by user
  ml(dam, sire, f, dii, n[0], g[0], 1);

  for(k = g[0]; k < n[0]; ++k){  // iterate through each individual (not groups)
    // check to see if k has 2 phantom parents from same genetic group
    sj = sire[k];
    dj = dam[k];
    if((sj == dj) && (sj < g[0])){
      xA[pA[k]] += 1.0;       // k,k
      istart = pA[sj];
      xA[istart] += 1.0;      // sire,sire (same as dam,dam)
      h = istart;
      lb = pA[sj + 1] - 1 - h;
      while(lb > 0){
        step = lb/2;
        it = h + step;
        if(iA[it] < k){
          h = ++it;
          lb -= step + 1;
        }
        else lb = step;
      }  // end while
      if(iA[h] == k) xA[h] += -1.0;      //sire/dam,k
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
              h = istart;
              lb = pA[sj+1] - 1 - h;
              while(lb > 0){
                step = lb/2;
                it = h + step;
                if(iA[it] < dj){
                  h = ++it;
                  lb -= step+1;
                }
                else lb = step;
              }
              if(iA[h] == dj) xA[h] += alphai;
            }
         }  // end if sj <= dj
         // sire,k
         h = istart;
         lb = pA[sj+1] - 1 - h;
         while(lb > 0){
           step = lb/2;
           it = h + step;
           if(iA[it] < k){
             h = ++it;
             lb -= step+1;
           }
           else lb = step;
         }
         if(iA[h] == k) xA[h] += alphai * -2.0;
      } // end if sj not missing
      if(dj != n[0]){
         istart = pA[dj];
         // dam,dam
         xA[istart] += alphai;
         // sire
         // dam,k
         h = istart;
         lb = pA[dj+1] - 1 - h;
         while(lb > 0){
           step = lb/2;
           it = h + step;
           if(iA[it] < k){
             h = ++it;
             lb -= step+1;
           }
           else lb = step;
         } // end while
         if(iA[h] == k) xA[h] += alphai * -2.0;
         // dam,sire
         if(dj <= sj){
            if(sj != n[0]){
              h = istart;
              lb = pA[dj+1] - 1 - h;
              while(lb > 0){
                step = lb/2;
                it = h + step;
                if(iA[it] < sj){
                  h = ++it;
                  lb -= step+1;
                }
                else lb = step;
              } // end while
              if(iA[h] == sj) xA[h] += alphai;
            } // end if sj not missing
         } // end if dj <= sj
       } // end if dj not missing
    } // end else k doesn't have 2 phantom parents from same group
  } // end for k

}
}








///////////////////////////////////////////////////////////////////////////
//   basically the same M&L 1992 algorithm as 'ainvml', but
//   with fuzzy classification of genetic groups
// since nadiv > v2.14.3 uses lower_bound algorithm for matrix lookup
//// based on c++ <algorithm>std::lower_bound 
extern "C"{  

void ainvfuzz(
        int *dam,       
        int *sire,         
        int *phdam,       
        int *phsire,         
        double *f, 
        double *dii,    
        int *n,
	int *g,
	double *xF,
	int *iF,
	int *pF,
        double *xA,
	int *iA,
	int *pA
){         

  int     lb, step, it, h, i, j, k, s, d, cnt, sj, dj, mj, mk, sp, dp, fistart, aistart;
  double  ai, alphai, pij, pijp, pik, pikp;
  double  *AN = new double[2*n[0]];
  double  *li = new double[n[0]];

  for(k = g[0]; k < n[0]; k++){
     li[k]=0.0;               // set l to zero
  }
  for(k = g[0]; k < n[0]; k++){
     AN[k] = -1;               // set AN to "zero" 
                               //// (since an ID is 0, make 1 less than lowest ID)
  }

  for(k = g[0]; k<n[0]; k++){  // iterate through each row of l 
    dii[k] = 0.5-0.25*(f[dam[k]]+f[sire[k]]);
    if((k > 0) && (phdam[k] == phdam[k-1]) && (phsire[k] == phsire[k-1])){
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

    ///////////////////////////////////
    // number of phantom parents of k
    if(sire[k] == n[0]){
      mj=1;
      sj=phsire[k];
    }
    else {
      mj=0;
      sj=sire[k];
    } 
    if(dam[k] == n[0]){
      mk=1;
      dj=phdam[k];
    }
    else {
      mk=0;
      dj=dam[k];
    } 

    alphai = 1.0/(dii[k] * 4.0);
    // bottom-left triangle of matrix (with diagonal included)
    //// always use the column of the lesser ID and then find the row of the greater ID
    aistart = pA[k];
    // k,k
    xA[aistart] += alphai * 4.0;

    //////////////////////////////////////////////////////
    /// 1 of 4: both parents are phantom parents
    //////////////////////////////////////////////////////
    if((mj == 1) && (mk == 1)){
      ////// the sire contributions to k
      for(s=0; s<g[0]; s++){                   // for each sire genetic group
        fistart = pF[s];
        pij=0.0;
        for(j=fistart; j<pF[s+1]; j++){
          if(iF[j] == sj){
            pij += xF[j];
            // k,sire-group s
            aistart = pA[s];
            h = aistart;
            lb = pA[s+1] - 1 - h;
            while(lb > 0){
              step = lb/2;
              it = h + step;
              if(iA[it] < k){
                h=++it;
                lb-=step+1;
              }
              else lb = step;
            }
            if(iA[h] == k) xA[h] += alphai * -2.0 * pij;

            for(sp=0; sp<=s; sp++){
              fistart = pF[sp];
              pijp=0.0; 
              for(h=fistart; h<pF[sp+1]; h++){
                if(iF[h] == sj){
                  pijp += xF[h];
                  // group, k's sire's-group sp
                  aistart = pA[sp];
                  i = aistart;
                  lb = pA[sp+1] - 1 - i;
                  while(lb > 0){
                    step = lb/2;
                    it = i + step;
                    if(iA[it] < s){
                      i=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iA[i] == s) xA[i] += alphai * pij * pijp;
                  break;      			// break out of h for loop
                }  // end if(iF[h] == sj)
              }  // end h for loop

            }  // end sp for loop
            break;      			// break out of j for loop
          }  // end if(iF[j] == sj)
        }  // end j for loop

      }  // end s for loop

      /////////////////////////////////
      ////// the dam contributions to k
      for(d=0; d<g[0]; d++){                   // for each dam genetic group
        fistart = pF[d];
        pik=0.0;
        for(j=fistart; j<pF[d+1]; j++){
          if(iF[j] == dj){
            pik += xF[j];
            // k,dam-group d
            aistart = pA[d];
            h  =aistart;
            lb = pA[d+1] - 1 - h;
            while(lb > 0){
              step = lb/2;
              it = h + step;
              if(iA[it] < k){
                h=++it;
                lb-=step+1;
              }
              else lb  =step;
            }
            if(iA[h] == k) xA[h] += alphai * -2.0 * pik;

            for(dp=0; dp<=d; dp++){
              fistart = pF[dp];
              pikp=0.0; 
              for(h=fistart; h<pF[dp+1]; h++){
                if(iF[h] == dj){
                  pikp += xF[h];
                  // group, k's dam's-group dp
                  aistart = pA[dp];
                  i  =aistart;
                  lb = pA[dp+1] - 1 - i;
                  while(lb > 0){
                    step = lb/2;
                    it = i + step;
                    if(iA[it] < d){
                      i=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iA[i] == d) xA[i] += alphai * pik * pikp;
                  break;        		// break out of h for loop
                }
              }  // end h for loop

            }  // end dp for loop
            break;				// break out of j for loop
          }
        }  // end j for loop

      }  // end d for loop

      ////////////////////////////////////////
      ////// the sire contributions to the dam
      for(s=0; s<g[0]; s++){                   // for each sire genetic group
        fistart = pF[s];
        pij=0.0;
        for(j=fistart; j<pF[s+1]; j++){
          if(iF[j] == sj){
            pij += xF[j];
            for(d=0; d<g[0]; d++){
              fistart = pF[d];
              pik=0.0;
              for(h=fistart; h<pF[d+1]; h++){
                if(iF[h] == dj){
                  pik += xF[h];
                  // k's sire's group, k's dam's group 
                  aistart = pA[min(s, d)];   
                  if(s != d){
                    i = aistart;
                    lb = pA[min(s, d)+1] - 1 - i;
                    while(lb > 0){
                      step = lb/2;
                      it = i + step;
                      if(iA[it] < max(s, d)){
                        i=++it;
                        lb-=step+1;
                      }
                      else lb = step;
                    }
                    if(iA[i] == max(s, d)) xA[i] += alphai * pij * pik;
                  }
                  else {
                    xA[aistart] += alphai * 2.0 * pij * pik;
                  }
                  break;  			// break out of h for loop
                }
              }  // end h for loop

            }  // end d for loop
            break;  				// break out of j for loop
          }
        }  // end j for loop

      }  // end s for loop
    }  // end if both parents phantoms





    //////////////////////////////////////////////////////
    /// 2 of 4: phantom sire, but known dam
    //////////////////////////////////////////////////////
    if((mj == 1) && (mk == 0)){
      ////// the sire contributions to k
      for(s=0; s<g[0]; s++){                   // for each sire genetic group
        fistart = pF[s];
        pij=0.0;
        for(j=fistart; j<pF[s+1]; j++){
          if(iF[j] == sj){
            pij += xF[j];
            aistart = pA[s];
            for(h=aistart; h<pA[s+1]; h++){  // NOTE: not using lower_bound algorithm
              // dam of k,sire-group s
              if(iA[h] == dj){
                xA[h] += alphai * pij;
              }
              // k,sire-group s
              if(iA[h] == k){
                xA[h] += alphai * -2.0 * pij;
                break;  			// break out of h for loop
              }
            }  // end h for loop

            for(sp=0; sp<=s; sp++){
              fistart = pF[sp];
              pijp=0.0; 
              for(h=fistart; h<pF[sp+1]; h++){
                if(iF[h] == sj){
                  pijp += xF[h];
                  // group, k's sire's-group sp
                  aistart = pA[sp];
                  i = aistart;
                  lb = pA[sp+1] - 1 - i;
                  while(lb > 0){
                    step = lb/2;
                    it = i + step;
                    if(iA[it] < s){
                      i=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iA[i] == s) xA[i] += alphai * pij * pijp;
                  break;  			// break out of h for loop
                }
              }  // end h for loop

            }  // end sp for loop
            break;  				// break out of j for loop
          }
        }  // end j for loop

      }  // end s for loop

      //////////////////////////////////
      ////// the dam contributions to k
      aistart = pA[dj];
      xA[aistart] += alphai;
      j = aistart;
      lb = pA[dj+1] - 1 - j;
      while(lb > 0){
        step = lb/2;
        it = j + step;
        if(iA[it] < k){
          j=++it;
          lb-=step+1;
        }
        else lb = step;
      }
      if(iA[j] == k) xA[j] += alphai * -2.0;
    }  // end if phantom sire and known dam






    //////////////////////////////////////////////////////
    /// 3 of 4: known sire, phantom dam
    //////////////////////////////////////////////////////
    if((mj == 0) && (mk == 1)){
      ////// the sire contributions to k
      aistart = pA[sj];
      xA[aistart] += alphai;
      j = aistart;
      lb = pA[sj+1] - 1 - j;
      while(lb > 0){
        step = lb/2;
        it = j + step;
        if(iA[it] < k){
          j=++it;
          lb-=step+1;
        }
        else lb = step;
      }
      if(iA[j] == k) xA[j] += alphai * -2.0;
      //////////////////////////////////
      ////// the dam contributions to k
      for(d=0; d<g[0]; d++){                   // for each dam genetic group
        fistart = pF[d];
        pik=0.0;
        for(j=fistart; j<pF[d+1]; j++){
          if(iF[j] == dj){
            pik += xF[j];
            aistart = pA[d];
            for(h=aistart; h<pA[d+1]; h++){  // NOTE: not using lower_bound algorithm
              // sire of k,dam-group d
              if(iA[h] == sj){
                xA[h] += alphai * pik;
              }
              // k,dam-group d
              if(iA[h] == k){
                xA[h] += alphai * -2.0 * pik;
                break;  			// break out of h for loop
              }
            }  // end h for loop
            for(dp=0; dp<=d; dp++){
              fistart = pF[dp];
              pikp=0.0; 
              for(h=fistart; h<pF[dp+1]; h++){
                if(iF[h] == dj){
                  pikp += xF[h];
                  // k's dam, k's dam's-group dp
                  aistart = pA[dp];
                  i = aistart;
                  lb = pA[dp+1] - 1 - i;
                  while(lb > 0){
                    step = lb/2;
                    it= i + step;
                    if(iA[it] < d){
                      i=++it;
                      lb-=step+1;
                    }
                    else lb = step;
                  }
                  if(iA[i] == d) xA[i] += alphai * pik * pikp;
                  break;  			// break out of h for loop
                }
              }  // end h for loop

            }  // end dp for loop
            break;  				// break out of j for loop
          }
        }  // end j for loop

      }  // end d for loop
    }  // end if known sire, phantom dam






    //////////////////////////////////////////////////////
    /// 4 of 4: known sire & dam
    //////////////////////////////////////////////////////
    if((mj == 0) && (mk == 0)){
      aistart = pA[sj];
      // sire,sire
      xA[aistart] += alphai;
      for(j=aistart; j<pA[sj+1]; j++){  // NOTE: not using lower_bound algorithm
        // dam,sire
        if(iA[j] == dj){
          xA[j] += alphai;
        }        
        // k,sire
        if(iA[j] == k){
          xA[j] += alphai * -2.0;
          break;  				// break out of j for loop
        }
      }  // end j for loop

      aistart = pA[dj];
      // dam,dam
      xA[aistart] += alphai;
      for(j=aistart; j<pA[dj+1]; j++){  // NOTE: not using lower_bound algorithm
        //sire,dam
        if(iA[j] == sj){
          xA[j] += alphai;
        }
        // k,dam
        if(iA[j] == k){
          xA[j] += alphai * -2.0;
          break;  				// break out of j for loop
        }
      }  // end j for loop

    }  // end if known sire & dam
  } // end k for loop
  delete[] AN;
  delete[] li;
}
}


