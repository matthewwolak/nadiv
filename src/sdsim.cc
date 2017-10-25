#include "nadivcc.h"

extern "C"{  

void sdsim(int *da,  // N dam alleles
	int *sa,     // N sire alleles
        int *eN,     // N (number of replications)
	int *en,     // n pedigree size
	int *dam,    // dam number IDs
	int *sire,   // sire number IDs
        int *sex,    // sex or really number of homogametic sex chromosomes
	int *Sdi,    // i slot of sex-chromosome dominance relatedness matrix
	int *Sdp,    // p slot of matrix
	int *sdij    // x slot of matrix
){         

  int i, j, k, l, m, n, r;
  int mi, si, adij, cdama, csirea, rdama, rsirea;
  int p = 0;
  int c = 0;
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
       }  // end for j
     }

     if(sex[i] == 1){
       if(si != -999){
         k = i*eN[0];
         l = si*eN[0];
         for(j = 0; j < eN[0]; j++){
           sa[k] += da[l];
           k++;
           l++;
         }  // end for j
       }
    }  // end if sex for sire alleles of females (XX)/homogametic sex
  }  // end for i  
  PutRNGstate();


  //     *************      //
  for(m = 0; m < en[0]; m++){
    if(sex[m] == 1){
      Sdp[c] += p;
      c++;
      r = 0;                      // Need a row number of just sex==1
      for(n = 0; n < m+1; n++){
        if(sex[n] == 1){
          adij = 0;
          k = m*eN[0];
          l = n*eN[0];
          for(j = 0; j < eN[0]; j++){
            cdama = da[k];               // 'column' individual's dam allele
            csirea = sa[k];
            rdama = da[l];               // 'row' individual's dam allele
	    rsirea = sa[l];

            if(cdama == rdama){
              if(csirea == rsirea){
                adij += 1;
              }
            } 
            else{
              if(cdama == rsirea){
                if(csirea == rdama){
                  adij += 1;
                }
              }
            }  // end else

            k++;
            l++;
          }  // end for j



          if(adij > 0){
            Sdi[p] += r;
            sdij[p] += adij;
            p++;
          }  // end if adij > 0
          r++; 
        }  // end if sex[n]==1
      }  // end for n
    }  // end if sex[m]==1  
  }  // end for m               
  Sdp[c] += p;
}
}
