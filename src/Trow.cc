#include "nadivcc.h"

extern "C"{  

void Trow(
	int     *dam,
	int     *sire,
	double  *x,
	int     *i,
	int     *p,
	int     *n
	
){         

  int     ncol, n0, nend, c, pc, cnt, cdam, csire, r, cnt2, used, o;
  ncol = n[0];
  n0 = n[1];
  nend = n[2]; 

  for(c = n0; c < nend; c++){
    pc = p[c];
    cdam = dam[c];
    csire = sire[c];
    cnt = 0;
    if(dam[c-1] == cdam && sire[c-1] == csire){
      // don't include diagonal if c-1 has a 1.0 there (stop before cnt2)
      if(c-1 < ncol) cnt2 = pc - 1; else cnt2 = pc; 
      for(r = p[c-1]; r < cnt2; r++){
        i[pc + cnt] = i[r];
        x[pc + cnt] += x[r];
        cnt++;
      }  // end for r
      cnt2 = 0;  // reset  
    } else{

      cnt = 0;
      if(cdam != -999){
        for(r = p[cdam]; r < p[cdam + 1]; r++){
          i[pc + cnt] = i[r];
          x[pc + cnt] += 0.5 * x[r];
          cnt++;
        }
      }
    
      cnt2 = 0;
      if(csire != -999){
        for(r = p[csire]; r < p[csire + 1]; r++){
          // see if this row has been used (i.e., dam had non-zero ancestry too)
          used = 0;
          for(o = 0; o < cnt; o++){
            if(i[pc + o] == i[r]){
              used++;
              x[pc + o] += 0.5 * x[r];
              break;
            }  // end if
          }  // end for o
          if(used == 0){
            i[pc + cnt + cnt2] = i[r];
            x[pc + cnt + cnt2] += 0.5 * x[r];
            cnt2++;
          }  // end if
        }  //  end for r
      }  // end if csire
    }  // end if/else same dam and sire as previous c
    
    // if diagonal of entire T is in subset, add a 1 to location of the diagonal  
    if(c < ncol){
      i[pc + cnt + cnt2] = c;
      x[pc + cnt + cnt2] += 1.0;
      p[c + 1] = pc + cnt + cnt2 + 1;
    } else{
        // otherwise "c" individual does not get a 1 on the diagonal
        //// diagonal is not part of Trow subset
        p[c + 1] = pc + cnt + cnt2;
      }  // end if/else
  } // end for c
  
}
}

