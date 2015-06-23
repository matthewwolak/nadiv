#include "dfc.h"

extern "C"{  

void dfc(
	int *dam,
	int *sire,
	int *i,
	int *j,
	int *Ni,
	int *exct
){         

  int     k, idam, isire, jdam, jsire, idgd, idgs, isgd, isgs, jdgd, jdgs, jsgd, jsgs, dfcfbfs, dfcosfs;

  for(k = 0; k < Ni[0]; k++){
     dfcfbfs = 0;
     dfcosfs = 0;
     idam = dam[i[k]];
     isire = sire[i[k]];
     jdam = dam[j[k]];
     jsire = sire[j[k]];
     if((idam == -999) || (isire == -999) || (jdam == -999) || (jsire == -999) || (idam == jdam) || (isire == jsire)){
        i[k] = 0;
     }
        else{
           idgd = dam[idam];
	   idgs = sire[idam];
	   isgd = dam[isire];
	   isgs = sire[isire];
           jdgd = dam[jdam];
	   jdgs = sire[jdam];
	   jsgd = dam[jsire];
	   jsgs = sire[jsire];
           if((idgd == -999) || (idgs == -999) || (isgd == -999) || (isgs == -999) || (jdgd == -999) || (jdgs == -999) || (jsgd == -999) || (jsgs == -999)){
              i[k] = 0;
           }
              else{
                 if((idgd == jdgd) && (idgs == jdgs) && (isgd == jsgd) && (isgs == jsgs)){
                    dfcfbfs = 1;
                    if(exct[0] == 1){
                       if(dam[idgd] != -999){
                          if(dam[idgd] == dam[jdgd]){
                             dfcfbfs = 0;
                          }
                       }
                       if(sire[idgd] != -999){
                          if(sire[idgd] == sire[jdgd]){
                             dfcfbfs = 0;
                          }
                       }

                       if(dam[idgs] != -999){
                          if(dam[idgs] == dam[jdgs]){
                             dfcfbfs = 0;
                          }
                       }
                       if(sire[idgs] != -999){
                          if(sire[idgs] == sire[jdgs]){
                             dfcfbfs = 0;
                          }
                       }

                       if(dam[isgd] != -999){
                          if(dam[isgd] == dam[jsgd]){
                             dfcfbfs = 0;
                          }
                       }
                       if(sire[isgd] != -999){
                          if(sire[isgd] == sire[jsgd]){
                             dfcfbfs = 0;
                          }
                       }

                       if(dam[isgs] != -999){
                          if(dam[isgs] == dam[jsgs]){
                             dfcfbfs = 0;
                          }
                       }
                       if(sire[isgs] != -999){
                          if(sire[isgs] == sire[jsgs]){
                             dfcfbfs = 0;
                          }
                       }
                    }
                 }
                 if((idgd == jsgd) && (idgs == jsgs) && (isgd == jdgd) && (isgs == jdgs)){
                    dfcosfs = 1;
                    if(exct[0] == 1){
                       if(dam[idgd] != -999){
                          if(dam[idgd] == dam[jsgd]){
                             dfcosfs = 0;
                          }
                       }
                       if(sire[idgd] != -999){
                          if(sire[idgd] == sire[jsgd]){
                             dfcosfs = 0;
                          }
                       }

                       if(dam[idgs] != -999){
                          if(dam[idgs] == dam[jsgs]){
                             dfcosfs = 0;
                          }
                       }
                       if(sire[idgs] != -999){
                          if(sire[idgs] == sire[jsgs]){
                             dfcosfs = 0;
                          }
                       }

                       if(dam[isgd] != -999){
                          if(dam[isgd] == dam[jdgd]){
                             dfcosfs = 0;
                          }
                       }
                       if(sire[isgd] != -999){
                          if(sire[isgd] == sire[jdgd]){
                             dfcosfs = 0;
                          }
                       }

                       if(dam[isgs] != -999){
                          if(dam[isgs] == dam[jdgs]){
                             dfcosfs = 0;   
                          }   
                       }
                       if(sire[isgs] != -999){
                          if(sire[isgs] == sire[jdgs]){
                             dfcosfs = 0;
                          }
                       }
                    }  
                 }
                 if((dfcfbfs == 1) || (dfcosfs == 1)){
                    i[k] = 1;
                 }
                    else{
                       i[k] = 0;
                    }
              }
        }

  }
        
}
}
