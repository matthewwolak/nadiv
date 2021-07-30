#define _NADIV_H
#include "cs.h"
/* #include "R.h" included by cs.h */ 
#include "Rmath.h" 

//   M&L 1992 algorithm (for `ainvml`)
//   as presented in Mrode 2005
//// replaces elements of f and dii with calculated values in place
void ml(int *dam, int *sire,
	double *f, double *dii,
	int n, int g, int fmiss);

