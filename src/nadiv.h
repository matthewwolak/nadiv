#define _NADIV_H
#include "cs.h"
/* #include "R.h" included by cs.h */ 
#include "Rmath.h" 

#ifdef __cplusplus
extern "C" {
#endif


cs *cs_cbind(const cs *A, const cs *B);
/* Returns the two matrices A and B column bound*/
void cs_cov2cor(const cs *A);
/* transforms a dense covariance matrix A into a correlation matrix  */
double cs_dcmvnorm(const cs *beta,  const cs *mu, const cs *M, int *keep, int nkeep, int *cond, int ncond);
/* log-density of beta[keep] in the normal with mean mu and covariance M conditional on beta[cond] */
cs *cs_directsum(cs **KGinv, int nG, int nGR);
/* returns the direct sum of KGinv[1] KGinv[2] ... KGinv[nG]*/
void cs_directsumupdate(cs **KGinv, int nG, int nGR, const cs *C);
/* overwrites C with the direct sum of KGinv[1] KGinv[2] ... KGinv[nG]*/
double cs_dmvnorm(const cs *beta,  const cs *mu, double ldet, const cs *Minv);
/* log-density of beta in the normal with mean mu and inverse covariance Minv*/
cs *cs_inv(const cs *C);
/* returns the inverse of the dense matrix C*/
double cs_invR(const cs *C, const cs *A);
/* overwrites A with the inverse of the dense matrix C*/
cs *cs_kroneckerA(const cs *G, const cs *A);
/* forms the kronecker product of G and A*/
void cs_kroneckerAupdate(const cs *G, const cs *A, const cs *C);
/* overwrites C with the kronecker product of G and A*/
cs *cs_kroneckerI(const cs *A, int nI);
/* forms the kronecker product of the dense matrix A and an identity matrix with dimension nI*/
void cs_kroneckerIupdate(const cs *A, int nI, const cs*C);
/* overwrites C with the kronecker product of the dense matrix A and an identity matrix with dimension nI*/
cs *cs_kroneckerSI(const cs *A, int nI);
/* forms the kronecker product of the sparse matrix A and an identity matrix with dimension nI*/
void cs_kroneckerSIupdate(const cs *A, int nI, const cs*C);
/* overwrites C with the kronecker product of the sparse matrix A and an identity matrix with dimension nI*/
cs *cs_kroneckerD(const cs *A, int nI, double *diag, int reciprocal);
/* forms the kronecker product of the dense matrix A and a diagonal matrix with dimension nI and diag along the diagonal*/
void cs_kroneckerDupdate(const cs *A, int nI, double *diag, const cs *C, int reciprocal);
/* overwrites C with the kronecker product of the dense matrix A and a diagonal matrix with dimension nI and diag along the diagonal*/
cs *cs_kroneckerDI(double *D, int n, int nI);
/* forms the kronecker product of a diagonal matrix (with diagonal elements D) and a diagonal matrix of dimension nI*/
cs *cs_kroneckerDA(double *D, int n, const cs *A);
/* forms the kronecker product of a diagonal matrix (with diagonal elements D) and a sparse matrix A*/
cs *cs_omega(cs **KGinv, int nG, cs *pvB);
/* returns the direct sum of pvB KGinv[1] KGinv[2] ... KGinv[nG] */
void cs_omegaupdate(cs **KGinv, int nG, cs *pvB, const cs *C);
/* overwrites C with the direct sum of pvB KGinv[1] KGinv[2] ... KGinv[nG] */
cs *cs_rCinvwishart(const cs *A, double nu, int split, const cs *CM);
/* samples from the conditional inverse Wishart given *inverse* scale matrix A*/
cs *cs_rSinvwishart(const cs *A, double nu, int split);
/* samples upper sub-matrix (of a us-identity direct sum matrix) from the inverse Wishart given *inverse* scale matrix A*/
cs *cs_rinvwishart(const cs *A, double nu, const css *As);
/* samples from the inverse Wishart given *inverse* scale matrix A*/
cs *cs_rRsubinvwishart(const cs *A, double nu, int split, double nuR, const cs *pG, const cs *CM);
/* samples a correlation sub-matrix, then conditions on it and samples from inverse-Wishart given scale matrix A (CM is the correlation sub-matrix from the previous iteration)*/
cs *cs_rR(const cs *A, double nu, double nuR, const css *As, const cs *Roldinv, double Roldldet, const cs *pG);
/* samples a correlation matrix given scale matrix A */
cs *cs_rwishart(const cs *A, double nu, const css *As);
/* samples from the Wishart*/
void cs_sortdv(const cs *A);
/* sorts a dense vector*/
double dcutpoints(const cs *liab, double *yP, int *observed, int start,int finish, double *oldcutopints, double *newcutopints, int stcutpoints, int ncutpoints, double sdcp, double sdl);
/* log MH ratio of new cutpoints to old cutpoints given liabilities liab and all yP between start and finish for which observed==1. sdcp and sdl are the proposal standard deviations for cutpoints and thresholds. */
double pcmvnorm(const cs *predi, const cs *linki, const cs *G, int keep, double lower, double upper);
/* log cumlative distribution function between lower and upper for variable[keep] conditioning on linki[-keep]
linki is multivariate normal with mean predi and covariance G */
double rtcmvnorm(const cs *predi, const cs *linki, const cs *G, int keep, double lower, double upper);
/* sample variable[keep] truncated between lower and upper conditioning on linki[-keep]
linki is multivariate normal with mean predi and covariance G */
double rtnorm(double mu, double sd, double lower, double upper);
/* sample form the truncated (between lower and upper) normal with mean mu and standard devaition sd */
cs *cs_initialize(double *x, int *p, int *i, int n, int m, int nzmax);
/* allocate and fill a cs sparse matrix */
cs *cs_rAnte(const cs *location, int start, int dimG, int nlGR, int nk, const cs *pmuAnte, const cs *pvAnte, const cs *Ainv, int Aterm, double *ivar, int cvar, const cs *pG, double pnG);
/* sample an antedependnce structure */
cs *cs_schur(const cs *A,  int split, const cs *beta);
/* forms the Schur complement for the dense mxm matrix: A_22-A_21%*%solve(A_11)%*%A_12 where submatrices are defined by split. Also overwrites beta_rr with A_21%*%solve(A_11) */

#ifdef __cplusplus
}
#endif

