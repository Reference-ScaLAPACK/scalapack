#ifdef BTCINTFACE
#include "Cbt.h"

void blacs_gridinit_(ConTxt, order, nprow, npcol)
Int *ConTxt;
char *order;
Int *nprow;
Int  *npcol;
{
   void Cblacs_gridinit();

   Cblacs_gridinit(ConTxt, order, *nprow, *npcol);
}

void blacs_setup_(mypnum, nprocs)
Int *mypnum;
Int  *nprocs;
{
   void Cblacs_setup();
   Cblacs_setup(mypnum, nprocs);
}

void blacs_pinfo_(mypnum, nprocs)
Int *mypnum;
Int  *nprocs;
{
   void Cblacs_pinfo();
   Cblacs_pinfo(mypnum, nprocs);
}

void blacs_gridmap_(ConTxt, usermap, ldup, nprow, npcol)
Int *ConTxt;
Int *usermap;
Int *ldup;
Int *nprow;
Int  *npcol;
{
   void Cblacs_gridmap();
   Cblacs_gridmap(ConTxt, usermap, *ldup, *nprow, *npcol);
}

void blacs_gridexit_(ConTxt)
Int  *ConTxt;
{
   void Cblacs_gridexit();
   Cblacs_gridexit(*ConTxt);
}

void blacs_abort_(ConTxt, ErrNo)
Int *ConTxt;
Int  *ErrNo;
{
   void Cblacs_abort();
   Cblacs_abort(*ConTxt, *ErrNo);
}

void blacs_exit_(NotDone)
Int  *NotDone;
{
   void Cblacs_exit();
   Cblacs_exit(*NotDone);
}

void blacs_freebuff_(ConTxt, Wait)
Int *ConTxt;
Int  *Wait;
{
   void Cblacs_freebuff();
   Cblacs_freebuff(*ConTxt, *Wait);
}

void blacs_gridinfo_(ConTxt, nprow, npcol, myrow, mycol)
Int *ConTxt;
Int *nprow;
Int *npcol;
Int *myrow;
Int  *mycol;
{
   void Cblacs_gridinfo();
   Cblacs_gridinfo(*ConTxt, nprow, npcol, myrow, mycol);
}

void blacs_barrier_(ConTxt, scope)
Int *ConTxt;
char  *scope;
{
   void Cblacs_barrier();
   Cblacs_barrier(*ConTxt, scope);
}

Int blacs_pnum_(ConTxt, prow, pcol)
Int *ConTxt;
Int *prow;
Int  *pcol;
{
   Int Cblacs_pnum();
   return( Cblacs_pnum(*ConTxt, *prow, *pcol) );
}

void blacs_pcoord_(ConTxt, nodenum, prow, pcol)
Int *ConTxt;
Int *nodenum;
Int *prow;
Int  *pcol;
{
   void Cblacs_pcoord();
   Cblacs_pcoord(*ConTxt, *nodenum, prow, pcol);
}

void blacs_get_(ConTxt, what, I)
Int *ConTxt;
Int *what;
Int  *I;
{
   void Cblacs_get();
   Cblacs_get(*ConTxt, *what, I);
}

void blacs_set_(ConTxt, what, I)
Int *ConTxt;
Int *what;
Int  *I;
{
   void Cblacs_set();
   Cblacs_set(*ConTxt, *what, I);
}


void igesd2d_(ConTxt, m, n, A, lda, rdest, cdest)
Int *ConTxt;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cigesd2d();
   Cigesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void igerv2d_(ConTxt, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cigerv2d();
   Cigerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void igebs2d_(ConTxt, scope, top, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
Int *A;
Int  *lda;
{
   void Cigebs2d();
   Cigebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void igebr2d_(ConTxt, scope, top, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cigebr2d();
   Cigebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void itrsd2d_(ConTxt, uplo, diag, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Citrsd2d();
   Citrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void itrrv2d_(ConTxt, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Citrrv2d();
   Citrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void itrbs2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
Int *A;
Int  *lda;
{
   void Citrbs2d();
   Citrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void itrbr2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Citrbr2d();
   Citrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void igsum2d_(ConTxt, scope, top, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cigsum2d();
   Cigsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void igamx2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Cigamx2d();
   Cigamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void igamn2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
Int *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Cigamn2d();
   Cigamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void dgesd2d_(ConTxt, m, n, A, lda, rdest, cdest)
Int *ConTxt;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cdgesd2d();
   Cdgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void dgerv2d_(ConTxt, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cdgerv2d();
   Cdgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void dgebs2d_(ConTxt, scope, top, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int  *lda;
{
   void Cdgebs2d();
   Cdgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void dgebr2d_(ConTxt, scope, top, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cdgebr2d();
   Cdgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void dtrsd2d_(ConTxt, uplo, diag, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cdtrsd2d();
   Cdtrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void dtrrv2d_(ConTxt, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cdtrrv2d();
   Cdtrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void dtrbs2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int  *lda;
{
   void Cdtrbs2d();
   Cdtrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void dtrbr2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cdtrbr2d();
   Cdtrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void dgsum2d_(ConTxt, scope, top, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cdgsum2d();
   Cdgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void dgamx2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Cdgamx2d();
   Cdgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void dgamn2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Cdgamn2d();
   Cdgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void sgesd2d_(ConTxt, m, n, A, lda, rdest, cdest)
Int *ConTxt;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Csgesd2d();
   Csgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void sgerv2d_(ConTxt, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Csgerv2d();
   Csgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void sgebs2d_(ConTxt, scope, top, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int  *lda;
{
   void Csgebs2d();
   Csgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void sgebr2d_(ConTxt, scope, top, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Csgebr2d();
   Csgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void strsd2d_(ConTxt, uplo, diag, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cstrsd2d();
   Cstrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void strrv2d_(ConTxt, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cstrrv2d();
   Cstrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void strbs2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int  *lda;
{
   void Cstrbs2d();
   Cstrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void strbr2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cstrbr2d();
   Cstrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void sgsum2d_(ConTxt, scope, top, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Csgsum2d();
   Csgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void sgamx2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Csgamx2d();
   Csgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void sgamn2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Csgamn2d();
   Csgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void cgesd2d_(ConTxt, m, n, A, lda, rdest, cdest)
Int *ConTxt;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Ccgesd2d();
   Ccgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void cgerv2d_(ConTxt, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Ccgerv2d();
   Ccgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void cgebs2d_(ConTxt, scope, top, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int  *lda;
{
   void Ccgebs2d();
   Ccgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void cgebr2d_(ConTxt, scope, top, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Ccgebr2d();
   Ccgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void ctrsd2d_(ConTxt, uplo, diag, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cctrsd2d();
   Cctrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void ctrrv2d_(ConTxt, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cctrrv2d();
   Cctrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void ctrbs2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int  *lda;
{
   void Cctrbs2d();
   Cctrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void ctrbr2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cctrbr2d();
   Cctrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void cgsum2d_(ConTxt, scope, top, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Ccgsum2d();
   Ccgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void cgamx2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Ccgamx2d();
   Ccgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void cgamn2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
float *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Ccgamn2d();
   Ccgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void zgesd2d_(ConTxt, m, n, A, lda, rdest, cdest)
Int *ConTxt;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Czgesd2d();
   Czgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void zgerv2d_(ConTxt, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Czgerv2d();
   Czgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void zgebs2d_(ConTxt, scope, top, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int  *lda;
{
   void Czgebs2d();
   Czgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void zgebr2d_(ConTxt, scope, top, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Czgebr2d();
   Czgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void ztrsd2d_(ConTxt, uplo, diag, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Cztrsd2d();
   Cztrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void ztrrv2d_(ConTxt, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cztrrv2d();
   Cztrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void ztrbs2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int  *lda;
{
   void Cztrbs2d();
   Cztrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void ztrbr2d_(ConTxt, scope, top, uplo, diag, m, n, A, lda, rsrc, csrc)
Int *ConTxt;
char *scope;
char *top;
char *uplo;
char *diag;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rsrc;
Int  *csrc;
{
   void Cztrbr2d();
   Cztrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void zgsum2d_(ConTxt, scope, top, m, n, A, lda, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rdest;
Int  *cdest;
{
   void Czgsum2d();
   Czgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void zgamx2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Czgamx2d();
   Czgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void zgamn2d_(ConTxt, scope, top, m, n, A, lda, rA, cA, ldia, rdest, cdest)
Int *ConTxt;
char *scope;
char *top;
Int *m;
Int *n;
double *A;
Int *lda;
Int *rA;
Int *cA;
Int *ldia;
Int *rdest;
Int  *cdest;
{
   void Czgamn2d();
   Czgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}
#endif
