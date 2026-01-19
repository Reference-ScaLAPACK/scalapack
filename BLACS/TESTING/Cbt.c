#ifdef BTCINTFACE
#include "Cbt.h"
extern void Cblacs_gridinit(Int*, char*, Int, Int);
extern void Cblacs_setup(Int*, Int*);
extern void Cblacs_pinfo(Int*, Int*);
extern void Cblacs_gridmap(Int*, Int*, Int, Int, Int);
extern void Cblacs_gridexit(Int);
extern void Cblacs_abort(Int, Int);
extern void Cblacs_exit(Int);
extern void Cblacs_freebuff(Int, Int);
extern void Cblacs_gridinfo(Int, Int*, Int*, Int*, Int*);
extern void Cblacs_barrier(Int, char*);
extern Int Cblacs_pnum(Int, Int, Int);
extern void Cblacs_pcoord(Int, Int, Int*, Int*);
extern void Cblacs_get(Int, Int, Int*);
extern void Cblacs_set(Int, Int, Int*);
extern void Cigesd2d(Int, Int, Int, Int*, Int, Int, Int);
extern void Cigerv2d(Int, Int, Int, Int*, Int, Int, Int);
extern void Cigebs2d(Int, char*, char*, Int, Int, Int*, Int);
extern void Cigebr2d(Int, char*, char*, Int, Int, Int*, Int, Int, Int);
extern void Citrsd2d(Int, char*, char*, Int, Int, Int*, Int, Int, Int);
extern void Citrrv2d(Int, char*, char*, Int, Int, Int*, Int, Int, Int);
extern void Citrbs2d(Int, char*, char*, char*, char*, Int, Int, Int*, Int);
extern void Citrbr2d(Int, char*, char*, char*, char*, Int, Int, Int*, Int, Int, Int);
extern void Cigsum2d(Int, char*, char*, Int, Int, Int*, Int, Int, Int);
extern void Cigamx2d(Int, char*, char*, Int, Int, Int*, Int, Int*, Int*, Int, Int, Int);
extern void Cigamn2d(Int, char*, char*, Int, Int, Int*, Int, Int*, Int*, Int, Int, Int);
extern void Cdgesd2d(Int, Int, Int, double*, Int, Int, Int);
extern void Cdgerv2d(Int, Int, Int, double*, Int, Int, Int);
extern void Cdgebs2d(Int, char*, char*, Int, Int, double*, Int);
extern void Cdgebr2d(Int, char*, char*, Int, Int, double*, Int, Int, Int);
extern void Cdtrsd2d(Int, char*, char*, Int, Int, double*, Int, Int, Int);
extern void Cdtrrv2d(Int, char*, char*, Int, Int, double*, Int, Int, Int);
extern void Cdtrbs2d(Int, char*, char*, char*, char*, Int, Int, double*, Int);
extern void Cdtrbr2d(Int, char*, char*, char*, char*, Int, Int, double*, Int, Int, Int);
extern void Cdgsum2d(Int, char*, char*, Int, Int, double*, Int, Int, Int);
extern void Cdgamx2d(Int, char*, char*, Int, Int, double*, Int, Int*, Int*, Int, Int, Int);
extern void Cdgamn2d(Int, char*, char*, Int, Int, double*, Int, Int*, Int*, Int, Int, Int);
extern void Csgesd2d(Int, Int, Int, float*, Int, Int, Int);
extern void Csgerv2d(Int, Int, Int, float*, Int, Int, Int);
extern void Csgebs2d(Int, char*, char*, Int, Int, float*, Int);
extern void Csgebr2d(Int, char*, char*, Int, Int, float*, Int, Int, Int);
extern void Cstrsd2d(Int, char*, char*, Int, Int, float*, Int, Int, Int);
extern void Cstrrv2d(Int, char*, char*, Int, Int, float*, Int, Int, Int);
extern void Cstrbs2d(Int, char*, char*, char*, char*, Int, Int, float*, Int);
extern void Cstrbr2d(Int, char*, char*, char*, char*, Int, Int, float*, Int, Int, Int);
extern void Csgsum2d(Int, char*, char*, Int, Int, float*, Int, Int, Int);
extern void Csgamx2d(Int, char*, char*, Int, Int, float*, Int, Int*, Int*, Int, Int, Int);
extern void Csgamn2d(Int, char*, char*, Int, Int, float*, Int, Int*, Int*, Int, Int, Int);
void blacs_gridinit_(Int *ConTxt, char *order, Int *nprow, Int *npcol)
{
Cblacs_gridinit(ConTxt, order, *nprow, *npcol);
}
void blacs_setup_(Int *mypnum, Int *nprocs)
{
Cblacs_setup(mypnum, nprocs);
}
void blacs_pinfo_(Int *mypnum, Int *nprocs)
{
Cblacs_pinfo(mypnum, nprocs);
}
void blacs_gridmap_(Int *ConTxt, Int *usermap, Int *ldup, Int *nprow, Int *npcol)
{
Cblacs_gridmap(ConTxt, usermap, *ldup, *nprow, *npcol);
}
void blacs_gridexit_(Int *ConTxt)
{
Cblacs_gridexit(*ConTxt);
}
void blacs_abort_(Int *ConTxt, Int *ErrNo)
{
Cblacs_abort(*ConTxt, *ErrNo);
}
void blacs_exit_(Int *NotDone)
{
Cblacs_exit(*NotDone);
}
void blacs_freebuff_(Int *ConTxt, Int *Wait)
{
Cblacs_freebuff(*ConTxt, *Wait);
}
void blacs_gridinfo_(Int *ConTxt, Int *nprow, Int *npcol, Int *myrow, Int *mycol)
{
Cblacs_gridinfo(*ConTxt, nprow, npcol, myrow, mycol);
}
void blacs_barrier_(Int *ConTxt, char *scope)
{
Cblacs_barrier(*ConTxt, scope);
}
Int blacs_pnum_(Int *ConTxt, Int *prow, Int *pcol)
{
return( Cblacs_pnum(*ConTxt, *prow, *pcol) );
}
void blacs_pcoord_(Int *ConTxt, Int *nodenum, Int *prow, Int *pcol)
{
Cblacs_pcoord(*ConTxt, *nodenum, prow, pcol);
}
void blacs_get_(Int *ConTxt, Int *what, Int *I)
{
Cblacs_get(*ConTxt, *what, I);
}
void blacs_set_(Int *ConTxt, Int *what, Int *I)
{
Cblacs_set(*ConTxt, *what, I);
}
void igesd2d_(Int *ConTxt, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
Cigesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}
void igerv2d_(Int *ConTxt, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
Cigerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}
void igebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda)
{
Cigebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}
void igebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
Cigebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}
void itrsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
Citrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}
void itrrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
Citrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void itrbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda)
{
Citrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}
void itrbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
Citrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void igsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
Cigsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}
void igamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Cigamx2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
void igamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Cigamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
void dgesd2d_(Int *ConTxt, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
Cdgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}
void dgerv2d_(Int *ConTxt, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
Cdgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}
void dgebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda)
{
Cdgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}
void dgebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
Cdgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}
void dtrsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
Cdtrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}
void dtrrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
Cdtrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void dtrbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda)
{
Cdtrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}
void dtrbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
Cdtrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void dgsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
Cdgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}
void dgamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Cdgamx2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
void dgamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Cdgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
void sgesd2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
Csgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}
void sgerv2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
Csgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}
void sgebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda)
{
Csgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}
void sgebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
Csgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}
void strsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
Cstrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}
void strrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
Cstrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void strbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda)
{
Cstrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}
void strbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
Cstrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
void sgsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
Csgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}
void sgamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Csgamx2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
void sgamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
Csgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia, *rdest, *cdest);
}
#endif



