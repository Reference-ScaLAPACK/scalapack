#ifdef BTCINTFACE
#include "Cbt.h"

void blacs_gridinit_(Int *ConTxt, char *order, Int *nprow, Int *npcol)
{
   void Cblacs_gridinit( Int* context, char* order, Int nprow, Int npcol );

   Cblacs_gridinit(ConTxt, order, *nprow, *npcol);
}

void blacs_setup_(Int *mypnum, Int *nprocs)
{
   void Cblacs_setup( Int* mypnum, Int* nprocs );
   Cblacs_setup(mypnum, nprocs);
}

void blacs_pinfo_(Int *mypnum, Int *nprocs)
{
   void Cblacs_pinfo( Int* mypnum, Int* nprocs );
   Cblacs_pinfo(mypnum, nprocs);
}

void blacs_gridmap_(Int *ConTxt, Int *usermap, Int *ldup, Int *nprow, Int *npcol)
{
   void Cblacs_gridmap( Int* context, Int* usermap, Int ldumap, Int nprow, Int npcol );
   Cblacs_gridmap(ConTxt, usermap, *ldup, *nprow, *npcol);
}

void blacs_gridexit_(Int *ConTxt)
{
   void Cblacs_gridexit( Int context );
   Cblacs_gridexit(*ConTxt);
}

void blacs_abort_(Int *ConTxt, Int *ErrNo)
{
   void Cblacs_abort(Int ConTxt, Int ErrNo);
   Cblacs_abort(*ConTxt, *ErrNo);
}

void blacs_exit_(Int *NotDone)
{
   void Cblacs_exit(Int NotDone);
   Cblacs_exit(*NotDone);
}

void blacs_freebuff_(Int *ConTxt, Int *Wait)
{
   void Cblacs_freebuff(Int ConTxt, Int Wait);
   Cblacs_freebuff(*ConTxt, *Wait);
}

void blacs_gridinfo_(Int *ConTxt, Int *nprow, Int *npcol, Int *myrow, Int *mycol)
{
   void Cblacs_gridinfo( Int context, Int* nprow, Int* npcol, Int* myrow, Int* mycol );
   Cblacs_gridinfo(*ConTxt, nprow, npcol, myrow, mycol);
}

void blacs_barrier_(Int *ConTxt, char *scope)
{
   void Cblacs_barrier(Int ConTxt, char *scope);
   Cblacs_barrier(*ConTxt, scope);
}

Int blacs_pnum_(Int *ConTxt, Int *prow, Int *pcol)
{
   Int Cblacs_pnum( Int context, Int prow, Int pcol );
   return( Cblacs_pnum(*ConTxt, *prow, *pcol) );
}

void blacs_pcoord_(Int *ConTxt, Int *nodenum, Int *prow, Int *pcol)
{
   void Cblacs_pcoord( Int context, Int pnum, Int* prow, Int* pcol );
   Cblacs_pcoord(*ConTxt, *nodenum, prow, pcol);
}

void blacs_get_(Int *ConTxt, Int *what, Int *I)
{
   void Cblacs_get( Int context, Int what, Int* I );
   Cblacs_get(*ConTxt, *what, I);
}

void blacs_set_(Int *ConTxt, Int *what, Int *I)
{
   void Cblacs_set(Int ConTxt, Int what, Int *I);
   Cblacs_set(*ConTxt, *what, I);
}


void igesd2d_(Int *ConTxt, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cigesd2d(Int ConTxt, Int m, Int n, Int *A, Int lda, Int rdest, Int cdest);
   Cigesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void igerv2d_(Int *ConTxt, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cigerv2d(Int ConTxt, Int m, Int n, Int *A, Int lda, Int rsrc, Int csrc);
   Cigerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void igebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda)
{
   void Cigebs2d(Int ConTxt, char *scope, char *top, Int m, Int n, Int *A, Int lda);
   Cigebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void igebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cigebr2d(Int ConTxt, char *scope, char *top, Int m, Int n, Int *A, Int lda, Int rsrc, Int csrc);
   Cigebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void itrsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
   void Citrsd2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, Int *A, Int lda, Int rdest, Int cdest);
   Citrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void itrrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Citrrv2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, Int *A, Int lda, Int rsrc, Int csrc);
   Citrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void itrbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda)
{
   void Citrbs2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, Int *A, Int lda);
   Citrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void itrbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, Int *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Citrbr2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, Int *A, Int lda, Int rsrc, Int csrc);
   Citrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void igsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cigsum2d(Int ConTxt, char *scope, char *top, Int m, Int n, Int *A, Int lda, Int rdest, Int cdest);
   Cigsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void igamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Cigamx2d(Int ConTxt, char *scope, char *top, Int m, Int n, Int *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Cigamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void igamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, Int *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Cigamn2d(Int ConTxt, char *scope, char *top, Int m, Int n, Int *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Cigamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void dgesd2d_(Int *ConTxt, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cdgesd2d(Int ConTxt, Int m, Int n, double *A, Int lda, Int rdest, Int cdest);
   Cdgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void dgerv2d_(Int *ConTxt, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cdgerv2d(Int ConTxt, Int m, Int n, double *A, Int lda, Int rsrc, Int csrc);
   Cdgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void dgebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda)
{
   void Cdgebs2d(Int ConTxt, char *scope, char *top, Int m, Int n, double *A, Int lda);
   Cdgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void dgebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cdgebr2d(Int ConTxt, char *scope, char *top, Int m, Int n, double *A, Int lda, Int rsrc, Int csrc);
   Cdgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void dtrsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cdtrsd2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, double *A, Int lda, Int rdest, Int cdest);
   Cdtrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void dtrrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cdtrrv2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, double *A, Int lda, Int rsrc, Int csrc);
   Cdtrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void dtrbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda)
{
   void Cdtrbs2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, double *A, Int lda);
   Cdtrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void dtrbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, double *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cdtrbr2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, double *A, Int lda, Int rsrc, Int csrc);
   Cdtrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void dgsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cdgsum2d(Int ConTxt, char *scope, char *top, Int m, Int n, double *A, Int lda, Int rdest, Int cdest);
   Cdgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void dgamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Cdgamx2d(Int ConTxt, char *scope, char *top, Int m, Int n, double *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Cdgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void dgamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, double *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Cdgamn2d(Int ConTxt, char *scope, char *top, Int m, Int n, double *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Cdgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void sgesd2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
   void Csgesd2d(Int ConTxt, Int m, Int n, float *A, Int lda, Int rdest, Int cdest);
   Csgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void sgerv2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Csgerv2d(Int ConTxt, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Csgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void sgebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda)
{
   void Csgebs2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda);
   Csgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void sgebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Csgebr2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Csgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void strsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cstrsd2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rdest, Int cdest);
   Cstrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void strrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cstrrv2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Cstrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void strbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda)
{
   void Cstrbs2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, float *A, Int lda);
   Cstrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void strbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cstrbr2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Cstrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void sgsum2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
   void Csgsum2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda, Int rdest, Int cdest);
   Csgsum2d(*ConTxt, scope, top, *m, *n, A, *lda, *rdest, *cdest);
}

void sgamx2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Csgamx2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Csgamx2d(*ConTxt, scope, top, *m, *n, A, *lda,  rA, cA, *ldia,
            *rdest, *cdest);
}

void sgamn2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rA, Int *cA, Int *ldia, Int *rdest, Int *cdest)
{
   void Csgamn2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest);
   Csgamn2d(*ConTxt, scope, top, *m, *n, A, *lda, rA, cA, *ldia,
            *rdest, *cdest);
}

void cgesd2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
   void Ccgesd2d(Int ConTxt, Int m, Int n, float *A, Int lda, Int rdest, Int cdest);
   Ccgesd2d(*ConTxt, *m, *n, A, *lda, *rdest, *cdest);
}

void cgerv2d_(Int *ConTxt, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Ccgerv2d(Int ConTxt, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Ccgerv2d(*ConTxt, *m, *n, A, *lda, *rsrc, *csrc);
}

void cgebs2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda)
{
   void Ccgebs2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda);
   Ccgebs2d(*ConTxt, scope, top, *m, *n, A, *lda);
}

void cgebr2d_(Int *ConTxt, char *scope, char *top, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Ccgebr2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Ccgebr2d(*ConTxt, scope, top, *m, *n, A, *lda, *rsrc, *csrc);
}

void ctrsd2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rdest, Int *cdest)
{
   void Cctrsd2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rdest, Int cdest);
   Cctrsd2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rdest, *cdest);
}

void ctrrv2d_(Int *ConTxt, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cctrrv2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Cctrrv2d(*ConTxt, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}

void ctrbs2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda)
{
   void Cctrbs2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, float *A, Int lda);
   Cctrbs2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda);
}

void ctrbr2d_(Int *ConTxt, char *scope, char *top, char *uplo, char *diag, Int *m, Int *n, float *A, Int *lda, Int *rsrc, Int *csrc)
{
   void Cctrbr2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag, Int m, Int n, float *A, Int lda, Int rsrc, Int csrc);
   Cctrbr2d(*ConTxt, scope, top, uplo, diag, *m, *n, A, *lda, *rsrc, *csrc);
}
#endif
