#include "Bdef.h"

#if (INTFACE == C_CALL)
void Csgerv2d(int ConTxt, int m, int n, float *A, int lda, int rsrc, int csrc)
#else
F_VOID_FUNC sgerv2d_(int *ConTxt, int *m, int *n, float *A, int *lda,
                     int *rsrc, int *csrc)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Globally-blocking point to point general real receive.
 *
 *  Arguments
 *  =========
 *
 *  ConTxt  (input) Ptr to int
 *          Index into MyConTxts00 (my contexts array).
 *
 *  M       (input) Ptr to int
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) Ptr to int
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (output) Ptr to real two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *
 *  LDA     (input) Ptr to int
 *          The leading dimension of the array A.  LDA >= M.
 *
 *
 *  RSRC    (input) Ptr to int
 *          The process row of the source of the matrix.
 *
 *  CSRC    (input) Ptr to int
 *          The process column of the source of the matrix.
 *
 *
 * ------------------------------------------------------------------------
 */
{
/*
 *  Prototypes and variable declarations
 */
   void BI_ArgCheck(int, int, char *, char, char, char, int, int, int, int,
                    int *, int *);
   MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *, int, int, int,
                                   MPI_Datatype, int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Srecv(BLACSCONTEXT *, int, int, BLACBUFF *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);
   int tlda;
   int ierr;
   MPI_Datatype MatTyp;
   BLACSCONTEXT *ctxt;
   extern BLACBUFF BI_AuxBuff, *BI_ActiveQ;

   MGetConTxt(Mpval(ConTxt), ctxt);
#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_RV, __FILE__, 'a', 'u', 'u', Mpval(m),
               Mpval(n), Mpval(lda), 1, Mpaddress(rsrc), Mpaddress(csrc));
#endif
   if (Mpval(lda) < Mpval(m)) tlda = Mpval(m);
   else tlda = Mpval(lda);
   ctxt->scp = &ctxt->pscp;

   MatTyp = BI_GetMpiGeType(ctxt, Mpval(m), Mpval(n), tlda,
                            MPI_FLOAT, &BI_AuxBuff.N);
   BI_AuxBuff.Buff = (char *) A;
   BI_AuxBuff.dtype = MatTyp;
   BI_Srecv(ctxt, Mkpnum(ctxt, Mpval(rsrc), Mpval(csrc)), PT2PTID, &BI_AuxBuff);
   ierr=BI_MPI_TYPE_FREE(&MatTyp);
   if (BI_ActiveQ) BI_UpdateBuffs(NULL);
}
