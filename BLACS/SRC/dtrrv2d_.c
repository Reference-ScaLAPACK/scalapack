#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cdtrrv2d(Int ConTxt, char *uplo, char *diag, Int m, Int n, double *A,
              Int lda, Int rsrc, Int csrc)
#else
F_VOID_FUNC dtrrv2d_(Int *ConTxt, F_CHAR uplo, F_CHAR diag, Int *m, Int *n,
                     double *A, Int *lda, Int *rsrc, Int *csrc)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Globally-blocking point to point trapezoidal double precision receive.
 *
 *  Arguments
 *  =========
 *
 *  ConTxt  (input) Ptr to Int
 *          Index into MyConTxts00 (my contexts array).
 *
 *  UPLO    (input) Ptr to char
 *          Specifies the part of the matrix to be sent.
 *          = 'U':      Upper trapezoidal part
 *          ELSE :      Lower trapezoidal part
 *
 *  DIAG    (input) Ptr to char
 *          Specifies whether the matrix is unit diagonal or not.
 *          = 'U':      Matrix is unit diagonal, diagonal not communicated.
 *          ELSE :      Matrix is not unit diagonal, diagonal is communicated.
 *
 *  M       (input) Ptr to Int
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) Ptr to Int
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (output) Ptr to double precision two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *          If UPLO = 'U', only the upper trapezoid is accessed;
 *          if UPLO = 'L', only the lower trapezoid is accessed.
 *
 *  LDA     (input) Ptr to Int
 *          The leading dimension of the array A.  LDA >= M.
 *
 *
 *  RSRC    (input) Ptr to Int
 *          The process row of the source of the matrix.
 *
 *  CSRC    (input) Ptr to Int
 *          The process column of the source of the matrix.
 *
 *
 * ------------------------------------------------------------------------
 */
{
/*
 *  Prototypes and variable declarations
 */
   void BI_ArgCheck(Int, Int, char *, char, char, char, Int, Int, Int, Int,
                    Int *, Int *);
   MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *, char, char, Int, Int, Int,
                                   MPI_Datatype, Int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Srecv(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(Int);
   Int BI_BuffIsFree(BLACBUFF *, Int);
   Int tuplo, tdiag, tlda;
   Int ierr, length;
   BLACBUFF *bp;
   MPI_Datatype MatTyp;
   BLACSCONTEXT *ctxt;
   extern BLACBUFF BI_AuxBuff, *BI_ActiveQ;

   MGetConTxt(Mpval(ConTxt), ctxt);
   tdiag = F2C_CharTrans(diag);
   tuplo = F2C_CharTrans(uplo);
   tdiag = Mlowcase(tdiag);
   tuplo = Mlowcase(tuplo);

#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_RV, __FILE__, 'a', tuplo, tdiag, Mpval(m),
               Mpval(n), Mpval(lda), 1, Mpaddress(rsrc), Mpaddress(csrc));
#endif
   if (Mpval(lda) < Mpval(m)) tlda = Mpval(m);
   else tlda = Mpval(lda);
   ctxt->scp = &ctxt->pscp;

   MatTyp = BI_GetMpiTrType(ctxt, tuplo, tdiag, Mpval(m), Mpval(n), tlda,
                            MPI_DOUBLE, &BI_AuxBuff.N);
   BI_AuxBuff.Buff = (char *) A;
   BI_AuxBuff.dtype = MatTyp;
   BI_Srecv(ctxt, Mkpnum(ctxt, Mpval(rsrc), Mpval(csrc)), PT2PTID, &BI_AuxBuff);
   ierr=BI_MPI_TYPE_FREE(&MatTyp);
   if (BI_ActiveQ) BI_UpdateBuffs(NULL);
}
