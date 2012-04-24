#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cztrrv2d(int ConTxt, char *uplo, char *diag, int m, int n, double *A,
              int lda, int rsrc, int csrc)
#else
F_VOID_FUNC ztrrv2d_(int *ConTxt, F_CHAR uplo, F_CHAR diag, int *m, int *n,
                     double *A, int *lda, int *rsrc, int *csrc)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Globally-blocking point to point trapezoidal double complex receive.
 *
 *  Arguments
 *  =========
 *
 *  ConTxt  (input) Ptr to int
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
 *  M       (input) Ptr to int
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) Ptr to int
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (output) Ptr to double complex two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *          If UPLO = 'U', only the upper trapezoid is accessed;
 *          if UPLO = 'L', only the lower trapezoid is accessed.
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
   MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *, char, char, int, int, int,
                                   MPI_Datatype, int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Srecv(BLACSCONTEXT *, int, int, BLACBUFF *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);
   int tuplo, tdiag, tlda;
   int ierr, length;
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
                            MPI_DOUBLE_COMPLEX, &BI_AuxBuff.N);
   BI_AuxBuff.Buff = (char *) A;
   BI_AuxBuff.dtype = MatTyp;
   BI_Srecv(ctxt, Mkpnum(ctxt, Mpval(rsrc), Mpval(csrc)), PT2PTID, &BI_AuxBuff);
   ierr=BI_MPI_TYPE_FREE(&MatTyp);
   if (BI_ActiveQ) BI_UpdateBuffs(NULL);
}
