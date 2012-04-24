#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cztrsd2d(int ConTxt, char *uplo, char *diag, int m, int n, double *A,
              int lda, int rdest, int cdest)
#else
F_VOID_FUNC ztrsd2d_(int *ConTxt, F_CHAR uplo, F_CHAR diag, int *m, int *n,
                     double *A, int *lda, int *rdest, int *cdest)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Locally-blocking point-to-point trapezoidal double complex send.
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
 *  A       (input) Ptr to double complex two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *          If UPLO = 'U', only the upper trapezoid is accessed;
 *          if UPLO = 'L', only the lower trapezoid is accessed.
 *
 *  LDA     (input) Ptr to int
 *          The leading dimension of the array A.  LDA >= M.
 *
 *  RDEST   (input) Ptr to int
 *          The process row of the destination process.
 *
 *  CDEST   (input) Ptr to int
 *          The process column of the destination process.
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_ArgCheck(int, int, char *, char, char, char, int, int, int, int,
                    int *, int *);
   MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *, char, char, int, int, int,
                                   MPI_Datatype, int *);
   BLACBUFF *BI_Pack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Ssend(BLACSCONTEXT *, int, int, BLACBUFF *);
   void BI_Asend(BLACSCONTEXT *, int, int, BLACBUFF *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);

   char tuplo, tdiag;
   int dest, length, tlda, ierr;
   BLACBUFF *bp;
   BLACSCONTEXT *ctxt;
   MPI_Datatype MatTyp;
   extern BLACBUFF BI_AuxBuff, *BI_ActiveQ;

   MGetConTxt(Mpval(ConTxt), ctxt);
   tuplo = F2C_CharTrans(uplo);
   tdiag = F2C_CharTrans(diag);
   tuplo = Mlowcase(tuplo);
   tdiag = Mlowcase(tdiag);

#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_SD, "ZTRSD2D", 'a', tuplo, tdiag, Mpval(m),
               Mpval(n), Mpval(lda), 1, Mpaddress(rdest), Mpaddress(cdest));
#endif
   if (Mpval(lda) < Mpval(m)) tlda = Mpval(m);
   else tlda = Mpval(lda);
   dest = Mvkpnum(ctxt, Mpval(rdest), Mpval(cdest));
   ctxt->scp = &ctxt->pscp;

   MatTyp = BI_GetMpiTrType(ctxt, tuplo, tdiag, Mpval(m), Mpval(n), tlda,
                            MPI_DOUBLE_COMPLEX, &BI_AuxBuff.N);
#ifdef SndIsLocBlk
   BI_AuxBuff.Buff = (char *) A;
   BI_AuxBuff.dtype = MatTyp;
   BI_Ssend(ctxt, dest, PT2PTID, &BI_AuxBuff);
#else
   bp = BI_Pack(ctxt, (BVOID *) A, NULL, MatTyp);
   BI_Asend(ctxt, Mkpnum(ctxt, Mpval(rdest), Mpval(cdest)), PT2PTID, bp);
#endif
   ierr=BI_MPI_TYPE_FREE(&MatTyp);

/*
 * Having started the async send, update the buffers (reform links, check if
 * active buffers have become inactive, etc.)
 */
#ifdef SndIsLocBlk
   if (BI_ActiveQ) BI_UpdateBuffs(NULL);
#else
   BI_UpdateBuffs(bp);
#endif
}  /* end of ztrsd2d */
