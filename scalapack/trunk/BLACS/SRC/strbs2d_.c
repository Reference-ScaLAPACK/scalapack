#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cstrbs2d(int ConTxt, char *scope, char *top, char *uplo, char *diag,
              int m, int n, float *A, int lda)
#else
F_VOID_FUNC strbs2d_(int *ConTxt, F_CHAR scope, F_CHAR top, F_CHAR uplo,
                     F_CHAR diag, int *m, int *n, float *A, int *lda)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Broadcast/send for trapezoidal real arrays.
 *
 *  Arguments
 *  =========
 *
 *  ConTxt  (input) Ptr to int
 *          Index into MyConTxts00 (my contexts array).
 *
 *  SCOPE   (input) Ptr to char
 *          Limit the scope of the operation.
 *          = 'R' :   Operation is performed by a process row.
 *          = 'C' :   Operation is performed by a process column.
 *          = 'A' :   Operation is performed by all processes in grid.
 *
 *  TOP     (input) Ptr to char
 *          Controls fashion in which messages flow within the operation.
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
 *  A       (input) Ptr to real two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *          If UPLO = 'U', only the upper trapezoid is accessed;
 *          if UPLO = 'L', only the lower trapezoid is accessed.
 *
 *  LDA     (input) Ptr to int
 *          The leading dimension of the array A.  LDA >= M.
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_ArgCheck(int, int, char *, char, char, char, int, int, int, int,
                    int *, int *);
   int BI_HypBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR);
   void BI_IdringBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, int);
   void BI_SringBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR);
   void BI_MpathBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, int);
   void BI_TreeBS(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, int);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);
   MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *, char, char, int, int, int,
                                   MPI_Datatype, int *);
   BLACBUFF *BI_Pack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Ssend(BLACSCONTEXT *, int, int, BLACBUFF *);
   void BI_Asend(BLACSCONTEXT *, int, int, BLACBUFF *);

   char ttop, tscope, tuplo, tdiag;
   int error, tlda;
   MPI_Datatype MatTyp;
   SDRVPTR send;
   BLACBUFF *bp;
   BLACSCONTEXT *ctxt;
   extern BLACBUFF BI_AuxBuff, *BI_ActiveQ;
/*
 * get context, lowcase char variables, and perform parameter checking
 */
   MGetConTxt(Mpval(ConTxt), ctxt);
   ttop = F2C_CharTrans(top);
   ttop = Mlowcase(ttop);
   tscope = F2C_CharTrans(scope);
   tscope = Mlowcase(tscope);
   tuplo = F2C_CharTrans(uplo);
   tuplo = Mlowcase(tuplo);
   tdiag = F2C_CharTrans(diag);
   tdiag = Mlowcase(tdiag);
#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_BS, __FILE__, 'a', tuplo, tdiag, Mpval(m),
               Mpval(n), Mpval(lda), 0, NULL, NULL);
#endif
/*
 *  If the user has set the default broadcast topology, use it instead of
 *  BLACS default
 */
#ifdef DefBSTop
   if (ttop == ' ') ttop = DefBSTop;
#endif
   if (Mpval(lda) < Mpval(m)) tlda = Mpval(m);
   else tlda = Mpval(lda);

   switch(tscope)
   {
   case 'r':
      ctxt->scp = &ctxt->rscp;
      break;
   case 'c':
      ctxt->scp = &ctxt->cscp;
      break;
   case 'a':
      ctxt->scp = &ctxt->ascp;
      break;
   default:
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown scope '%c'",
                  tscope);
   }
   MatTyp = BI_GetMpiTrType(ctxt, tuplo, tdiag, Mpval(m), Mpval(n), tlda,
                            MPI_FLOAT, &BI_AuxBuff.N);
/*
 * If using default topology, use MPI native broadcast
 */
   if (ttop == ' ')
   {
      error=MPI_Bcast(A, BI_AuxBuff.N, MatTyp, ctxt->scp->Iam, ctxt->scp->comm);
      error=BI_MPI_TYPE_FREE(&MatTyp);
      if (BI_ActiveQ) BI_UpdateBuffs(NULL);
      return;
   }
/*
 * If MPI handles non-contiguous buffering well, always use MPI data types
 * instead of packing
 */
#ifdef MpiBuffGood
   send = BI_Ssend;
   BI_AuxBuff.Buff = (char *) A;
   BI_AuxBuff.dtype = MatTyp;
   bp = &BI_AuxBuff;
#endif
/*
 * Pack and use non-blocking sends for broadcast if MPI's data types aren't
 * more efficient
 */
#ifndef MpiBuffGood
   send = BI_Asend;
   bp = BI_Pack(ctxt, (BVOID *) A, NULL, MatTyp);
#endif

/*
 * Call correct topology for BS/BR
 */
   switch(ttop)
   {
   case 'h':
      error = BI_HypBS(ctxt, bp, send);
      if (error == NPOW2) BI_TreeBS(ctxt, bp, send, 2);
      break;
   case '1':
   case '2':
   case '3':
   case '4':
   case '5':
   case '6':
   case '7':
   case '8':
   case '9':
      BI_TreeBS(ctxt, bp, send, ttop-47);
      break;
   case 't':
      BI_TreeBS(ctxt, bp, send, ctxt->Nb_bs);
      break;
   case 'i':
      BI_IdringBS(ctxt, bp, send, 1);
      break;
   case 'd':
      BI_IdringBS(ctxt, bp, send, -1);
      break;
   case 's':
      BI_SringBS(ctxt, bp, send);
      break;
   case 'f':
      BI_MpathBS(ctxt, bp, send, FULLCON);
      break;
   case 'm':
      BI_MpathBS(ctxt, bp, send, ctxt->Nr_bs);
      break;
   default :
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown topology '%c'",
                  ttop);
   }

   error=BI_MPI_TYPE_FREE(&MatTyp);
   if (bp == &BI_AuxBuff)
   {
      if (BI_ActiveQ) BI_UpdateBuffs(NULL);
   }
   else BI_UpdateBuffs(bp);
}  /* end  strbs2d_  */
