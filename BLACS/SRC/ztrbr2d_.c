#include "Bdef.h"

#if (INTFACE == C_CALL)
void Cztrbr2d(Int ConTxt, char *scope, char *top, char *uplo, char *diag,
              Int m, Int n, double *A, Int lda, Int rsrc, Int csrc)
#else
F_VOID_FUNC ztrbr2d_(Int *ConTxt, F_CHAR scope, F_CHAR top, F_CHAR uplo,
                     F_CHAR diag, Int *m, Int *n, double *A, Int *lda,
                     Int *rsrc, Int *csrc)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Broadcast/receive for trapezoidal double complex arrays.
 *
 *  Arguments
 *  =========
 *
 *  ConTxt  (input) Ptr to Int
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
 *  M       (input) Ptr to Int
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) Ptr to Int
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (output) Ptr to double complex two dimensional array
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
   void BI_ArgCheck(Int, Int, char *, char, char, char, Int, Int, Int, Int,
                    Int *, Int *);
   Int BI_HypBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int);
   void BI_IdringBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int, Int);
   void BI_SringBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int);
   void BI_MpathBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int, Int);
   void BI_TreeBR(BLACSCONTEXT *, BLACBUFF *, SDRVPTR, Int, Int);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(Int);
   Int BI_BuffIsFree(BLACBUFF *, Int);
   MPI_Datatype BI_GetMpiTrType(BLACSCONTEXT *, char, char, Int, Int, Int,
                                   MPI_Datatype, Int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_Ssend(BLACSCONTEXT *, Int, Int, BLACBUFF *);
   void BI_Asend(BLACSCONTEXT *, Int, Int, BLACBUFF *);

   BLACSCONTEXT *ctxt;
   BLACBUFF *bp=NULL;
   SDRVPTR send;
   MPI_Datatype MatTyp;
   Int length, src, tlda, error, one=1;
   char ttop, tscope, tuplo, tdiag;
   extern BLACBUFF BI_AuxBuff, *BI_ActiveQ;

   MGetConTxt(Mpval(ConTxt), ctxt);
   ttop = F2C_CharTrans(top);
   ttop = Mlowcase(ttop);
   tscope = F2C_CharTrans(scope);
   tscope = Mlowcase(tscope);
   tdiag = F2C_CharTrans(diag);
   tdiag = Mlowcase(tdiag);
   tuplo = F2C_CharTrans(uplo);
   tuplo = Mlowcase(tuplo);
#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_BR, __FILE__, tscope, tuplo, tdiag, Mpval(m),
               Mpval(n), Mpval(lda), 1, Mpaddress(rsrc), Mpaddress(csrc));
#endif
/*
 *  If the user has set the default broadcast topology, use it instead of
 *  BLACS default
 */
#ifdef DefBSTop
   if (ttop == ' ') ttop = DefBSTop;
#endif
   if (Mpval(m) <= Mpval(lda)) tlda = Mpval(lda);
   else tlda = Mpval(m);

   switch(tscope)
   {
   case 'r':
      ctxt->scp = &ctxt->rscp;
      src = Mpval(csrc);
      break;
   case 'c':
      ctxt->scp = &ctxt->cscp;
      src = Mpval(rsrc);
      break;
   case 'a':
      ctxt->scp = &ctxt->ascp;
      src = Mvkpnum(ctxt, Mpval(rsrc), Mpval(csrc));
      break;
   default:
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown scope '%c'",
                  tscope);
   }

   MatTyp = BI_GetMpiTrType(ctxt, tuplo, tdiag, Mpval(m), Mpval(n), tlda,
                            MPI_DOUBLE_COMPLEX, &BI_AuxBuff.N);
/*
 * If using default topology, use MPI native broadcast
 */
   if (ttop == ' ')
   {
      error=MPI_Bcast(A, BI_AuxBuff.N, MatTyp, src, ctxt->scp->comm);
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
#else

   send = BI_Asend;
   MPI_Pack_size(one, MatTyp, ctxt->scp->comm, &length);
   bp = BI_GetBuff(length);
   bp->N = length;
   bp->dtype = MPI_PACKED;
#if ZeroByteTypeBug
   if (MatTyp == MPI_BYTE)
   {
      send = BI_Ssend;
      bp->N = 0;
      bp->dtype = MPI_BYTE;
   }
#endif

#endif

   switch(ttop)
   {
   case 'h':
      error = BI_HypBR(ctxt, bp, send, src);
      if (error == NPOW2) BI_TreeBR(ctxt, bp, send, src, 2);
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
      BI_TreeBR(ctxt, bp, send, src, ttop-47);
      break;
   case 't':
      BI_TreeBR(ctxt, bp, send, src, ctxt->Nb_bs);
      break;
   case 'i':
      BI_IdringBR(ctxt, bp, send, src, 1);
      break;
   case 'd':
      BI_IdringBR(ctxt, bp, send, src, -1);
      break;
   case 's':
      BI_SringBR(ctxt, bp, send, src);
      break;
   case 'm':
      BI_MpathBR(ctxt, bp, send, src, ctxt->Nr_bs);
      break;
   case 'f':
      BI_MpathBR(ctxt, bp, send, src, FULLCON);
      break;
   default :
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown topology '%c'",
                  ttop);
   }

#ifdef MpiBuffGood
   error=BI_MPI_TYPE_FREE(&MatTyp);
   if (BI_ActiveQ) BI_UpdateBuffs(NULL);
#endif
#ifndef MpiBuffGood
   BI_Unpack(ctxt, (BVOID *) A, bp, MatTyp);
   BI_UpdateBuffs(bp);
#endif
}
