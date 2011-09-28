#include "Bdef.h"


#if (INTFACE == C_CALL)
void Csgsum2d(int ConTxt, char *scope, char *top, int m, int n, float *A,
              int lda, int rdest, int cdest)
#else
F_VOID_FUNC sgsum2d_(int *ConTxt, F_CHAR scope, F_CHAR top, int *m, int *n,
                     float *A, int *lda, int *rdest, int *cdest)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Combine sum operation for real rectangular matrices.
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
 *  RDEST   (input) Ptr to int
 *          The process row of the destination of the sum.
 *          If rdest == -1, then result is left on all processes in scope.
 *
 *  CDEST   (input) Ptr to int
 *          The process column of the destination of the sum.
 *          If rdest == -1, then CDEST ignored.
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_ArgCheck(int, int, char *, char, char, char, int, int, int, int,
                    int *, int *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(int);
   int BI_BuffIsFree(BLACBUFF *, int);
   MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *, int, int, int,
                                   MPI_Datatype, int *);
   BLACBUFF *BI_Pack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *, int, int, int,
                                   MPI_Datatype, int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_MringComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, int, VVFUNPTR,
                     int, int);
   void BI_TreeComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, int, VVFUNPTR,
                    int, int);
   void BI_BeComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, int, VVFUNPTR);
   void BI_svvsum(int, char *, char *);
/*
 *  Variable Declarations
 */
   BLACBUFF *bp, *bp2;
   BLACSCONTEXT *ctxt;
   char ttop, tscope;
   int N, length, dest, tlda, trdest, ierr;
   extern BLACBUFF *BI_ActiveQ;
   extern BLACBUFF BI_AuxBuff;

   MGetConTxt(Mpval(ConTxt), ctxt);
   ttop = F2C_CharTrans(top);
   ttop = Mlowcase(ttop);
   tscope = F2C_CharTrans(scope);
   tscope = Mlowcase(tscope);
/*
 *  If the user has set the default combine topology, use it instead of
 *  BLACS default
 */
#ifdef DefCombTop
   if (ttop == ' ') ttop = DefCombTop;
#endif
   if (Mpval(cdest) == -1) trdest = -1;
   else trdest = Mpval(rdest);
#if (BlacsDebugLvl > 0)
   BI_ArgCheck(Mpval(ConTxt), RT_COMB, __FILE__, tscope, 'u', 'u', Mpval(m),
               Mpval(n), Mpval(lda), 1, &trdest, Mpaddress(cdest));
#endif
   if (Mpval(lda) >= Mpval(m)) tlda = Mpval(lda);
   else tlda = Mpval(m);
   switch(tscope)
   {
   case 'r':
      ctxt->scp = &ctxt->rscp;
      if (trdest == -1) dest = -1;
      else dest = Mpval(cdest);
      break;
   case 'c':
      ctxt->scp = &ctxt->cscp;
      dest = trdest;
      break;
   case 'a':
      ctxt->scp = &ctxt->ascp;
      if (trdest == -1) dest = -1;
      else dest = Mvkpnum(ctxt, trdest, Mpval(cdest));
      break;
   default:
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown scope '%c'",
                  tscope);
   }


/*
 * It's not defined how MPI reacts to 0 element reductions, so use BLACS 1-tree
 * topology if we've got one.  Also, we can't use MPI functions if we need to
 * guarantee repeatability.
 */
   if (ttop == ' ')
      if ( (Mpval(m) < 1) || (Mpval(n) < 1) || (ctxt->TopsRepeat) ) ttop = '1';
   N = Mpval(m) * Mpval(n);
   length = N * sizeof(float);
/*
 * If A is contiguous, we can use it as one of the buffers
 */
   if ( (Mpval(m) == tlda) || (Mpval(n) == 1) )
   {
      bp = &BI_AuxBuff;
      bp->Buff = (char *) A;
      bp2 = BI_GetBuff(length);
   }
/*
 * Otherwise, we must allocate both buffers
 */
   else
   {
      bp = BI_GetBuff(length*2);
      bp2 = &BI_AuxBuff;
      bp2->Buff = &bp->Buff[length];
      BI_smvcopy(Mpval(m), Mpval(n), A, tlda, bp->Buff);
   }
   bp->dtype = bp2->dtype = MPI_FLOAT;
   bp->N = bp2->N = N;

   switch(ttop)
   {
   case ' ':         /* use MPI's reduction by default */
      if (dest != -1)
      {
         ierr=MPI_Reduce(bp->Buff, bp2->Buff, bp->N, bp->dtype, MPI_SUM,
                       dest, ctxt->scp->comm);
         if (ctxt->scp->Iam == dest)
	    BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp2->Buff);
      }
      else
      {
         ierr=MPI_Allreduce(bp->Buff, bp2->Buff, bp->N, bp->dtype, MPI_SUM,
		          ctxt->scp->comm);
	 BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp2->Buff);
      }
      if (BI_ActiveQ) BI_UpdateBuffs(NULL);
      return;
      break;
   case 'i':
      BI_MringComb(ctxt, bp, bp2, N, BI_svvsum, dest, 1);
      break;
   case 'd':
      BI_MringComb(ctxt, bp, bp2, N, BI_svvsum, dest, -1);
      break;
   case 's':
      BI_MringComb(ctxt, bp, bp2, N, BI_svvsum, dest, 2);
      break;
   case 'm':
      BI_MringComb(ctxt, bp, bp2, N, BI_svvsum, dest, ctxt->Nr_co);
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
      BI_TreeComb(ctxt, bp, bp2, N, BI_svvsum, dest, ttop-47);
      break;
   case 'f':
      BI_TreeComb(ctxt, bp, bp2, N, BI_svvsum, dest, FULLCON);
      break;
   case 't':
      BI_TreeComb(ctxt, bp, bp2, N, BI_svvsum, dest, ctxt->Nb_co);
      break;
   case 'h':
/*
 *    Use bidirectional exchange if everyone wants answer
 */
      if ( (trdest == -1) && !(ctxt->TopsCohrnt) )
         BI_BeComb(ctxt, bp, bp2, N, BI_svvsum);
      else
         BI_TreeComb(ctxt, bp, bp2, N, BI_svvsum, dest, 2);
      break;
   default :
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown topology '%c'",
                  ttop);
   }

/*
 * If I am selected to receive answer
 */
   if (bp != &BI_AuxBuff)
   {
      if ( (ctxt->scp->Iam == dest) || (dest == -1) )
         BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp->Buff);
      BI_UpdateBuffs(bp);
   }
   else
   {
      if (BI_ActiveQ) BI_UpdateBuffs(NULL);
      BI_BuffIsFree(bp, 1);
   }
}
