#include "Bdef.h"



#if (INTFACE == C_CALL)
void Csgamn2d(Int ConTxt, char *scope, char *top, Int m, Int n, float *A,
              Int lda, Int *rA, Int *cA, Int ldia, Int rdest, Int cdest)
#else
F_VOID_FUNC sgamn2d_(Int *ConTxt, F_CHAR scope, F_CHAR top, Int *m, Int *n,
                     float *A, Int *lda, Int *rA, Int *cA, Int *ldia,
                     Int *rdest, Int *cdest)
#endif
/*
 *  -- V1.1 BLACS routine --
 *  University of Tennessee, May 1, 1996
 *  Written by Clint Whaley.
 *
 *  Purpose
 *  =======
 *  Combine amn operation for real rectangular matrices.
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
 *  M       (input) Ptr to Int
 *          The number of rows of the matrix A.  M >= 0.
 *
 *  N       (input) Ptr to Int
 *          The number of columns of the matrix A.  N >= 0.
 *
 *  A       (output) Ptr to real two dimensional array
 *          The m by n matrix A.  Fortran77 (column-major) storage
 *          assumed.
 *
 *  LDA     (input) Ptr to Int
 *          The leading dimension of the array A.  LDA >= M.
 *
 *  RA      (output) Integer Array, dimension (LDIA, N)
 *          Contains process row that the amn of each element
 *          of A was found on: i.e., rA(1,2) contains the process
 *          row that the amn of A(1,2) was found on.
 *          Values are left on process {rdest, cdest} only, others
 *          may be modified, but not left with interesting data.
 *          If rdest == -1, then result is left on all processes in scope.
 *          If LDIA == -1, this array is not accessed, and need not exist.
 *
 *  CA      (output) Integer Array, dimension (LDIA, N)
 *          Contains process column that the amn of each element
 *          of A was found on: i.e., cA(1,2) contains the process
 *          column that the max/min of A(1,2) was found on.
 *          Values are left on process {rdest, cdest} only, others
 *          may be modified, but not left with interesting data.
 *          If rdest == -1, then result is left on all processes in scope.
 *          If LDIA == -1, this array is not accessed, and need not exist.
 *
 *  LDIA    (input) Ptr to Int
 *          If (LDIA == -1), then the arrays RA and CA are not accessed.
 *          ELSE leading dimension of the arrays RA and CA.  LDIA >= M.
 *
 *  RDEST   (input) Ptr to Int
 *          The process row of the destination of the amn.
 *          If rdest == -1, then result is left on all processes in scope.
 *
 *  CDEST   (input) Ptr to Int
 *          The process column of the destination of the amn.
 *          If rdest == -1, then CDEST ignored.
 *
 * ------------------------------------------------------------------------
 */
{
   void BI_ArgCheck(Int, Int, char *, char, char, char, Int, Int, Int, Int,
                    Int *, Int *);
   void BI_UpdateBuffs(BLACBUFF *);
   BLACBUFF *BI_GetBuff(Int);
   Int BI_BuffIsFree(BLACBUFF *, Int);
   MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *, Int, Int, Int,
                                   MPI_Datatype, Int *);
   BLACBUFF *BI_Pack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   MPI_Datatype BI_GetMpiGeType(BLACSCONTEXT *, Int, Int, Int,
                                   MPI_Datatype, Int *);
   void BI_Unpack(BLACSCONTEXT *, BVOID *, BLACBUFF *, MPI_Datatype);
   void BI_MringComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, Int, VVFUNPTR,
                     Int, Int);
   void BI_TreeComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, Int, VVFUNPTR,
                    Int, Int);
   void BI_BeComb(BLACSCONTEXT *, BLACBUFF *, BLACBUFF *, Int, VVFUNPTR);
   void BI_svvamn(Int, char *, char *);
   void BI_svvamn2(Int, char *, char *);
   void BI_sMPI_amn(void *, void *, MpiInt *, MPI_Datatype *);
   void BI_sMPI_amn2(void *, void *, MpiInt *, MPI_Datatype *);
/*
 *  Variable Declarations
 */
   VVFUNPTR vvop;
   BLACBUFF *bp, *bp2;
   BLACSCONTEXT *ctxt;
   char ttop, tscope;
   Int i, j, N, dest, idist, length, tlda, tldia, trdest, ierr;
   MpiInt len[2];
   MPI_Aint disp[2];
   MPI_Datatype dtypes[2];
   MPI_Op BlacComb;
   MPI_Datatype IntTyp, MyType;
   BI_DistType *dist, mydist;
   extern BLACBUFF *BI_ActiveQ;
   extern BLACBUFF BI_AuxBuff;

   MPI_Type_match_size(MPI_TYPECLASS_INTEGER, sizeof(Int), &IntTyp);

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
   if (Mpval(ldia) < Mpval(m))
   {
      if (Mpval(ldia) != -1)
         BI_BlacsWarn(Mpval(ConTxt), __LINE__, __FILE__,
                      "LDIA too small (LDIA=%d, but M=%d)", Mpval(ldia),
                      Mpval(m));
   }
#endif
   if (Mpval(lda) >= Mpval(m)) tlda = Mpval(lda);
   else tlda = Mpval(m);
   if (Mpval(ldia) < Mpval(m)) tldia = Mpval(m);
   else tldia = Mpval(ldia);
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
 * topology if we've got one
 */
   if (ttop == ' ')
      if ( (Mpval(m) < 1) || (Mpval(n) < 1) || (ctxt->TopsRepeat) ) ttop = '1';
   N = Mpval(m) * Mpval(n);
/*
 * If process who has amn is to be communicated, must set up distance
 * vector after value vector
 */
   if (Mpval(ldia) != -1)
   {
      vvop = BI_svvamn;
      length = N * sizeof(float);
      i = length % sizeof(BI_DistType);  /* ensure dist vec aligned correctly */
      if (i) length += sizeof(BI_DistType) - i;
      idist = length;
      length += N * sizeof(BI_DistType);
/*
 *    For performance, insist second buffer is at least 8-byte aligned
 */
      j = 8;
      if (sizeof(float) > j) j = sizeof(float);
      i = length % j;
      if (i) length += j - i;
      i = 2 * length;

      bp = BI_GetBuff(i);
      bp2 = &BI_AuxBuff;
      bp2->Buff = &bp->Buff[length];
      BI_smvcopy(Mpval(m), Mpval(n), A, tlda, bp->Buff);
/*
 *    Fill in distance vector
 */
      if (dest == -1) mydist = ctxt->scp->Iam;
      else mydist = (ctxt->scp->Np + ctxt->scp->Iam - dest) % ctxt->scp->Np;
      dist = (BI_DistType *) &bp->Buff[idist];
      for (i=0; i < N; i++) dist[i] = mydist;

/*
 *    Create the MPI datatype holding both user's buffer and distance vector
 */
      len[0] = len[1] = N;
      disp[0] = 0;
      disp[1] = idist;
      dtypes[0] = MPI_FLOAT;
      dtypes[1] = BI_MpiDistType;
#ifdef ZeroByteTypeBug
      if (N > 0)
      {
#endif
      i = 2;
      ierr=MPI_Type_create_struct(i, len, disp, dtypes, &MyType);
      ierr=MPI_Type_commit(&MyType);
      bp->N = bp2->N = 1;
      bp->dtype = bp2->dtype = MyType;
#ifdef ZeroByteTypeBug
      }
      else
      {
         bp->N = bp2->N = 0;
         bp->dtype = bp2->dtype = IntTyp;
      }
#endif
   }
   else
   {
      vvop = BI_svvamn2;
      length = N * sizeof(float);
/*
 *    If A is contiguous, we can use it as one of our buffers
 */
      if ( (Mpval(m) == tlda) || (Mpval(n) == 1) )
      {
         bp = &BI_AuxBuff;
         bp->Buff = (char *) A;
         bp2 = BI_GetBuff(length);
      }
      else
      {
         bp = BI_GetBuff(length*2);
         bp2 = &BI_AuxBuff;
         bp2->Buff = &bp->Buff[length];
         BI_smvcopy(Mpval(m), Mpval(n), A, tlda, bp->Buff);
      }
      bp->N = bp2->N = N;
      bp->dtype = bp2->dtype = MPI_FLOAT;
   }

   switch(ttop)
   {
   case ' ':         /* use MPI's reduction by default */
      i = 1;
      if (Mpval(ldia) == -1)
      {
         ierr=MPI_Op_create(BI_sMPI_amn2, i, &BlacComb);
      }
      else
      {
         ierr=MPI_Op_create(BI_sMPI_amn, i, &BlacComb);
         BI_AuxBuff.Len = N;  /* set this up for the MPI OP wrappers */
      }

      if (trdest != -1)
      {
         ierr=MPI_Reduce(bp->Buff, bp2->Buff, bp->N, bp->dtype, BlacComb, dest,
	 	       ctxt->scp->comm);
         if (ctxt->scp->Iam == dest)
	 {
	    BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp2->Buff);
	    if (Mpval(ldia) != -1)
               BI_TransDist(ctxt, tscope, Mpval(m), Mpval(n), rA, cA, tldia,
                            (BI_DistType *) &bp2->Buff[idist],
			    trdest, Mpval(cdest));
	 }
      }
      else
      {
         ierr=MPI_Allreduce(bp->Buff, bp2->Buff, bp->N, bp->dtype, BlacComb,
		          ctxt->scp->comm);
	 BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp2->Buff);
         if (Mpval(ldia) != -1)
            BI_TransDist(ctxt, tscope, Mpval(m), Mpval(n), rA, cA, tldia,
                         (BI_DistType *) &bp2->Buff[idist],
                         trdest, Mpval(cdest));
      }
      ierr=MPI_Op_free(&BlacComb);
      if (Mpval(ldia) != -1)
#ifdef ZeroByteTypeBug
         if (N > 0)
#endif
         ierr=BI_MPI_TYPE_FREE(&MyType);
      if (BI_ActiveQ) BI_UpdateBuffs(NULL);
      return;
      break;
   case 'i':
      BI_MringComb(ctxt, bp, bp2, N, vvop, dest, 1);
      break;
   case 'd':
      BI_MringComb(ctxt, bp, bp2, N, vvop, dest, -1);
      break;
   case 's':
      BI_MringComb(ctxt, bp, bp2, N, vvop, dest, 2);
      break;
   case 'm':
      BI_MringComb(ctxt, bp, bp2, N, vvop, dest, ctxt->Nr_co);
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
      BI_TreeComb(ctxt, bp, bp2, N, vvop, dest, ttop-47);
      break;
   case 'f':
      BI_TreeComb(ctxt, bp, bp2, N, vvop, dest, FULLCON);
      break;
   case 't':
      BI_TreeComb(ctxt, bp, bp2, N, vvop, dest, ctxt->Nb_co);
      break;
   case 'h':
/*
 *    Use bidirectional exchange if everyone wants answer
 */
      if ( (trdest == -1) && !(ctxt->TopsCohrnt) )
         BI_BeComb(ctxt, bp, bp2, N, vvop);
      else
         BI_TreeComb(ctxt, bp, bp2, N, vvop, dest, 2);
      break;
   default :
      BI_BlacsErr(Mpval(ConTxt), __LINE__, __FILE__, "Unknown topology '%c'",
                  ttop);
   }

   if (Mpval(ldia) != -1)
#ifdef ZeroByteTypeBug
      if (N > 0)
#endif
      ierr=BI_MPI_TYPE_FREE(&MyType);
/*
 * If I am selected to receive answer
 */
   if ( (ctxt->scp->Iam == dest) || (dest == -1) )
   {
/*
 *    Translate the distances stored in the latter part of bp->Buff into
 *    process grid coordinates, and output these coordinates in the
 *    arrays rA and cA.
 */
      if (Mpval(ldia) != -1)
         BI_TransDist(ctxt, tscope, Mpval(m), Mpval(n), rA, cA, tldia,
                      dist, trdest, Mpval(cdest));
/*
 *    Unpack the amn array
 */
      if (bp != &BI_AuxBuff) BI_svmcopy(Mpval(m), Mpval(n), A, tlda, bp->Buff);
   }
}
