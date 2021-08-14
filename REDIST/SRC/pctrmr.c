#include "redist.h"
/** $Id: pctrmr.c,v 1.1.1.1 2000/02/15 18:04:09 susan Exp $
  ------------------------------------------------------------------------

    -- ScaLAPACK routine (version 1.7) --
       Oak Ridge National Laboratory, Univ. of Tennessee, and Univ. of
       California, Berkeley.
       October 31, 1994.

      SUBROUTINE PCTRMR2D(UPLO, DIAG, M, N,
     $                     A, IA, JA, ADESC,
     $                     B, IB, JB, BDESC,
     $                     CTXT)
  ------------------------------------------------------------------------
    Purpose
    =======

    PCTRMR2D copies a submatrix of A on a submatrix of B.
    A and B can have different distributions: they can be on different
    processor grids, they can have different blocksizes, the beginning
    of the area to be copied can be at a different places on A and B.

    The parameters can be confusing when the grids of A and B are
    partially or completly disjoint, in the case a processor calls
    this routines but is either not in the A context or B context, the
    ADESC[CTXT] or BDESC[CTXT] must be equal to -1, to ensure the
    routine recognise this situation.
    To summarize the rule:
    - If a processor is in A context, all parameters related to A must be valid.
    - If a processor is in B context, all parameters related to B must be valid.
    -  ADESC[CTXT] and BDESC[CTXT] must be either valid contexts or equal to -1.
    - M and N must be valid for everyone.
    - other parameters are not examined.

    The submatrix to be copied is assumed to be trapezoidal. So only
    the upper or the lower part will be copied. The other part is
    unchanged.

    Notes
    =====

    A description vector is associated with each 2D block-cyclicly dis-
    tributed matrix.  This vector stores the information required to
    establish the mapping between a matrix entry and its corresponding
    process and memory location.

    In the following comments, the character _ should be read as
    "of the distributed matrix".  Let A be a generic term for any 2D
    block cyclicly distributed matrix.  Its description vector is DESC_A:

   NOTATION        STORED IN      EXPLANATION
   --------------- -------------- --------------------------------------
   DT_A   (global) DESCA( DT_ )   The descriptor type.
   CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                                  the BLACS process grid A is distribu-
                                  ted over. The context itself is glo-
                                  bal, but the handle (the integer
                                  value) may vary.
   M_A    (global) DESCA( M_ )    The number of rows in the distributed
                                  matrix A.
   N_A    (global) DESCA( N_ )    The number of columns in the distri-
                                  buted matrix A.
   MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                                  the rows of A.
   NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                                  the columns of A.
   RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                                  row of the matrix A is distributed.
   CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                                  first column of A is distributed.
   LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                                  array storing the local blocks of the
                                  distributed matrix A.
                                  LLD_A >= MAX(1,LOCp(M_A)).



    Important notice
    ================
     The parameters of the routine have changed in April 1996
     There is a new last argument. It must be a context englobing
     all processors involved in the initial and final distribution.

     Be aware that all processors  included in this
      context must call the redistribution routine.

    Parameters
    ==========

    UPLO     (input) CHARACTER*1.
             On entry, UPLO specifies whether we should copy the upper
             part of the lower part of the defined submatrix:
               UPLO = 'U' or 'u' copy the upper triangular part.
               UPLO = 'L' or 'l' copy the lower triangular part.
             Unchanged on exit.

    DIAG     (input) CHARACTER*1.
             On entry, DIAG specifies whether we should copy the diagonal.
                DIAG = 'U' or 'u' do NOT copy the diagonal of the submatrix.
                DIAG = 'N' or 'n' DO copy the diagonal of the submatrix.
             Unchanged on exit.

    M        (input) INTEGER.
             On entry, M specifies the number of rows of the
             submatrix to be copied.  M must be at least zero.
             Unchanged on exit.

    N        (input) INTEGER.
             On entry, N specifies the number of cols of the submatrix
             to be redistributed.rows of B.  M must be at least zero.
             Unchanged on exit.

    A        (input) COMPLEX
             On entry, the source matrix.
             Unchanged on exit.

    IA,JA    (input) INTEGER
             On entry,the coordinates of the beginning of the submatrix
             of A to copy.
             1 <= IA <= M_A - M + 1,1 <= JA <= N_A - N + 1,
             Unchanged on exit.

    ADESC    (input) A description vector (see Notes above)
             If the current processor is not part of the context of A
             the ADESC[CTXT] must be equal to -1.


    B        (output) COMPLEX
             On entry, the destination matrix.
             The portion corresponding to the defined submatrix are updated.

    IB,JB    (input) INTEGER
             On entry,the coordinates of the beginning of the submatrix
             of B that will be updated.
             1 <= IB <= M_B - M + 1,1 <= JB <= N_B - N + 1,
             Unchanged on exit.

    BDESC    (input) B description vector (see Notes above)
             For processors not part of the context of B
             BDESC[CTXT] must be equal to -1.

    CTXT     (input) a context englobing at least all processors included
                in either A context or B context



   Memory requirement :
   ====================

   for the processors belonging to grid 0, one buffer of size block 0
   and for the processors belonging to grid 1, also one buffer of size
   block 1.

   ============================================================
   Created March 1993 by B. Tourancheau (See sccs for modifications).
   Modifications by Loic PRYLLI 1995
   ============================================================ */
#define static2 static
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2d pctrmr2do_
#define fortran_mr2dnew pctrmr2d_
#elif defined(UpCase)
#define fortran_mr2dnew PCTRMR2D
#define fortran_mr2d PCTRMR2DO
#define ccopy_ CCOPY
#define clacpy_ CLACPY
#else
#define fortran_mr2d pctrmr2do
#define fortran_mr2dnew pctrmr2d
#define ccopy_ ccopy
#define clacpy_ clacpy
#endif
#define Clacpy Cctrlacpy
void  Clacpy();
typedef struct {
  float r, i;
}     complex;
typedef struct {
  Int   desctype;
  Int   ctxt;
  Int   m;
  Int   n;
  Int   nbrow;
  Int   nbcol;
  Int   sprow;
  Int   spcol;
  Int   lda;
}     MDESC;
#define BLOCK_CYCLIC_2D 1
typedef struct {
  Int   gstart;
  Int   len;
}     IDESC;
#define SHIFT(row,sprow,nbrow) ((row)-(sprow)+ ((row) >= (sprow) ? 0 : (nbrow)))
#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)>(B)?(B):(A))
#define DIVUP(a,b) ( ((a)-1) /(b)+1)
#define ROUNDUP(a,b) (DIVUP(a,b)*(b))
#ifdef MALLOCDEBUG
#define malloc mymalloc
#define free myfree
#define realloc myrealloc
#endif
/* Cblacs */
extern void Cblacs_pcoord();
extern Int Cblacs_pnum();
extern void Csetpvmtids();
extern void Cblacs_get();
extern void Cblacs_pinfo();
extern void Cblacs_gridinfo();
extern void Cblacs_gridinit();
extern void Cblacs_exit();
extern void Cblacs_gridexit();
extern void Cblacs_setup();
extern void Cigebs2d();
extern void Cigebr2d();
extern void Cigesd2d();
extern void Cigerv2d();
extern void Cigsum2d();
extern void Cigamn2d();
extern void Cigamx2d();
extern void Ccgesd2d();
extern void Ccgerv2d();
/* lapack */
void  clacpy_();
/* aux fonctions */
extern Int localindice();
extern void *mr2d_malloc();
extern Int ppcm();
extern Int localsize();
extern Int memoryblocksize();
extern Int changeorigin();
extern void paramcheck();
/* tools and others function */
#define scanD0 ctrscanD0
#define dispmat ctrdispmat
#define setmemory ctrsetmemory
#define freememory ctrfreememory
#define scan_intervals ctrscan_intervals
extern void scanD0();
extern void dispmat();
extern void setmemory();
extern void freememory();
extern Int scan_intervals();
extern void Cpctrmr2do();
extern void Cpctrmr2d();
/* some defines for Cpctrmr2do */
#define SENDBUFF 0
#define RECVBUFF 1
#define SIZEBUFF 2
#if 0
#define DEBUG
#endif
#ifndef DEBUG
#define NDEBUG
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define DESCLEN 9
void 
fortran_mr2d(char *uplo, char *diag, Int *m, Int *n, complex *A, Int *ia, Int *ja, Int desc_A[DESCLEN],
	     complex *B, Int *ib, Int *jb, Int desc_B[DESCLEN])
{
  Cpctrmr2do(uplo, diag, *m, *n, A, *ia, *ja, (MDESC *) desc_A,
	     B, *ib, *jb, (MDESC *) desc_B);
  return;
}
void 
fortran_mr2dnew(char *uplo, char *diag, Int *m, Int *n, complex *A, Int *ia, Int *ja, Int desc_A[DESCLEN],
		complex *B, Int *ib, Int *jb, Int desc_B[DESCLEN], Int *gcontext)
{
  Cpctrmr2d(uplo, diag, *m, *n, A, *ia, *ja, (MDESC *) desc_A,
	    B, *ib, *jb, (MDESC *) desc_B, *gcontext);
  return;
}
static2 void init_chenille();
static2 Int inter_len();
static2 Int block2buff();
static2 void buff2block();
static2 void gridreshape();
void
Cpctrmr2do(uplo, diag, m, n,
	   ptrmyblock, ia, ja, ma,
	   ptrmynewblock, ib, jb, mb)
  char *uplo, *diag;
  complex *ptrmyblock, *ptrmynewblock;
/* pointers to the memory location of the matrix and the redistributed matrix */
  MDESC *ma;
  MDESC *mb;
  Int   ia, ja, ib, jb, m, n;
{
  Int   dummy, nprocs;
  Int   gcontext;
  /* first we initialize a global grid which serve as a reference to
   * communicate from grid a to grid b */
  Cblacs_pinfo(&dummy, &nprocs);
  Cblacs_get(0, 0, &gcontext);
  Cblacs_gridinit(&gcontext, "R", (Int)1, nprocs);
  Cpctrmr2d(uplo, diag, m, n, ptrmyblock, ia, ja, ma,
	    ptrmynewblock, ib, jb, mb, gcontext);
  Cblacs_gridexit(gcontext);
}
#define NBPARAM 20	/* p0,q0,p1,q1, puis ma,na,mba,nba,rowa,cola puis
			 * idem B puis ia,ja puis ib,jb */
#define MAGIC_MAX 100000000
void
Cpctrmr2d(uplo, diag, m, n,
	  ptrmyblock, ia, ja, ma,
	  ptrmynewblock, ib, jb, mb, globcontext)
  char *uplo, *diag;
  complex *ptrmyblock, *ptrmynewblock;
/* pointers to the memory location of the matrix and the redistributed matrix */
  MDESC *ma;
  MDESC *mb;
  Int   ia, ja, ib, jb, m, n, globcontext;
{
  complex *ptrsendbuff, *ptrrecvbuff, *ptrNULL = 0;
  complex *recvptr;
  MDESC newa, newb;
  Int  *proc0, *proc1, *param;
  Int   mypnum, myprow0, mypcol0, myprow1, mypcol1, nprocs;
  Int   i, j;
  Int   nprow, npcol, gcontext;
  Int   recvsize, sendsize;
  IDESC *h_inter;	/* to store the horizontal intersections */
  IDESC *v_inter;	/* to store the vertical intersections */
  Int   hinter_nb, vinter_nb;	/* number of intrsections in both directions */
  Int   dummy;
  Int   p0, q0, p1, q1;
  Int  *ra, *ca;
  /* end of variables */
  /* To simplify further calcul we change the matrix indexation from
   * 1..m,1..n (fortran) to 0..m-1,0..n-1 */
  if (m == 0 || n == 0)
    return;
  ia -= 1;
  ja -= 1;
  ib -= 1;
  jb -= 1;
  Cblacs_gridinfo(globcontext, &nprow, &npcol, &dummy, &mypnum);
  gcontext = globcontext;
  nprocs = nprow * npcol;
  /* if the global context that is given to us has not the shape of a line
   * (nprow != 1), create a new context.  TODO: to be optimal, we should
   * avoid this because it is an uncessary synchronisation */
  if (nprow != 1) {
    gridreshape(&gcontext);
    Cblacs_gridinfo(gcontext, &dummy, &dummy, &dummy, &mypnum);
  }
  Cblacs_gridinfo(ma->ctxt, &p0, &q0, &myprow0, &mypcol0);
  /* compatibility T3D, must check myprow  and mypcol are within bounds */
  if (myprow0 >= p0 || mypcol0 >= q0)
    myprow0 = mypcol0 = -1;
  assert((myprow0 < p0 && mypcol0 < q0) || (myprow0 == -1 && mypcol0 == -1));
  Cblacs_gridinfo(mb->ctxt, &p1, &q1, &myprow1, &mypcol1);
  if (myprow1 >= p1 || mypcol1 >= q1)
    myprow1 = mypcol1 = -1;
  assert((myprow1 < p1 && mypcol1 < q1) || (myprow1 == -1 && mypcol1 == -1));
  /* exchange the missing parameters among the processors: shape of grids and
   * location of the processors */
  param = (Int *) mr2d_malloc(3 * (nprocs * 2 + NBPARAM) * sizeof(Int));
  ra = param + nprocs * 2 + NBPARAM;
  ca = param + (nprocs * 2 + NBPARAM) * 2;
  for (i = 0; i < nprocs * 2 + NBPARAM; i++)
    param[i] = MAGIC_MAX;
  proc0 = param + NBPARAM;
  proc1 = param + NBPARAM + nprocs;
  /* we calulate proc0 and proc1 that will give the number of a proc in
   * respectively a or b in the global context */
  if (myprow0 >= 0) {
    proc0[myprow0 * q0 + mypcol0] = mypnum;
    param[0] = p0;
    param[1] = q0;
    param[4] = ma->m;
    param[5] = ma->n;
    param[6] = ma->nbrow;
    param[7] = ma->nbcol;
    param[8] = ma->sprow;
    param[9] = ma->spcol;
    param[10] = ia;
    param[11] = ja;
  }
  if (myprow1 >= 0) {
    proc1[myprow1 * q1 + mypcol1] = mypnum;
    param[2] = p1;
    param[3] = q1;
    param[12] = mb->m;
    param[13] = mb->n;
    param[14] = mb->nbrow;
    param[15] = mb->nbcol;
    param[16] = mb->sprow;
    param[17] = mb->spcol;
    param[18] = ib;
    param[19] = jb;
  }
  Cigamn2d(gcontext, "All", "H", 2 * nprocs + NBPARAM, (Int)1, param, 2 * nprocs + NBPARAM,
	   ra, ca, 2 * nprocs + NBPARAM, (Int)-1, (Int)-1);
  newa = *ma;
  newb = *mb;
  ma = &newa;
  mb = &newb;
  if (myprow0 == -1) {
    p0 = param[0];
    q0 = param[1];
    ma->m = param[4];
    ma->n = param[5];
    ma->nbrow = param[6];
    ma->nbcol = param[7];
    ma->sprow = param[8];
    ma->spcol = param[9];
    ia = param[10];
    ja = param[11];
  }
  if (myprow1 == -1) {
    p1 = param[2];
    q1 = param[3];
    mb->m = param[12];
    mb->n = param[13];
    mb->nbrow = param[14];
    mb->nbcol = param[15];
    mb->sprow = param[16];
    mb->spcol = param[17];
    ib = param[18];
    jb = param[19];
  }
  for (i = 0; i < NBPARAM; i++) {
    if (param[i] == MAGIC_MAX) {
      fprintf(stderr, "xxGEMR2D:something wrong in the parameters\n");
      exit(1);
    }
  }
#ifndef NDEBUG
  for (i = 0; i < p0 * q0; i++)
    assert(proc0[i] >= 0 && proc0[i] < nprocs);
  for (i = 0; i < p1 * q1; i++)
    assert(proc1[i] >= 0 && proc1[i] < nprocs);
#endif
  /* check the validity of the parameters */
  paramcheck(ma, ia, ja, m, n, p0, q0, gcontext);
  paramcheck(mb, ib, jb, m, n, p1, q1, gcontext);
  /* we change the problem so that ia < a->nbrow ... andia + m = a->m ... */
  {
    Int   decal;
    ia = changeorigin(myprow0, ma->sprow, p0,
		      ma->nbrow, ia, &decal, &ma->sprow);
    ptrmyblock += decal;
    ja = changeorigin(mypcol0, ma->spcol, q0,
		      ma->nbcol, ja, &decal, &ma->spcol);
    ptrmyblock += decal * ma->lda;
    ma->m = ia + m;
    ma->n = ja + n;
    ib = changeorigin(myprow1, mb->sprow, p1,
		      mb->nbrow, ib, &decal, &mb->sprow);
    ptrmynewblock += decal;
    jb = changeorigin(mypcol1, mb->spcol, q1,
		      mb->nbcol, jb, &decal, &mb->spcol);
    ptrmynewblock += decal * mb->lda;
    mb->m = ib + m;
    mb->n = jb + n;
    if (p0 == 1)
      ma->nbrow = ma->m;
    if (q0 == 1)
      ma->nbcol = ma->n;
    if (p1 == 1)
      mb->nbrow = mb->m;
    if (q1 == 1)
      mb->nbcol = mb->n;
#ifndef NDEBUG
    paramcheck(ma, ia, ja, m, n, p0, q0, gcontext);
    paramcheck(mb, ib, jb, m, n, p1, q1, gcontext);
#endif
  }
  /* We compute the size of the memory buffer ( we choose the worst case,
   * when the buffer sizes == the memory block sizes). */
  if (myprow0 >= 0 && mypcol0 >= 0) {
    /* Initialize pointer variables */
    setmemory(&ptrsendbuff, memoryblocksize(ma));
  };	/* if (mypnum < p0 * q0) */
  if (myprow1 >= 0 && mypcol1 >= 0) {
    /* Initialize pointer variables */
    setmemory(&ptrrecvbuff, memoryblocksize(mb));
  };	/* if (mypnum < p1 * q1) */
  /* allocing room for the tabs, alloc for the worst case,local_n or local_m
   * intervals, in fact the worst case should be less, perhaps half that,I
   * should think of that one day. */
  h_inter = (IDESC *) mr2d_malloc(DIVUP(ma->n, q0 * ma->nbcol) *
				  ma->nbcol * sizeof(IDESC));
  v_inter = (IDESC *) mr2d_malloc(DIVUP(ma->m, p0 * ma->nbrow)
				  * ma->nbrow * sizeof(IDESC));
  /* We go for the scanning of indices. For each processor including mypnum,
   * we fill the sendbuff buffer (scanD0(SENDBUFF)) and when it is done send
   * it. Then for each processor, we compute the size of message to be
   * receive scanD0(SIZEBUFF)), post a receive and then allocate the elements
   * of recvbuff the right place (scanD)(RECVBUFF)) */
  recvptr = ptrrecvbuff;
  {
    Int   tot, myrang, step, sens;
    Int  *sender, *recver;
    Int   mesending, merecving;
    tot = max(p0 * q0, p1 * q1);
    init_chenille(mypnum, nprocs, p0 * q0, proc0, p1 * q1, proc1,
		  &sender, &recver, &myrang);
    if (myrang == -1)
      goto after_comm;
    mesending = myprow0 >= 0;
    assert(sender[myrang] >= 0 || !mesending);
    assert(!mesending || proc0[sender[myrang]] == mypnum);
    merecving = myprow1 >= 0;
    assert(recver[myrang] >= 0 || !merecving);
    assert(!merecving || proc1[recver[myrang]] == mypnum);
    step = tot - 1 - myrang;
    do {
      for (sens = 0; sens < 2; sens++) {
	/* be careful here, when we communicating with ourselves, we must
	 * send first (myrang > step == 0) */
	if (mesending && recver[step] >= 0 &&
	    (sens == 0)) {
	  i = recver[step] / q1;
	  j = recver[step] % q1;
	  vinter_nb = scan_intervals('r', ia, ib, m, ma, mb, p0, p1, myprow0, i,
				     v_inter);
	  hinter_nb = scan_intervals('c', ja, jb, n, ma, mb, q0, q1, mypcol0, j,
				     h_inter);
	  scanD0(uplo, diag, SENDBUFF, ptrsendbuff, &sendsize,
		 m, n, ma, ia, ja, p0, q0, mb, ib, jb, p1, q1,
		 v_inter, vinter_nb, h_inter, hinter_nb,
		 ptrmyblock);
	}	/* if (mesending...) { */
	if (mesending && recver[step] >= 0 &&
	    (sens == myrang > step)) {
	  i = recver[step] / q1;
	  j = recver[step] % q1;
	  if (sendsize > 0
	      && (step != myrang || !merecving)
		) {
	    Ccgesd2d(gcontext, sendsize, (Int)1, ptrsendbuff, sendsize,
		     (Int)0, proc1[i * q1 + j]);
	  }	/* sendsize > 0 */
	}	/* if (mesending ... */
	if (merecving && sender[step] >= 0 &&
	    (sens == myrang <= step)) {
	  i = sender[step] / q0;
	  j = sender[step] % q0;
	  vinter_nb = scan_intervals('r', ib, ia, m, mb, ma, p1, p0, myprow1, i,
				     v_inter);
	  hinter_nb = scan_intervals('c', jb, ja, n, mb, ma, q1, q0, mypcol1, j,
				     h_inter);
	  scanD0(uplo, diag, SIZEBUFF, ptrNULL, &recvsize,
		 m, n, ma, ia, ja, p0, q0, mb, ib, jb, p1, q1,
		 v_inter, vinter_nb, h_inter, hinter_nb, ptrNULL);
	  if (recvsize > 0) {
	    if (step == myrang && mesending) {
	      Clacpy(recvsize, 1,
		     ptrsendbuff, recvsize,
		     ptrrecvbuff, recvsize);
	    } else {
	      Ccgerv2d(gcontext, recvsize, (Int)1, ptrrecvbuff, recvsize,
		       (Int)0, proc0[i * q0 + j]);
	    }
	  }	/* recvsize > 0 */
	}	/* if (merecving ...) */
	if (merecving && sender[step] >= 0 && sens == 1) {
	  scanD0(uplo, diag, RECVBUFF, ptrrecvbuff, &recvsize,
		 m, n, ma, ia, ja, p0, q0, mb, ib, jb, p1, q1,
		 v_inter, vinter_nb, h_inter, hinter_nb, ptrmynewblock);
	}	/* if (merecving...)  */
      }	/* for (sens = 0) */
      step -= 1;
      if (step < 0)
	step = tot - 1;
    } while (step != tot - 1 - myrang);
after_comm:
    free(sender);
  }	/* { int tot,nr,ns ...} */
  /* don't forget to clean up things! */
  if (myprow1 >= 0 && mypcol1 >= 0) {
    freememory((char *) ptrrecvbuff);
  };
  if (myprow0 >= 0 && mypcol0 >= 0) {
    freememory((char *) ptrsendbuff);
  };
  if (nprow != 1)
    Cblacs_gridexit(gcontext);
  free(v_inter);
  free(h_inter);
  free(param);
}/* distrib */
static2 void 
init_chenille(Int mypnum, Int nprocs, Int n0, Int *proc0, Int n1, Int *proc1, Int **psend, Int **precv, Int *myrang)
{
  Int   ns, nr, i, tot;
  Int  *sender, *recver, *g0, *g1;
  tot = max(n0, n1);
  sender = (Int *) mr2d_malloc((nprocs + tot) * sizeof(Int) * 2);
  recver = sender + tot;
  *psend = sender;
  *precv = recver;
  g0 = recver + tot;
  g1 = g0 + nprocs;
  for (i = 0; i < nprocs; i++) {
    g0[i] = -1;
    g1[i] = -1;
  }
  for (i = 0; i < tot; i++) {
    sender[i] = -1;
    recver[i] = -1;
  }
  for (i = 0; i < n0; i++)
    g0[proc0[i]] = i;
  for (i = 0; i < n1; i++)
    g1[proc1[i]] = i;
  ns = 0;
  nr = 0;
  *myrang = -1;
  for (i = 0; i < nprocs; i++)
    if (g0[i] >= 0 && g1[i] >= 0) {
      if (i == mypnum)
	*myrang = nr;
      sender[ns] = g0[i];
      ns += 1;
      recver[nr] = g1[i];
      nr += 1;
      assert(ns <= n0 && nr <= n1 && nr == ns);
    }
  for (i = 0; i < nprocs; i++)
    if (g0[i] >= 0 && g1[i] < 0) {
      if (i == mypnum)
	*myrang = ns;
      sender[ns] = g0[i];
      ns += 1;
      assert(ns <= n0);
    }
  for (i = 0; i < nprocs; i++)
    if (g1[i] >= 0 && g0[i] < 0) {
      if (i == mypnum)
	*myrang = nr;
      recver[nr] = g1[i];
      nr += 1;
      assert(nr <= n1);
    }
}
void 
Clacpy(Int m, Int n, complex *a, Int lda, complex *b, Int ldb)
{
  Int   i, j;
  lda -= m;
  ldb -= m;
  assert(lda >= 0 && ldb >= 0);
  for (j = 0; j < n; j++) {
    for (i = 0; i < m; i++)
      *b++ = *a++;
    b += ldb;
    a += lda;
  }
}
static2 void 
gridreshape(Int *ctxtp)
{
  Int   ori, final;	/* original context, and new context created, with
			 * line form */
  Int   nprow, npcol, myrow, mycol;
  Int  *usermap;
  Int   i, j;
  ori = *ctxtp;
  Cblacs_gridinfo(ori, &nprow, &npcol, &myrow, &mycol);
  usermap = mr2d_malloc(sizeof(Int) * nprow * npcol);
  for (i = 0; i < nprow; i++)
    for (j = 0; j < npcol; j++) {
      usermap[i + j * nprow] = Cblacs_pnum(ori, i, j);
    }
  /* Cblacs_get(0, 0, &final); */
  Cblacs_get(ori, (Int)10, &final);
  Cblacs_gridmap(&final, usermap, (Int)1, (Int)1, nprow * npcol);
  *ctxtp = final;
  free(usermap);
}
