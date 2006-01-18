#include "redist.h"
/* $Id: pctrmrdrv.c,v 1.1.1.1 2000/02/15 18:04:11 susan Exp $
 * 
 * pctrmrdrv.c :
 * 
 * 
 * PURPOSE:
 * 
 * this driver is testing the PCTRMR2D routine. It calls it to obtain a new
 * scattered block data decomposition of a distributed COMPLEX (block
 * scattered) matrix. Then it calls PCTRMR2D for the inverse redistribution
 * and checks the results with the initial data.
 * 
 * Data are going from a Block Scattered nbrow0 x nbcol0 decomposition on the
 * processor grid p0 x q0, to data distributed in a BS nbrow1 x nbcol1 on the
 * processor grid p1 x q1, then back to the BS nbrow0 x nbcol0 decomposition
 * on the processor grid p0 x q0.
 * 
 * See pctrmr.c file for detailed info on the PCTRMR2D function.
 * 
 * 
 * The testing parameters are read from the file TRMR2D.dat, see the file in the
 * distribution to have an example.
 * 
 * created by Bernard Tourancheau in April 1994.
 * 
 * modifications : see sccs history
 * 
 * ===================================
 * 
 * 
 * NOTE :
 * 
 * - the matrix elements are COMPLEX
 * 
 * - memory requirements : this procedure requires approximately 3 times the
 * memory space of the initial data block in grid 0 (initial block, copy for
 * test and second redistribution result) and 1 time the memory space of the
 * result data block in grid 1. with  the element size = sizeof(complex)
 * bytes,
 * 
 * 
 * - use the procedures of the files:
 * 
 * pctrmr.o pctrmr2.o pctrmraux.o
 * 
 * 
 * ======================================
 * 
 * WARNING ASSUMPTIONS :
 * 
 * 
 * ========================================
 * 
 * 
 * Planned changes:
 * 
 * 
 * 
 * ========================================= */
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
  int   desctype;
  int   ctxt;
  int   m;
  int   n;
  int   nbrow;
  int   nbcol;
  int   sprow;
  int   spcol;
  int   lda;
}     MDESC;
#define BLOCK_CYCLIC_2D 1
typedef struct {
  int   gstart;
  int   len;
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
extern int Cblacs_pnum();
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
extern int localindice();
extern void *mr2d_malloc();
extern int ppcm();
extern int localsize();
extern int memoryblocksize();
extern int changeorigin();
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
extern int scan_intervals();
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
#include <string.h>
#include <ctype.h>
#include <assert.h>
/* initblock: intialize the local part of a matrix with random data (well,
 * not very random) */
static2 void
initblock(block, m, n)
  complex *block;
  int   m, n;
{
  complex *pdata;
  int   i;
  pdata = block;
  for (i = 0; i < m * n; i++, pdata++) {
    (*pdata).r = i;
  };
}
/* getparam:read from a file a list of integer parameters, the end of the
 * parameters to read is given by a NULL at the end of the args list */
#ifdef __STDC__
#include <stdarg.h>
static void 
getparam(FILE * f,...)
{
#else
#include <varargs.h>
static void 
getparam(va_alist)
va_dcl
{
  FILE *f;
#endif
  va_list ap;
  int   i;
  static int nbline;
  char *ptr, *next;
  int  *var;
  static char buffer[200];
#ifdef __STDC__
  va_start(ap, f);
#else
  va_start(ap);
  f = va_arg(ap, FILE *);
#endif
  do {
    next = fgets(buffer, 200, f);
    if (next == NULL) {
      fprintf(stderr, "bad configuration driver file:after line %d\n", nbline);
      exit(1);
    }
    nbline += 1;
  } while (buffer[0] == '#');
  ptr = buffer;
  var = va_arg(ap, int *);
  while (var != NULL) {
    *var = strtol(ptr, &next, 10);
    if (ptr == next) {
      fprintf(stderr, "bad configuration driver file:error line %d\n", nbline);
      exit(1);
    }
    ptr = next;
    var = va_arg(ap, int *);
  }
  va_end(ap);
}
void 
initforpvm(argc, argv)
  int   argc;
  char *argv[];
{
  int   pnum, nproc;
  Cblacs_pinfo(&pnum, &nproc);
  if (nproc < 1) {	/* we are with PVM */
    if (pnum == 0) {
      if (argc < 2) {
	fprintf(stderr, "usage with PVM:xctrmr nbproc\n\
\t where nbproc is the number of nodes to initialize\n");
	exit(1);
      }
      nproc = atoi(argv[1]);
    }
    Cblacs_setup(&pnum, &nproc);
  }
}
int
main(argc, argv)
  int   argc;
  char *argv[];
{
  /* We initialize the data-block on the current processor, then redistribute
   * it, and perform the inverse redistribution  to compare the local memory
   * with the initial one. */
  /* Data file */
  FILE *fp;
  int   nbre, nbremax;
  /* Data distribution 0 parameters */
  int   p0,	/* # of rows in the processor grid */
        q0;	/* # of columns in the processor grid */
  /* Data distribution 1 parameters */
  int   p1, q1;
  /* # of parameter to be read on the keyboard */
#define nbparameter 24
  /* General variables */
  int   blocksize0;
  int   mypnum, nprocs;
  int   parameters[nbparameter], nberrors;
  int   i;
  int   ia, ja, ib, jb, m, n;
  int   gcontext, context0, context1;
  int   myprow1, myprow0, mypcol0, mypcol1;
  int   dummy;
  MDESC ma, mb;
  char *uplo, *diag;
  complex *ptrmyblock, *ptrsavemyblock, *ptrmyblockcopy, *ptrmyblockvide;
#ifdef UsingMpiBlacs
   MPI_Init(&argc, &argv);
#endif
  setvbuf(stdout, NULL, _IOLBF, 0);
  setvbuf(stderr, NULL, _IOLBF, 0);
#ifdef T3D
  free(malloc(14000000));
#endif
  initforpvm(argc, argv);
  /* Read physical parameters */
  Cblacs_pinfo(&mypnum, &nprocs);
  /* initialize BLACS for the parameter communication */
  Cblacs_get(0, 0, &gcontext);
  Cblacs_gridinit(&gcontext, "R", nprocs, 1);
  Cblacs_gridinfo(gcontext, &dummy, &dummy, &mypnum, &dummy);
  if (mypnum == 0) {
    if ((fp = fopen("TRMR2D.dat", "r")) == NULL) {
      fprintf(stderr, "Can't open TRMR2D.dat\n");
      exit(1);
    };
    printf("\n// CTRMR2D TESTER for COMPLEX //\n");
    getparam(fp, &nbre, NULL);
    printf("////////// %d tests \n\n", nbre);
    parameters[0] = nbre;
    Cigebs2d(gcontext, "All", "H", 1, 1, parameters, 1);
  } else {
    Cigebr2d(gcontext, "All", "H", 1, 1, parameters, 1, 0, 0);
    nbre = parameters[0];
  };
  if (mypnum == 0) {
    printf("\n  m   n  m0  n0  sr0 sc0 i0  j0  p0  q0 nbr0 nbc0 \
m1  n1  sr1 sc1 i1  j1  p1  q1 nbr1 nbc1\n\n");
  };
  /****** TEST LOOP *****/
  /* Here we are in grip 1xnprocs */
  nbremax = nbre;
#ifdef DEBUG
  fprintf(stderr, "bonjour,je suis le noeud %d\n", mypnum);
#endif
  while (nbre-- != 0) {	/* Loop on the serie of tests */
    /* All the processors read the parameters so we have to be in a 1xnprocs
     * grid at each iteration */
    /* Read processors grid and matrices parameters */
    if (mypnum == 0) {
      int   u, d;
      getparam(fp,
	       &m, &n,
	       &ma.m, &ma.n, &ma.sprow, &ma.spcol,
	       &ia, &ja, &p0, &q0, &ma.nbrow, &ma.nbcol,
	       &mb.m, &mb.n, &mb.sprow, &mb.spcol,
	       &ib, &jb, &p1, &q1, &mb.nbrow, &mb.nbcol,
	       &u, &d,
	       NULL);
      uplo = u ? "UPPER" : "LOWER";
      diag = d ? "UNIT" : "NONUNIT";
      printf("\t\t************* TEST # %d **********\n",
	     nbremax - nbre);
      printf(" %3d %3d %3d %3d %3d %3d %3d %3d \
%3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d %3d",
	     m, n,
	     ma.m, ma.n, ma.sprow, ma.spcol,
	     ia, ja, p0, q0, ma.nbrow, ma.nbcol,
	     mb.m, mb.n, mb.sprow, mb.spcol,
	     ib, jb, p1, q1, mb.nbrow, mb.nbcol);
      printf(" %s %s", toupper(*uplo) == 'U' ? "up" : "low",
	     toupper(*diag) == 'U' ? "unit" : "nonunit");
      printf("\n");
      if (p0 * q0 > nprocs || p1 * q1 > nprocs) {
	fprintf(stderr, "not enough nodes:%d processors required\n",
		max(p0 * q0, p1 * q1));
	exit(1);
      }
      parameters[0] = p0;
      parameters[1] = q0;
      parameters[2] = ma.nbrow;
      parameters[3] = ma.nbcol;
      parameters[4] = p1;
      parameters[5] = q1;
      parameters[6] = mb.nbrow;
      parameters[7] = mb.nbcol;
      parameters[8] = ma.m;
      parameters[9] = ma.n;
      parameters[10] = ma.sprow;
      parameters[11] = ma.spcol;
      parameters[12] = mb.sprow;
      parameters[13] = mb.spcol;
      parameters[14] = ia;
      parameters[15] = ja;
      parameters[16] = ib;
      parameters[17] = jb;
      parameters[18] = m;
      parameters[19] = n;
      parameters[20] = mb.m;
      parameters[21] = mb.n;
      parameters[22] = *uplo == 'U';
      parameters[23] = *diag == 'U';
      Cigebs2d(gcontext, "All", "H", 1, nbparameter, parameters, 1);
    } else {
      Cigebr2d(gcontext, "All", "H", 1, nbparameter, parameters, 1, 0, 0);
      p0 = parameters[0];
      q0 = parameters[1];
      ma.nbrow = parameters[2];
      ma.nbcol = parameters[3];
      p1 = parameters[4];
      q1 = parameters[5];
      mb.nbrow = parameters[6];
      mb.nbcol = parameters[7];
      ma.m = parameters[8];
      ma.n = parameters[9];
      ma.sprow = parameters[10];
      ma.spcol = parameters[11];
      mb.sprow = parameters[12];
      mb.spcol = parameters[13];
      ia = parameters[14];
      ja = parameters[15];
      ib = parameters[16];
      jb = parameters[17];
      m = parameters[18];
      n = parameters[19];
      mb.m = parameters[20];
      mb.n = parameters[21];
      ma.desctype = BLOCK_CYCLIC_2D;
      mb.desctype = BLOCK_CYCLIC_2D;
      uplo = parameters[22] ? "UPPER" : "LOWER";
      diag = parameters[23] ? "UNIT" : "NONUNIT";
    };
    Cblacs_get(0, 0, &context0);
    Cblacs_gridinit(&context0, "R", p0, q0);
    Cblacs_get(0, 0, &context1);
    Cblacs_gridinit(&context1, "R", p1, q1);
    Cblacs_gridinfo(context0, &dummy, &dummy, &myprow0, &mypcol0);
    if (myprow0 >= p0 || mypcol0 >= q0)
      myprow0 = mypcol0 = -1;
    Cblacs_gridinfo(context1, &dummy, &dummy, &myprow1, &mypcol1);
    if (myprow1 >= p1 || mypcol1 >= q1)
      myprow1 = mypcol1 = -1;
    assert((myprow0 < p0 && mypcol0 < q0) || (myprow0 == -1 && mypcol0 == -1));
    assert((myprow1 < p1 && mypcol1 < q1) || (myprow1 == -1 && mypcol1 == -1));
    ma.ctxt = context0;
    mb.ctxt = context1;
    /* From here, we are not assuming that only the processors working in the
     * redistribution are calling  xxMR2D, but the ones not concerned will do
     * nothing. */
    /* We compute the exact size of the local memory block for the memory
     * allocations */
    if (myprow0 >= 0 && mypcol0 >= 0) {
      blocksize0 = memoryblocksize(&ma);
      ma.lda = localsize(SHIFT(myprow0, ma.sprow, p0), p0, ma.nbrow, ma.m);
      setmemory(&ptrmyblock, blocksize0);
      initblock(ptrmyblock, 1, blocksize0);
      setmemory(&ptrmyblockcopy, blocksize0);
      memcpy((char *) ptrmyblockcopy, (char *) ptrmyblock,
	     blocksize0 * sizeof(complex));
      setmemory(&ptrmyblockvide, blocksize0);
      for (i = 0; i < blocksize0; i++)
	ptrmyblockvide[i].r = -1;
    };	/* if (mypnum < p0 * q0) */
    if (myprow1 >= 0 && mypcol1 >= 0) {
      setmemory(&ptrsavemyblock, memoryblocksize(&mb));
      mb.lda = localsize(SHIFT(myprow1, mb.sprow, p1), p1, mb.nbrow, mb.m);
    };	/* if (mypnum < p1 * q1)  */
    /* Redistribute the matrix from grid 0 to grid 1 (memory location
     * ptrmyblock to ptrsavemyblock) */
    Cpctrmr2d(uplo, diag, m, n,
	      ptrmyblock, ia, ja, &ma,
	      ptrsavemyblock, ib, jb, &mb, gcontext);
    /* Perform the inverse redistribution of the matrix from grid 1 to grid 0
     * (memory location ptrsavemyblock to ptrmyblockvide) */
    Cpctrmr2d(uplo, diag, m, n,
	      ptrsavemyblock, ib, jb, &mb,
	      ptrmyblockvide, ia, ja, &ma, gcontext);
    /* Check the differences */
    nberrors = 0;
    if (myprow0 >= 0 && mypcol0 >= 0) {
      /* only for the processors that do have data at the begining */
      for (i = 0; i < blocksize0; i++) {
	int   li, lj, gi, gj;
	int   in;
	in = 1;
	li = i % ma.lda;
	lj = i / ma.lda;
	gi = (li / ma.nbrow) * p0 * ma.nbrow +
	      SHIFT(myprow0, ma.sprow, p0) * ma.nbrow + li % ma.nbrow;
	gj = (lj / ma.nbcol) * q0 * ma.nbcol +
	      SHIFT(mypcol0, ma.spcol, q0) * ma.nbcol + lj % ma.nbcol;
	assert(gi < ma.m && gj < ma.n);
	gi -= (ia - 1);
	gj -= (ja - 1);
	if (gi < 0 || gj < 0 || gi >= m || gj >= n)
	  in = 0;
	else if (toupper(*uplo) == 'U')
	  in = (gi <= gj + max(0, m - n) - (toupper(*diag) == 'U'));
	else
	  in = (gi >= gj - max(0, n - m) + (toupper(*diag) == 'U'));
	if (!in) {
	  ptrmyblockcopy[i].r = -1;
	}
	if (ptrmyblockvide[i].r != ptrmyblockcopy[i].r) {
	  nberrors++;
	  printf("Proc %d : Error element number %d, value = %f , initvalue =%f \n"
		 ,mypnum, i,
		 ptrmyblockvide[i].r, ptrmyblockcopy[i].r);
	};
      };
      if (nberrors > 0) {
	printf("Processor %d, has tested  %d COMPLEX elements,\
Number of redistribution errors = %d \n",
	       mypnum, blocksize0, nberrors);
      }
    }
    /* Look at the errors on all the processors at this point. */
    Cigsum2d(gcontext, "All", "H", 1, 1, &nberrors, 1, 0, 0);
    if (mypnum == 0)
      if (nberrors)
	printf("  => Total number of redistribution errors = %d \n",
	       nberrors);
      else
	printf("TEST PASSED OK\n");
    /* release memory for the next iteration */
    if (myprow0 >= 0 && mypcol0 >= 0) {
      freememory((char *) ptrmyblock);
      freememory((char *) ptrmyblockvide);
      freememory((char *) ptrmyblockcopy);
    };	/* if (mypnum < p0 * q0) */
    /* release memory for the next iteration */
    if (myprow1 >= 0 && mypcol1 >= 0) {
      freememory((char *) ptrsavemyblock);
    };
    if (myprow0 >= 0)
      Cblacs_gridexit(context0);
    if (myprow1 >= 0)
      Cblacs_gridexit(context1);
  };	/* while nbre != 0 */
  if (mypnum == 0) {
    fclose(fp);
  };
  Cblacs_exit(0);
  return 0;
}/* main */
