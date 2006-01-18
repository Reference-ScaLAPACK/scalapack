#include "redist.h"
/* $Id: pgemraux.c,v 1.1.1.1 2000/02/15 18:04:10 susan Exp $
 * 
 * some functions used by the pigemr2d routine see file pigemr.c for more
 * documentation.
 * 
 * Created March 1993 by B. Tourancheau (See sccs for modifications). */
#define static2 static
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2d pigemr2do_
#define fortran_mr2dnew pigemr2d_
#elif defined(UpCase)
#define fortran_mr2dnew PIGEMR2D
#define fortran_mr2d PIGEMR2DO
#define icopy_ ICOPY
#define ilacpy_ ILACPY
#else
#define fortran_mr2d pigemr2do
#define fortran_mr2dnew pigemr2d
#define icopy_ icopy
#define ilacpy_ ilacpy
#endif
#define Clacpy Cigelacpy
void  Clacpy();
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
  int   lstart;
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
extern void Cigesd2d();
extern void Cigerv2d();
/* lapack */
void  ilacpy_();
/* aux fonctions */
extern int localindice();
extern void *mr2d_malloc();
extern int ppcm();
extern int localsize();
extern int memoryblocksize();
extern int changeorigin();
extern void paramcheck();
/* tools and others function */
#define scanD0 igescanD0
#define dispmat igedispmat
#define setmemory igesetmemory
#define freememory igefreememory
#define scan_intervals igescan_intervals
extern void scanD0();
extern void dispmat();
extern void setmemory();
extern void freememory();
extern int scan_intervals();
extern void Cpigemr2do();
extern void Cpigemr2d();
/* some defines for Cpigemr2do */
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
void *
mr2d_malloc(n)
  int   n;
{
  void *ptr;
  assert(n > 0);
  ptr = (void *) malloc(n);
  if (ptr == NULL) {
    fprintf(stderr, "xxmr2d:out of memory\n");
    exit(2);
  }
  return ptr;
}
int 
pgcd(a, b)
  int   a, b;
{
  int   aux;
  if (a < b)
    return pgcd(b, a);
  else {
    aux = a % b;
    if (aux == 0)
      return b;
    else
      return pgcd(b, aux);
  }
}
int 
ppcm(a, b)
  int   a, b;
{
  int   pg;
  pg = pgcd(a, b);
  return a * (b / pg);
}
/* localsize:return the number of rows on the local processor given by its
 * row number myprow, of a distributed matrix with m rows distributed of on a
 * grid of processors with p rows with blocksize nbrow : this procedure can
 * also be used to compute the number of cols by replacing rows by cols */
int 
localsize(myprow, p, nbrow, m)
  int   myprow, p, nbrow, m;
{
  int   templateheight, blockheight;
  templateheight = p * nbrow;
  if (m % templateheight != 0) {	/* not an exact boundary */
    if ((m % templateheight) > (nbrow * myprow)) {	/* processor
							 * (myprow,mypcol) has
							 * some elements in that
							 * incomplete template */
      if ((m % templateheight) >= (nbrow * (myprow + 1))) {	/* processor
								 * (myprow,mypcol)'s
								 * part is complete */
	blockheight = (m / templateheight) * nbrow + nbrow;
      } else {	/* processor (myprow,mypcol)'s part is not complete */
	blockheight = (m / templateheight) * nbrow + (m % nbrow);
      };	/* if ((m%templateheight) > (nbrow*(myprow+1))) */
    } else {	/* processor (myprow,mypcol) has no element in that
		 * incomplete template */
      blockheight = (m / templateheight) * nbrow;
    };	/* if ((m%templateheight) > (nbrow*myprow)) */
  } else {	/* exact boundary */
    blockheight = m / p;	/* (m/templateheight) * nbrow */
  };	/* if (m%templateheight !=0) */
  return blockheight;
}
/****************************************************************/
/* Returns the exact memory block size corresponding to the parameters */
int
memoryblocksize(a)
  MDESC *a;
{
  int   myprow, mypcol, p, q;
  /* Compute the (myprow,mypcol) indices of processor mypnum in P0xQ0 We
   * assume the row-major ordering of the BLACS */
  Cblacs_gridinfo(a->ctxt, &p, &q, &myprow, &mypcol);
  myprow = SHIFT(myprow, a->sprow, p);
  mypcol = SHIFT(mypcol, a->spcol, q);
  assert(myprow >= 0 && mypcol >= 0);
  return localsize(myprow, p, a->nbrow, a->m) *
	localsize(mypcol, q, a->nbcol, a->n);
}
void 
checkequal(ctxt, a)
  int   a, ctxt;
{
  int   np, dummy, nbrow, myp, b;
  Cblacs_gridinfo(ctxt, &nbrow, &np, &dummy, &myp);
  assert(nbrow == 1);
  if (np == 1)
    return;
  if (myp == 0) {
    Cigesd2d(ctxt, 1, 1, &a, 1, 0, 1);
    Cigerv2d(ctxt, 1, 1, &b, 1, 0, np - 1);
    assert(a == b);
  } else {
    Cigerv2d(ctxt, 1, 1, &b, 1, 0, myp - 1);
    assert(a == b);
    Cigesd2d(ctxt, 1, 1, &a, 1, 0, (myp + 1) % np);
  }
}
void 
paramcheck(a, i, j, m, n, p, q, gcontext)
  MDESC *a;
  int   i, j, m, n, p, q;
{
  int   p2, q2, myprow, mypcol;
#ifndef NDEBUG
  checkequal(gcontext, p);
  checkequal(gcontext, q);
  checkequal(gcontext, a->sprow);
  checkequal(gcontext, a->spcol);
  checkequal(gcontext, a->m);
  checkequal(gcontext, a->n);
  checkequal(gcontext, i);
  checkequal(gcontext, j);
  checkequal(gcontext, a->nbrow);
  checkequal(gcontext, a->nbcol);
#endif
  Cblacs_gridinfo(a->ctxt, &p2, &q2, &myprow, &mypcol);
  /* compatibility T3D, must check myprow  and mypcol are within bounds */
  if (myprow >= p2 || mypcol >= q2)
    myprow = mypcol = -1;
  if ((myprow >= 0 || mypcol >= 0) && (p2 != p && q2 != q)) {
    fprintf(stderr, "??MR2D:incoherent p,q parameters\n");
    exit(1);
  }
  assert(myprow < p && mypcol < q);
  if (a->sprow < 0 || a->sprow >= p || a->spcol < 0 || a->spcol >= q) {
    fprintf(stderr, "??MR2D:Bad first processor coordinates\n");
    exit(1);
  }
  if (i < 0 || j < 0 || i + m > a->m || j + n > a->n) {
    fprintf(stderr, "??MR2D:Bad submatrix:i=%d,j=%d,\
m=%d,n=%d,M=%d,N=%d\n",
	    i, j, m, n, a->m, a->n);
    exit(1);
  }
  if ((myprow >= 0 || mypcol >= 0) &&
      localsize(SHIFT(myprow, a->sprow, p), p, a->nbrow, a->m) > a->lda) {
    fprintf(stderr, "??MR2D:bad lda arg:row=%d,m=%d,p=%d,\
nbrow=%d,lda=%d,sprow=%d\n",
	    myprow, a->m, p, a->nbrow, a->lda, a->sprow);
    exit(1);
  }
}
/* to change from the submatrix beginning at line i to one beginning at line
 * i' with i'< blocksize return the line number on the local process where
 * the new matrix begin, the new process number, and i' */
int 
changeorigin(myp, sp, p, bs, i, decal, newsp)
  int   myp, sp, p, bs, i;
  int  *decal, *newsp;
{
  int   tempheight, firstblock, firsttemp;
  /* we begin by changing the parameters so that ia < templatewidth,... */
  tempheight = bs * p;
  firsttemp = i / tempheight;
  firstblock = (i / bs) % p;
  *newsp = (sp + firstblock) % p;
  if (myp >= 0)
    *decal = firsttemp * bs + (SHIFT(myp, sp, p) < firstblock ? bs : 0);
  else
    *decal = 0;
  return i % bs;
}
/******************************************************************/
/* Return the indice in local memory of element of indice a in the matrix */
int
localindice(ig, jg, templateheight, templatewidth, a)
  int   templateheight, templatewidth, ig, jg;
  MDESC *a;
/* Return the indice in local memory (scattered distribution) of the element
 * of indice a in global matrix */
{
  int   vtemp, htemp, vsubtemp, hsubtemp, il, jl;
  assert(ig >= 0 && ig < a->m && jg >= 0 && jg < a->n);
  /* coordinates in global matrix with the tests in intersect, ig MUST BE in
   * [0..m] and jg in [0..n] */
  /* coordinates of the template that "owns" the element */
  vtemp = ig / templateheight;
  htemp = jg / templatewidth;
  /* coordinates of the element in the subblock of the (vtemp, htemp)
   * template */
  vsubtemp = ig % a->nbrow;
  hsubtemp = jg % a->nbcol;
  /* coordinates of the element in the local block of the processor */
  il = a->nbrow * vtemp + vsubtemp;
  jl = a->nbcol * htemp + hsubtemp;
  assert(il < a->lda);
#ifndef NDEBUG
  {
    int   pr, pc, p, q, lp, lq;
    Cblacs_gridinfo(a->ctxt, &p, &q, &pr, &pc);
    p = templateheight / a->nbrow;
    q = templatewidth / a->nbcol;
    lp = ig % templateheight / a->nbrow;
    lq = jg % templatewidth / a->nbcol;
    assert(lp == SHIFT(pr, a->sprow, p));
    assert(lq == SHIFT(pc, a->spcol, q));
  }
#endif
  return (jl * a->lda + il);
}
