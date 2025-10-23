#include "redist.h"
#include <stddef.h>
/* $Id: pztrmr2.c,v 1.1.1.1 2000/02/15 18:04:10 susan Exp $
 * 
 * some functions used by the pztrmr2d routine see file pztrmr.c for more
 * documentation.
 * 
 * Created March 1993 by B. Tourancheau (See sccs for modifications). */
#define static2 static
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2d pztrmr2do_
#define fortran_mr2dnew pztrmr2d_
#elif defined(UpCase)
#define fortran_mr2dnew PZTRMR2D
#define fortran_mr2d PZTRMR2DO
#define zcopy_ ZCOPY
#define zlacpy_ ZLACPY
#else
#define fortran_mr2d pztrmr2do
#define fortran_mr2dnew pztrmr2d
#define zcopy_ zcopy
#define zlacpy_ zlacpy
#endif
typedef struct {
  double r, i;
}     dcomplex;
#define Clacpy Cztrlacpy
void  Clacpy( Int m, Int n, dcomplex *a, Int lda, dcomplex *b, Int ldb );
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
extern void Cblacs_pcoord( Int context, Int pnum, Int* prow, Int* pcol );
extern Int Cblacs_pnum( Int context, Int prow, Int pcol );
extern void Csetpvmtids();
extern void Cblacs_get( Int context, Int what, Int* val );
extern void Cblacs_pinfo( Int* mypnum, Int* nprocs );
extern void Cblacs_gridinfo( Int context, Int* nprow, Int* npcol, Int* myrow, Int* mycol );
extern void Cblacs_gridinit( Int* context, char* order, Int nprow, Int npcol );
extern void Cblacs_exit( Int continue_blacs );
extern void Cblacs_gridexit( Int context );
extern void Cblacs_setup( Int* mypnum, Int* nprocs );
extern void Cigebs2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda );
extern void Cigebr2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int rsrc, Int csrc );
extern void Cigesd2d( Int context, Int m, Int n, Int* A, Int lda, Int rdest, Int cdest );
extern void Cigerv2d( Int context, Int m, Int n, Int* A, Int lda, Int rsrc, Int csrc );
extern void Cigsum2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int rdest, Int cdest );
extern void Cigamn2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int* RA, Int* CA, Int rcflag, Int rdest, Int cdest );
extern void Cigamx2d( Int context, char* scope, char* top, Int m, Int n, Int* A, Int lda, Int* RA, Int* CA, Int rcflag, Int rdest, Int cdest );
extern void Czgesd2d( Int context, Int m, Int n, dcomplex* A, Int lda, Int rdest, Int cdest );
extern void Czgerv2d( Int context, Int m, Int n, dcomplex* A, Int lda, Int rsrc, Int csrc );
/* lapack */
void  zlacpy_();
/* aux fonctions */
extern Int localindice( Int ig, Int jg, Int templateheight, Int templatewidth, MDESC *a );
extern void *mr2d_malloc( size_t n );
extern Int ppcm( Int a, Int b );
extern Int localsize( Int myprow, Int p, Int nbrow, Int m );
extern Int memoryblocksize( MDESC *a );
extern Int changeorigin( Int myp, Int sp, Int p, Int bs, Int i, Int *decal, Int *newsp );
extern void paramcheck( MDESC *a, Int i, Int j, Int m, Int n, Int p, Int q, Int gcontext );
/* tools and others function */
#define scanD0 ztrscanD0
#define dispmat ztrdispmat
#define setmemory ztrsetmemory
#define freememory ztrfreememory
#define scan_intervals ztrscan_intervals
extern void scanD0( char* uplo, char* diag, Int action, dcomplex* ptrbuff, Int* ptrsizebuff, Int m, Int n, MDESC* ma, Int ia, Int ja, Int p0, Int q0, MDESC* mb, Int ib, Int jb, Int p1, Int q1, IDESC* v_inter, Int vinter_nb, IDESC* h_inter, Int hinter_nb, dcomplex* ptrblock );
extern void dispmat();
extern void setmemory( dcomplex** ptr, Int size );
extern void freememory( dcomplex* ptr );
extern Int scan_intervals( char type, Int ja, Int jb, Int n, MDESC *ma, MDESC *mb, Int q0, Int q1, Int col0, Int col1, IDESC *result );
extern void Cpztrmr2do();
extern void Cpztrmr2d();
/* some defines for Cpztrmr2do */
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
#include <assert.h>
#include <ctype.h>
/* Created March 1993 by B. Tourancheau (See sccs for modifications). */
/************************************************************************/
/* Set the memory space with the malloc function */
void
setmemory(dcomplex **adpointer, Int blocksize)
{
  assert(blocksize >= 0);
  if (blocksize == 0) {
    *adpointer = NULL;
    return;
  }
  *adpointer = (dcomplex *) mr2d_malloc(
					(size_t) blocksize * sizeof(dcomplex));
}
/******************************************************************/
/* Free the memory space after the malloc */
void
freememory(dcomplex *ptrtobefreed)
{
  if (ptrtobefreed == NULL)
    return;
  free((char *) ptrtobefreed);
}
/* extern functions for intersect() extern zcopy_(); */
/**************************************************************/
/* return the number of elements int the column after i and the distance of
 * the first one from i, i,j can be negative out of borns, the number of
 * elements returned can be negative (means 0) */
static2 Int
insidemat(uplo, diag, i, j, m, n, offset)
  Int   m, n, i, j;	/* coordonnees de depart, taille de la sous-matrice */
  char *uplo, *diag;
  Int  *offset;
{
  /* tests outside mxn */
  assert(j >= 0 && j < n);
  assert(i >= 0);
  if (toupper(*uplo) == 'U') {
    Int   nbline;	/* number of lines in the j_th column */
    Int   virtualnbline;	/* number of line if we were not limited by m */
    *offset = 0;
    virtualnbline = max(m - n, 0) + j + (toupper(*diag) == 'N');
    nbline = min(virtualnbline, m);
    return nbline - i;
  } else {
    Int   firstline;	/* first line in the j_th column */
    Int   diagcol;	/* column where the diag begin */
    Int   virtualline;	/* virtual first line if the matrix was extended with
			 * negative indices */
    Int   off;
    diagcol = max(n - m, 0);
    virtualline = j - diagcol + (toupper(*diag) == 'U');
    firstline = max(0, virtualline);
    off = max(firstline - i, 0);
    *offset = off;
    i += off;
    return m - i;
  }
}/* insidemat() */
/********************************************************************/
/* Execute an action on the local memories when an intersection occurs (the
 * action can be the filling of the memory buffer, the count of the memory
 * buffer size or the setting of the memory with the element received) */
static2 void
intersect(uplo, diag,
	  j, start, end,
	  action,
	  ptrsizebuff, pptrbuff, ptrblock,
	  m, n,
	  ma, ia, ja, templateheight0, templatewidth0,
	  mb, ib, jb, templateheight1, templatewidth1)
  Int   action, *ptrsizebuff;
  Int   j, start, end;
  dcomplex **pptrbuff, *ptrblock;
  Int   templateheight0, templatewidth0;
  Int   templateheight1, templatewidth1;
  MDESC *ma, *mb;
  Int   ia, ja, ib, jb, m, n;
  char *uplo, *diag;
/* Execute the action on the local memory for the current interval and
 * increment pptrbuff and ptrsizebuff of the intervalsize */
/* Notice that if the interval is contigous in the virtual matrice, it is
 * also contigous in the real one ! */
{
  /* int       un = 1; only when we use dcopy instead of memcpy */
  dcomplex *ptrstart;
  Int   offset, nbline;
  Int   intervalsize;
  assert(start < end);
  assert(j >= 0 && j < n);
  nbline =
	insidemat(uplo, diag, start, j, m, n, &offset);
  if (nbline <= 0)
    return;
  start += offset;
  if (start >= end)
    return;
  intervalsize = min(end - start, nbline);
  (*ptrsizebuff) += intervalsize;
  switch (action) {
  case SENDBUFF:	/* fill buff with local elements to be sent */
    ptrstart = ptrblock + localindice(start + ia, j + ja,
				      templateheight0, templatewidth0, ma);
    memcpy((char *) (*pptrbuff), (char *) ptrstart,
	   intervalsize * sizeof(dcomplex));
    /* zcopy_(&intervalsize, (char *) (ptrstart), &un, (char *) (*pptrbuff),
     * &un); */
    (*pptrbuff) += intervalsize;
    break;
  case RECVBUFF:	/* fill local memory with the values received */
    ptrstart = ptrblock + localindice(start + ib, j + jb,
				      templateheight1, templatewidth1, mb);
    memcpy((char *) ptrstart, (char *) (*pptrbuff),
	   intervalsize * sizeof(dcomplex));
    /* zcopy_(&intervalsize, (char *) (*pptrbuff), &un, (char *) (ptrstart),
     * &un); */
    (*pptrbuff) += intervalsize;
    break;
  case SIZEBUFF:	/* computation of sizebuff */
    break;
  default:
    printf("action is  %d outside the scope of the case [0..2] !! \n ", action);
    exit(0);
    break;
  }	/* switch (action) */
}/* intersect() */
/* scan_intervals: scans two distributions in one dimension, and compute the
 * intersections on the local processor. result must be long enough to
 * contains the result that are stocked in IDESC structure, the function
 * returns the number of intersections found */
Int 
scan_intervals(type, ja, jb, n, ma, mb, q0, q1, col0, col1,
	       result)
  char  type;
  Int   ja, jb, n, q0, q1, col0, col1;
  MDESC *ma, *mb;
  IDESC *result;
{
  Int   offset, j0, j1, templatewidth0, templatewidth1, nbcol0, nbcol1;
  Int   l;	/* local indice on the beginning of the interval */
  assert(type == 'c' || type == 'r');
  nbcol0 = (type == 'c' ? ma->nbcol : ma->nbrow);
  nbcol1 = (type == 'c' ? mb->nbcol : mb->nbrow);
  templatewidth0 = q0 * nbcol0;
  templatewidth1 = q1 * nbcol1;
  {
    Int   sp0 = (type == 'c' ? ma->spcol : ma->sprow);
    Int   sp1 = (type == 'c' ? mb->spcol : mb->sprow);
    j0 = SHIFT(col0, sp0, q0) * nbcol0 - ja;
    j1 = SHIFT(col1, sp1, q1) * nbcol1 - jb;
  }
  offset = 0;
  l = 0;
  /* a small check to verify that the submatrix begin inside the first block
   * of the original matrix, this done by a sort of coordinate change at the
   * beginning of the Cpztrmr2d */
  assert(j0 + nbcol0 > 0);
  assert(j1 + nbcol1 > 0);
  while ((j0 < n) && (j1 < n)) {
    Int   end0, end1;
    Int   start, end;
    end0 = j0 + nbcol0;
    end1 = j1 + nbcol1;
    if (end0 <= j1) {
      j0 += templatewidth0;
      l += nbcol0;
      continue;
    }
    if (end1 <= j0) {
      j1 += templatewidth1;
      continue;
    }
    /* compute the raw intersection */
    start = max(j0, j1);
    start = max(start, 0);
    /* the start is correct now, update the corresponding fields */
    result[offset].gstart = start;
    end = min(end0, end1);
    if (end0 == end) {
      j0 += templatewidth0;
      l += nbcol0;
    }
    if (end1 == end)
      j1 += templatewidth1;
    /* throw the limit if they go out of the matrix */
    end = min(end, n);
    assert(end > start);
    /* it is a bit tricky to see why the length is always positive after all
     * this min and max, first we have the property that every interval
     * considered is at least partly into the submatrix, second we arrive
     * here only if the raw intersection is non-void, if we remove a limit
     * that means the corresponding frontier is in both intervals which
     * proove the final interval is non-void, clear ?? */
    result[offset].len = end - start;
    offset += 1;
  }	/* while */
  return offset;
}
/*********************************************************************/
/* Do the scanning of intervals and the requested action */
void
scanD0(uplo, diag, action, ptrbuff, ptrsizebuff,
       m, n,
       ma, ia, ja, p0, q0,
       mb, ib, jb, p1, q1,
       v_inter, vinter_nb,
       h_inter, hinter_nb,
       ptrblock)
  Int   action,	/* # of the action done on the intersected intervals  */
       *ptrsizebuff;	/* size of the communication ptrbuffer (chosen to be
			 * an output parameter in every cases) */
  dcomplex *ptrbuff	/* address of the communication ptrbuffer (a
			 * suffisant memory space is supposed to be allocated
      before the call) */ , *ptrblock;
  Int   p0, q0, p1, q1;
  IDESC *v_inter, *h_inter;
  Int   vinter_nb, hinter_nb;
  Int   m, n;
  Int   ia, ja, ib, jb;
  MDESC *ma, *mb;
  char *uplo, *diag;
{/* Rmk: the a+au type addresses are strict bounds as a+au does not belong to
  * the [a..a+au-1] interval of length au */
  Int   templateheight1, templatewidth1;
  Int   templateheight0, templatewidth0;
  Int   h, v;	/* for scanning the intervals */
  /* initializations */
  templateheight1 = p1 * mb->nbrow;
  templateheight0 = p0 * ma->nbrow;
  templatewidth1 = q1 * mb->nbcol;
  templatewidth0 = q0 * ma->nbcol;
  /* we now will deal will logical grids, that's to say we change our
   * numbering of processors so that (0,0) begin on logical processor (0,0) */
  /* in case we will not enter the while loop */
  (*ptrsizebuff) = 0;
  for (h = 0; h < hinter_nb; h++)
    for (v = 0; v < vinter_nb; v++) {
      Int   j;
      for (j = 0; j < h_inter[h].len; j++)
	intersect(uplo, diag, j + h_inter[h].gstart,
		  v_inter[v].gstart, v_inter[v].gstart + v_inter[v].len,
		  action, ptrsizebuff, &ptrbuff, ptrblock, m, n,
		  ma, ia, ja, templateheight0, templatewidth0,
		  mb, ib, jb, templateheight1, templatewidth1);
    }
}/* scanD0() */
