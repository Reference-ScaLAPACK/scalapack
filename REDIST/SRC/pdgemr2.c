#include "redist.h"
#include <stddef.h>
/* $Id: pdgemr2.c,v 1.1.1.1 2000/02/15 18:04:09 susan Exp $
 * 
 * some functions used by the pdgemr2d routine see file pdgemr.c for more
 * documentation.
 * 
 * Created March 1993 by B. Tourancheau (See sccs for modifications). */
#define static2 static
#if defined(Add_) || defined(f77IsF2C)
#define fortran_mr2d pdgemr2do_
#define fortran_mr2dnew pdgemr2d_
#elif defined(UpCase)
#define fortran_mr2dnew PDGEMR2D
#define fortran_mr2d PDGEMR2DO
#define dcopy_ DCOPY
#define dlacpy_ DLACPY
#else
#define fortran_mr2d pdgemr2do
#define fortran_mr2dnew pdgemr2d
#define dcopy_ dcopy
#define dlacpy_ dlacpy
#endif
#define Clacpy Cdgelacpy
void  Clacpy( Int m, Int n, double *a, Int lda, double *b, Int ldb );
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
  Int   lstart;
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
extern void Cdgesd2d( Int context, Int m, Int n, double* A, Int lda, Int rdest, Int cdest );
extern void Cdgerv2d( Int context, Int m, Int n, double* A, Int lda, Int rsrc, Int csrc );
/* lapack */
void  dlacpy_();
/* aux fonctions */
extern Int localindice( Int ig, Int jg, Int templateheight, Int templatewidth, MDESC *a );
extern void *mr2d_malloc( size_t n );
extern Int ppcm( Int a, Int b );
extern Int localsize( Int myprow, Int p, Int nbrow, Int m );
extern Int memoryblocksize( MDESC *a );
extern Int changeorigin( Int myp, Int sp, Int p, Int bs, Int i, Int *decal, Int *newsp );
extern void paramcheck( MDESC *a, Int i, Int j, Int m, Int n, Int p, Int q, Int gcontext );
/* tools and others function */
#define scanD0 dgescanD0
#define dispmat dgedispmat
#define setmemory dgesetmemory
#define freememory dgefreememory
#define scan_intervals dgescan_intervals
extern void scanD0();
extern void dispmat();
extern void setmemory( double** ptr, Int size );
extern void freememory( double* ptr );
extern Int scan_intervals( char type, Int ja, Int jb, Int n, MDESC *ma, MDESC *mb, Int q0, Int q1, Int col0, Int col1, IDESC *result );
extern void Cpdgemr2do();
extern void Cpdgemr2d();
/* some defines for Cpdgemr2do */
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
setmemory(double **adpointer, Int blocksize)
{
  assert(blocksize >= 0);
  if (blocksize == 0) {
    *adpointer = NULL;
    return;
  }
  *adpointer = (double *) mr2d_malloc(
				      (size_t)blocksize * sizeof(double));
}
/******************************************************************/
/* Free the memory space after the malloc */
void
freememory(double *ptrtobefreed)
{
  if (ptrtobefreed == NULL)
    return;
  free((char *) ptrtobefreed);
}
/* extern functions for intersect() extern dcopy_(); */
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
   * beginning of the Cpdgemr2d */
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
    result[offset].lstart = l + start - j0;
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
