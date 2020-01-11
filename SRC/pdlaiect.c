

/* ---------------------------------------------------------------------
*
*  -- ScaLAPACK routine (version 1.5) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*  ---------------------------------------------------------------------
*/

/*
 * Include Files
 */
#include "pxsyevx.h"
#include "pblas.h"
#include <stdio.h>
#include <math.h>
#define  proto(x)	()


void pdlasnbt_( Int *ieflag )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdalsnbt finds the position of the signbit of a double
*  double precision floating point number. This routine assumes IEEE
*  arithmetic, and hence, tests only the 32nd and 64th bits as
*  possibilities for the sign bit.
*
*  Note : For this release, we assume that sizeof(int) is 4 bytes.
*
*  Note : If a compile time flag (NO_IEEE) indicates that the
*  machine does not have IEEE arithmetic, IEFLAG = 0 is returned.
*
*  Arguments
*  =========
*
*  IEFLAG   (output) INTEGER
*           This indicates the position of the signbit of any double
*           precision floating point number.
*           IEFLAG = 0 if the compile time flag, NO_IEEE, indicates
*           that the machine does not have IEEE  arithmetic, or if
*           sizeof(int) is different from 4 bytes.
*           IEFLAG = 1 indicates that the sign bit is the 32nd
*           bit ( Big Endian ).
*           IEFLAG = 2 indicates that the sign bit is the 64th
*           bit ( Little Endian ).
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double x;
   Int         negone=-1, errornum;
   unsigned Int *ix; 
/* ..
*  .. Executable Statements ..
*/

#ifdef NO_IEEE
   *ieflag = 0;
#else
   if(sizeof(Int) != 4){
      *ieflag = 0;
      return;
   }
   x = (double) -1.0;
   ix = (unsigned Int *) &x;
   if(( *ix == 0xbff00000) && ( *(ix+1) == 0x0) ) 
   {
      *ieflag = 1;
   } else if(( *(ix+1) == 0xbff00000) && ( *ix == 0x0) ) {
      *ieflag = 2;
   } else {
      *ieflag = 0; 
   }
#endif
}

void pdlaiectb_( double *sigma, Int *n, double *d, Int *count )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlaiectb computes the number of negative eigenvalues of (A- SIGMA I).
*  This implementation of the Sturm Sequence loop exploits IEEE Arithmetic
*  and has no conditionals in the innermost loop. To extract the signbit,
*  this routine assumes that the double precision word is stored in
*  "Big Endian" word order, i.e, the signbit is assumed to be bit 32.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  SIGMA    (input) DOUBLE PRECISION
*           The shift. pdlaiectb finds the number of eigenvalues
*           less than equal to SIGMA.
*
*  N        (input) INTEGER
*           The order of the tridiagonal matrix T. N >= 1.
*
*  D        (input) DOUBLE PRECISION array, dimension (2*N - 1)
*           Contains the diagonals and the squares of the off-diagonal
*           elements of the tridiagonal matrix T. These elements are
*           assumed to be interleaved in memory for better cache
*           performance. The diagonal entries of T are in the entries
*           D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*           entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*           matrix must be scaled so that its largest entry is no greater
*           than overflow**(1/2) * underflow**(1/4) in absolute value,
*           and for greatest accuracy, it should not be much smaller
*           than that.
*
*  COUNT    (output) INTEGER
*           The count of the number of eigenvalues of T less than or
*           equal to SIGMA.
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double      lsigma, tmp, *pd, *pe2;
   Int         i;
/* ..
*  .. Executable Statements ..
*/

   lsigma = *sigma;
   pd = d; pe2 = d+1;
   tmp = *pd - lsigma; pd += 2;
   *count = (*((Int *)&tmp) >> 31) & 1;
   for(i = 1;i < *n;i++){
      tmp = *pd - *pe2/tmp - lsigma;
      pd += 2; pe2 += 2;
      *count += ((*((Int *)&tmp)) >> 31) & 1;
   }
}

void pdlaiectl_( double *sigma, Int *n, double *d, Int *count )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlaiectl computes the number of negative eigenvalues of (A- SIGMA I).
*  This implementation of the Sturm Sequence loop exploits IEEE Arithmetic
*  and has no conditionals in the innermost loop. To extract the signbit,
*  this routine assumes that the double precision word is stored in
*  "Little Endian" word order, i.e, the signbit is assumed to be bit 64.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  SIGMA    (input) DOUBLE PRECISION
*           The shift. pdlaiectl finds the number of eigenvalues
*           less than equal to SIGMA.
*
*  N        (input) INTEGER
*           The order of the tridiagonal matrix T. N >= 1.
*
*  D        (input) DOUBLE PRECISION array, dimension (2*N - 1)
*           Contains the diagonals and the squares of the off-diagonal
*           elements of the tridiagonal matrix T. These elements are
*           assumed to be interleaved in memory for better cache
*           performance. The diagonal entries of T are in the entries
*           D(1),D(3),...,D(2*N-1), while the squares of the off-diagonal
*           entries are D(2),D(4),...,D(2*N-2). To avoid overflow, the
*           matrix must be scaled so that its largest entry is no greater
*           than overflow**(1/2) * underflow**(1/4) in absolute value,
*           and for greatest accuracy, it should not be much smaller
*           than that.
*
*  COUNT    (output) INTEGER
*           The count of the number of eigenvalues of T less than or
*           equal to SIGMA.
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double      lsigma, tmp, *pd, *pe2;
   Int         i;
/* ..
*  .. Executable Statements ..
*/

   lsigma = *sigma;
   pd = d; pe2 = d+1;
   tmp = *pd - lsigma; pd += 2;
   *count = (*(((Int *)&tmp)+1) >> 31) & 1;
   for(i = 1;i < *n;i++){
      tmp = *pd - *pe2/tmp - lsigma;
      pd += 2; pe2 += 2;
      *count += (*(((Int *)&tmp)+1) >> 31) & 1;
   }
}

void pdlachkieee_( Int *isieee, double *rmax, double *rmin )
{
/* 
*
*  Purpose
*  ======= 
*
*  pdlachkieee performs a simple check to make sure that the features
*  of the IEEE standard that we rely on are implemented.  In some
*  implementations, pdlachkieee may not return.
*
*  Note that all arguments are call-by-reference so that this routine
*  can be directly called from Fortran code.
*
*  This is a ScaLAPACK internal subroutine and arguments are not
*  checked for unreasonable values.
*
*  Arguments
*  =========
*
*  ISIEEE   (local output) INTEGER
*           On exit, ISIEEE = 1 implies that all the features of the
*           IEEE standard that we rely on are implemented.
*           On exit, ISIEEE = 0 implies that some the features of the
*           IEEE standard that we rely on are missing.
*
*  RMAX     (local input) DOUBLE PRECISION
*           The overflow threshold ( = DLAMCH('O') ).
*
*  RMIN     (local input) DOUBLE PRECISION
*           The underflow threshold ( = DLAMCH('U') ).
*
*  =====================================================================
*
*  .. Local Scalars ..
*/
   double x, pinf, pzero, ninf, nzero;
   Int         ieflag, *ix, sbit1, sbit2, negone=-1, errornum;
/* ..
*  .. Executable Statements ..
*/

   pdlasnbt_( &ieflag );

   pinf = *rmax / *rmin;
   pzero = 1.0 / pinf;
   pinf = 1.0 / pzero;

   if( pzero != 0.0 ){
      printf("pzero = %g should be zero\n",pzero);
      *isieee = 0; 
      return ;
   }
   if( ieflag == 1 ){
      sbit1 = (*((Int *)&pzero) >> 31) & 1;
      sbit2 = (*((Int *)&pinf) >> 31) & 1;
   }else if(ieflag == 2){
      sbit1 = (*(((Int *)&pzero)+1) >> 31) & 1;
      sbit2 = (*(((Int *)&pinf)+1) >> 31) & 1;
   }
   if( sbit1 == 1 ){
      printf("Sign of positive infinity is incorrect\n");
      *isieee = 0;
   }
   if( sbit2 == 1 ){
      printf("Sign of positive zero is incorrect\n");
      *isieee = 0;
   }

   ninf = -pinf;
   nzero = 1.0 / ninf;
   ninf = 1.0 / nzero;

   if( nzero != 0.0 ){
      printf("nzero = %g should be zero\n",nzero);
      *isieee = 0;
   }
   if( ieflag == 1 ){
      sbit1 = (*((Int *)&nzero) >> 31) & 1;
      sbit2 = (*((Int *)&ninf) >> 31) & 1;
   }else if(ieflag == 2){
      sbit1 = (*(((Int *)&nzero)+1) >> 31) & 1;
      sbit2 = (*(((Int *)&ninf)+1) >> 31) & 1;
   }
   if( sbit1 == 0 ){
      printf("Sign of negative infinity is incorrect\n");
      *isieee = 0;
   }
   if( sbit2 == 0 ){
      printf("Sign of negative zero is incorrect\n");
      *isieee = 0;
   }
}
