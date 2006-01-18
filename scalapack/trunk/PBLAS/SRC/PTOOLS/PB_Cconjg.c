/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "../pblas.h"
#include "../PBpblas.h"
#include "../PBtools.h"
#include "../PBblacs.h"
#include "../PBblas.h"

#ifdef __STDC__
void PB_Cconjg( PBTYP_T * TYPE, char * ALPHA, char * CALPHA )
#else
void PB_Cconjg( TYPE, ALPHA, CALPHA )
/*
*  .. Scalar Arguments ..
*/
   char           * ALPHA, * CALPHA;
   PBTYP_T        * TYPE;
#endif
{
/*
*  Purpose
*  =======
*
*  PB_Cconjg conjugates of the scalar alpha.
*
*  Arguments
*  =========
*
*  TYPE    (local input) pointer to a PBTYP_T structure
*          On entry,  TYPE  is a pointer to a structure of type PBTYP_T,
*          that contains type information (See pblas.h).
*
*  ALPHA   (local input) pointer to CHAR
*          On entry, ALPHA specifies the scalar alpha.
*
*  CALPHA  (local output) pointer to CHAR
*          On exit, CALPHA contains the conjugate of the scalar alpha.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/* ..
*  .. Executable Statements ..
*
*/
  switch( TYPE->type )
  {
     case SCPLX:
        ((float*)(CALPHA))[REAL_PART]  =  ((float*)(ALPHA))[REAL_PART];
        ((float*)(CALPHA))[IMAG_PART]  = -((float*)(ALPHA))[IMAG_PART];
        break;
     case DCPLX:
        ((double*)(CALPHA))[REAL_PART] =  ((double*)(ALPHA))[REAL_PART];
        ((double*)(CALPHA))[IMAG_PART] = -((double*)(ALPHA))[IMAG_PART];
        break;
     case SREAL:
        ((float*)(CALPHA))[REAL_PART]  =  ((float*)(ALPHA))[REAL_PART];
        break;
     case DREAL:
        ((double*)(CALPHA))[REAL_PART] =  ((double*)(ALPHA))[REAL_PART];
        break;
     case INT:
        *((int*)(CALPHA))              = *((int*)(ALPHA));
        break;
     default: ;
  }
/*
*  End of PB_Cconjg
*/
}
