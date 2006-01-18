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

PBTYP_T * PB_Citypeset()
{
/*
*  Purpose
*  =======
*
*  PB_Citypeset on the first call initializes a static structure contai-
*  ning typed information and returns a pointer to it.  The  other calls
*  to this routine just returns this pointer.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   static int     setup=0;
   static PBTYP_T TypeStruct;
   static int     zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = INT;
   TypeStruct.usiz = sizeof( int );
   TypeStruct.size = sizeof( int );
   zero   =  0;
   one    =  1;
   negone = -1;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cigesd2d;
   TypeStruct.Cgerv2d   = Cigerv2d;
   TypeStruct.Cgebs2d   = Cigebs2d;
   TypeStruct.Cgebr2d   = Cigebr2d;
   TypeStruct.Cgsum2d   = Cigsum2d;

   TypeStruct.Fmmadd    = immadd_;
   TypeStruct.Fmmcadd   = immadd_;
   TypeStruct.Fmmtadd   = immtadd_;
   TypeStruct.Fmmtcadd  = immtadd_;
   TypeStruct.Fmmdda    = immdda_;
   TypeStruct.Fmmddac   = immdda_;
   TypeStruct.Fmmddat   = immddat_;
   TypeStruct.Fmmddact  = immddat_;

   TypeStruct.Fcshft    = NULL;
   TypeStruct.Frshft    = NULL;

   TypeStruct.Fvvdotu   = NULL;
   TypeStruct.Fvvdotc   = NULL;

   TypeStruct.Fset      = NULL;

   TypeStruct.Ftzpad    = NULL;
   TypeStruct.Ftzpadcpy = NULL;

   TypeStruct.Ftzscal   = NULL;
   TypeStruct.Fhescal   = NULL;
   TypeStruct.Ftzcnjg   = NULL;

   TypeStruct.Faxpy     = NULL;
   TypeStruct.Fcopy     = NULL;
   TypeStruct.Fswap     = NULL;

   TypeStruct.Fgemv     = NULL;
   TypeStruct.Fsymv     = NULL;
   TypeStruct.Fhemv     = NULL;
   TypeStruct.Ftrmv     = NULL;
   TypeStruct.Ftrsv     = NULL;
   TypeStruct.Fagemv    = NULL;
   TypeStruct.Fasymv    = NULL;
   TypeStruct.Fahemv    = NULL;
   TypeStruct.Fatrmv    = NULL;

   TypeStruct.Fgerc     = NULL;
   TypeStruct.Fgeru     = NULL;
   TypeStruct.Fsyr      = NULL;
   TypeStruct.Fher      = NULL;
   TypeStruct.Fsyr2     = NULL;
   TypeStruct.Fher2     = NULL;

   TypeStruct.Fgemm     = NULL;
   TypeStruct.Fsymm     = NULL;
   TypeStruct.Fhemm     = NULL;
   TypeStruct.Fsyrk     = NULL;
   TypeStruct.Fherk     = NULL;
   TypeStruct.Fsyr2k    = NULL;
   TypeStruct.Fher2k    = NULL;
   TypeStruct.Ftrmm     = NULL;
   TypeStruct.Ftrsm     = NULL;

   return( &TypeStruct );
/*
*  End of PB_Citypeset
*/
}
