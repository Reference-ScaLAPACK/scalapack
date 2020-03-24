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

PBTYP_T * PB_Cstypeset()
{
/*
*  Purpose
*  =======
*
*  PB_Cstypeset on the first call initializes a static structure contai-
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
   static Int     setup=0;
   static PBTYP_T TypeStruct;
   static float   zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = SREAL;
   TypeStruct.usiz = sizeof( float );
   TypeStruct.size = sizeof( float );

   zero   = ZERO;
   one    =  ONE;
   negone = -ONE;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Csgesd2d;
   TypeStruct.Cgerv2d   = Csgerv2d;
   TypeStruct.Cgebs2d   = Csgebs2d;
   TypeStruct.Cgebr2d   = Csgebr2d;
   TypeStruct.Cgsum2d   = Csgsum2d;

   TypeStruct.Fmmadd    = smmadd_;
   TypeStruct.Fmmcadd   = smmcadd_;
   TypeStruct.Fmmtadd   = smmtadd_;
   TypeStruct.Fmmtcadd  = smmtcadd_;
   TypeStruct.Fmmdda    = smmdda_;
   TypeStruct.Fmmddac   = smmddac_;
   TypeStruct.Fmmddat   = smmddat_;
   TypeStruct.Fmmddact  = smmddact_;

   TypeStruct.Fcshft    = scshft_;
   TypeStruct.Frshft    = srshft_;

   TypeStruct.Fvvdotu   = svvdot_;
   TypeStruct.Fvvdotc   = svvdot_;

   TypeStruct.Fset      = sset_;

   TypeStruct.Ftzpad    = stzpad_;
   TypeStruct.Ftzpadcpy = stzpadcpy_;
   TypeStruct.Ftzscal   = stzscal_;
   TypeStruct.Fhescal   = stzscal_;
   TypeStruct.Ftzcnjg   = stzscal_;

   TypeStruct.Faxpy     = saxpy_;
   TypeStruct.Fcopy     = scopy_;
   TypeStruct.Fswap     = sswap_;

   TypeStruct.Fgemv     = sgemv_;
   TypeStruct.Fsymv     = ssymv_;
   TypeStruct.Fhemv     = ssymv_;
   TypeStruct.Ftrmv     = strmv_;
   TypeStruct.Ftrsv     = strsv_;
   TypeStruct.Fagemv    = sagemv_;
   TypeStruct.Fasymv    = sasymv_;
   TypeStruct.Fahemv    = sasymv_;
   TypeStruct.Fatrmv    = satrmv_;

   TypeStruct.Fgerc     = sger_;
   TypeStruct.Fgeru     = sger_;
   TypeStruct.Fsyr      = ssyr_;
   TypeStruct.Fher      = ssyr_;
   TypeStruct.Fsyr2     = ssyr2_;
   TypeStruct.Fher2     = ssyr2_;

   TypeStruct.Fgemm     = sgemm_;
   TypeStruct.Fsymm     = ssymm_;
   TypeStruct.Fhemm     = ssymm_;
   TypeStruct.Fsyrk     = ssyrk_;
   TypeStruct.Fherk     = ssyrk_;
   TypeStruct.Fsyr2k    = ssyr2k_;
   TypeStruct.Fher2k    = ssyr2k_;
   TypeStruct.Ftrmm     = strmm_;
   TypeStruct.Ftrsm     = strsm_;

   return( &TypeStruct );
/*
*  End of PB_Cstypeset
*/
}
