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

PBTYP_T * PB_Cctypeset(void)
{
/*
*  Purpose
*  =======
*
*  PB_Cctypeset on the first call initializes a static structure contai-
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
   static cmplx   zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = SCPLX;
   TypeStruct.usiz = sizeof( float );
   TypeStruct.size = sizeof( cmplx );

   zero  [REAL_PART] = ZERO;
   zero  [IMAG_PART] = ZERO;
   one   [REAL_PART] =  ONE;
   one   [IMAG_PART] = ZERO;
   negone[REAL_PART] = -ONE;
   negone[IMAG_PART] = ZERO;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Ccgesd2d;
   TypeStruct.Cgerv2d   = Ccgerv2d;
   TypeStruct.Cgebs2d   = Ccgebs2d;
   TypeStruct.Cgebr2d   = Ccgebr2d;
   TypeStruct.Cgsum2d   = Ccgsum2d;

   TypeStruct.Fmmadd    = cmmadd_;
   TypeStruct.Fmmcadd   = cmmcadd_;
   TypeStruct.Fmmtadd   = cmmtadd_;
   TypeStruct.Fmmtcadd  = cmmtcadd_;
   TypeStruct.Fmmdda    = cmmdda_;
   TypeStruct.Fmmddac   = cmmddac_;
   TypeStruct.Fmmddat   = cmmddat_;
   TypeStruct.Fmmddact  = cmmddact_;

   TypeStruct.Fcshft    = ccshft_;
   TypeStruct.Frshft    = crshft_;

   TypeStruct.Fvvdotu   = cvvdotu_;
   TypeStruct.Fvvdotc   = cvvdotc_;

   TypeStruct.Fset      = cset_;

   TypeStruct.Ftzpad    = ctzpad_;
   TypeStruct.Ftzpadcpy = ctzpadcpy_;
   TypeStruct.Ftzscal   = ctzscal_;
   TypeStruct.Fhescal   = chescal_;
   TypeStruct.Ftzcnjg   = ctzcnjg_;

   TypeStruct.Faxpy     = caxpy_;
   TypeStruct.Fcopy     = ccopy_;
   TypeStruct.Fswap     = cswap_;

   TypeStruct.Fgemv     = cgemv_;
   TypeStruct.Fsymv     = csymv_;
   TypeStruct.Fhemv     = chemv_;
   TypeStruct.Ftrmv     = ctrmv_;
   TypeStruct.Ftrsv     = ctrsv_;
   TypeStruct.Fagemv    = cagemv_;
   TypeStruct.Fasymv    = casymv_;
   TypeStruct.Fahemv    = cahemv_;
   TypeStruct.Fatrmv    = catrmv_;

   TypeStruct.Fgerc     = cgerc_;
   TypeStruct.Fgeru     = cgeru_;
   TypeStruct.Fsyr      = csyr_;
   TypeStruct.Fher      = cher_;
   TypeStruct.Fsyr2     = csyr2_;
   TypeStruct.Fher2     = cher2_;

   TypeStruct.Fgemm     = cgemm_;
   TypeStruct.Fsymm     = csymm_;
   TypeStruct.Fhemm     = chemm_;
   TypeStruct.Fsyrk     = csyrk_;
   TypeStruct.Fherk     = cherk_;
   TypeStruct.Fsyr2k    = csyr2k_;
   TypeStruct.Fher2k    = cher2k_;
   TypeStruct.Ftrmm     = ctrmm_;
   TypeStruct.Ftrsm     = ctrsm_;

   return( &TypeStruct );
/*
*  End of PB_Cctypeset
*/
}
