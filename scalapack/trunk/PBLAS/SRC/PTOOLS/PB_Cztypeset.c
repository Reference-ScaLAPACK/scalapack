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

PBTYP_T * PB_Cztypeset()
{
/*
*  Purpose
*  =======
*
*  PB_Cztypeset on the first call initializes a static structure contai-
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
   static cmplx16 zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = DCPLX;
   TypeStruct.usiz = sizeof( double  );
   TypeStruct.size = sizeof( cmplx16 );

   zero  [REAL_PART] = ZERO;
   zero  [IMAG_PART] = ZERO;
   one   [REAL_PART] =  ONE;
   one   [IMAG_PART] = ZERO;
   negone[REAL_PART] = -ONE;
   negone[IMAG_PART] = ZERO;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Czgesd2d;
   TypeStruct.Cgerv2d   = Czgerv2d;
   TypeStruct.Cgebs2d   = Czgebs2d;
   TypeStruct.Cgebr2d   = Czgebr2d;
   TypeStruct.Cgsum2d   = Czgsum2d;

   TypeStruct.Fmmadd    = zmmadd_;
   TypeStruct.Fmmcadd   = zmmcadd_;
   TypeStruct.Fmmtadd   = zmmtadd_;
   TypeStruct.Fmmtcadd  = zmmtcadd_;
   TypeStruct.Fmmdda    = zmmdda_;
   TypeStruct.Fmmddac   = zmmddac_;
   TypeStruct.Fmmddat   = zmmddat_;
   TypeStruct.Fmmddact  = zmmddact_;

   TypeStruct.Fcshft    = zcshft_;
   TypeStruct.Frshft    = zrshft_;

   TypeStruct.Fvvdotu   = zvvdotu_;
   TypeStruct.Fvvdotc   = zvvdotc_;

   TypeStruct.Fset      = zset_;

   TypeStruct.Ftzpad    = ztzpad_;
   TypeStruct.Ftzpadcpy = ztzpadcpy_;
   TypeStruct.Ftzscal   = ztzscal_;
   TypeStruct.Fhescal   = zhescal_;
   TypeStruct.Ftzcnjg   = ztzcnjg_;

   TypeStruct.Faxpy     = zaxpy_;
   TypeStruct.Fcopy     = zcopy_;
   TypeStruct.Fswap     = zswap_;

   TypeStruct.Fgemv     = zgemv_;
   TypeStruct.Fsymv     = zsymv_;
   TypeStruct.Fhemv     = zhemv_;
   TypeStruct.Ftrmv     = ztrmv_;
   TypeStruct.Ftrsv     = ztrsv_;
   TypeStruct.Fagemv    = zagemv_;
   TypeStruct.Fasymv    = zasymv_;
   TypeStruct.Fahemv    = zahemv_;
   TypeStruct.Fatrmv    = zatrmv_;

   TypeStruct.Fgerc     = zgerc_;
   TypeStruct.Fgeru     = zgeru_;
   TypeStruct.Fsyr      = zsyr_;
   TypeStruct.Fher      = zher_;
   TypeStruct.Fsyr2     = zsyr2_;
   TypeStruct.Fher2     = zher2_;

   TypeStruct.Fgemm     = zgemm_;
   TypeStruct.Fsymm     = zsymm_;
   TypeStruct.Fhemm     = zhemm_;
   TypeStruct.Fsyrk     = zsyrk_;
   TypeStruct.Fherk     = zherk_;
   TypeStruct.Fsyr2k    = zsyr2k_;
   TypeStruct.Fher2k    = zher2k_;
   TypeStruct.Ftrmm     = ztrmm_;
   TypeStruct.Ftrsm     = ztrsm_;

   return( &TypeStruct );
/*
*  End of PB_Cztypeset
*/
}
