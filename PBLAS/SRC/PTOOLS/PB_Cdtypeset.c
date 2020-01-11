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

PBTYP_T * PB_Cdtypeset()
{
/*
*  Purpose
*  =======
*
*  PB_Cdtypeset on the first call initializes a static structure contai-
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
   static double  zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = DREAL;
   TypeStruct.usiz = sizeof( double );
   TypeStruct.size = sizeof( double );

   zero   = ZERO;
   one    =  ONE;
   negone = -ONE;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cdgesd2d;
   TypeStruct.Cgerv2d   = Cdgerv2d;
   TypeStruct.Cgebs2d   = Cdgebs2d;
   TypeStruct.Cgebr2d   = Cdgebr2d;
   TypeStruct.Cgsum2d   = Cdgsum2d;

   TypeStruct.Fmmadd    = dmmadd_;
   TypeStruct.Fmmcadd   = dmmcadd_;
   TypeStruct.Fmmtadd   = dmmtadd_;
   TypeStruct.Fmmtcadd  = dmmtcadd_;
   TypeStruct.Fmmdda    = dmmdda_;
   TypeStruct.Fmmddac   = dmmddac_;
   TypeStruct.Fmmddat   = dmmddat_;
   TypeStruct.Fmmddact  = dmmddact_;

   TypeStruct.Fcshft    = dcshft_;
   TypeStruct.Frshft    = drshft_;

   TypeStruct.Fvvdotu   = dvvdot_;
   TypeStruct.Fvvdotc   = dvvdot_;

   TypeStruct.Fset      = dset_;

   TypeStruct.Ftzpad    = dtzpad_;
   TypeStruct.Ftzpadcpy = dtzpadcpy_;
   TypeStruct.Ftzscal   = dtzscal_;
   TypeStruct.Fhescal   = dtzscal_;
   TypeStruct.Ftzcnjg   = dtzscal_;

   TypeStruct.Faxpy     = daxpy_;
   TypeStruct.Fcopy     = dcopy_;
   TypeStruct.Fswap     = dswap_;

   TypeStruct.Fgemv     = dgemv_;
   TypeStruct.Fsymv     = dsymv_;
   TypeStruct.Fhemv     = dsymv_;
   TypeStruct.Ftrmv     = dtrmv_;
   TypeStruct.Ftrsv     = dtrsv_;
   TypeStruct.Fagemv    = dagemv_;
   TypeStruct.Fasymv    = dasymv_;
   TypeStruct.Fahemv    = dasymv_;
   TypeStruct.Fatrmv    = datrmv_;

   TypeStruct.Fgerc     = dger_;
   TypeStruct.Fgeru     = dger_;
   TypeStruct.Fsyr      = dsyr_;
   TypeStruct.Fher      = dsyr_;
   TypeStruct.Fsyr2     = dsyr2_;
   TypeStruct.Fher2     = dsyr2_;

   TypeStruct.Fgemm     = dgemm_;
   TypeStruct.Fsymm     = dsymm_;
   TypeStruct.Fhemm     = dsymm_;
   TypeStruct.Fsyrk     = dsyrk_;
   TypeStruct.Fherk     = dsyrk_;
   TypeStruct.Fsyr2k    = dsyr2k_;
   TypeStruct.Fher2k    = dsyr2k_;
   TypeStruct.Ftrmm     = dtrmm_;
   TypeStruct.Ftrsm     = dtrsm_;

   return( &TypeStruct );
/*
*  End of PB_Cdtypeset
*/
}
