/*
 *  This file includes the standard C libraries, as well as system dependant
 *  include files.  All BLACS routines include this file.
 */

#ifndef BCONFIG_H
#define BCONFIG_H 1

/*
 * Include files
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <mpi.h>

/*
 * Integer types used by BLACS
 */
#ifndef Int
#define Int int
#endif
#ifndef MpiInt
#define MpiInt int
#endif

/*
 * These macros define the naming strategy needed for a fortran
 * routine to call a C routine, and whether to build so they may be
 * called from C or fortran.  For the fortran call C interface, ADD_ assumes that
 * fortran calls expect C routines to have an underscore postfixed to the name
 * (Suns, and the Intel expect this).  NOCHANGE indicates that fortran expects
 * the name called by fortran to be identical to that compiled by C
 * (AIX does this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this).  The variable FORTRAN_CALL_C is always
 * set to one of these values.  If the BLACS will be called from C, we define
 * INTFACE to be CALL_C, otherwise, it is set to FORTRAN_CALL_C.
 */
#define ADD_     0
#define NOCHANGE 1
#define UPCASE   2
#define FCISF2C  3
#define C_CALL   4

#ifdef UpCase
#define FORTRAN_CALL_C UPCASE
#endif

#ifdef NoChange
#define FORTRAN_CALL_C NOCHANGE
#endif

#ifdef Add_
#define FORTRAN_CALL_C ADD_
#endif

#ifdef FortranIsF2C
#define FORTRAN_CALL_C FCISF2C
#endif

#ifndef FORTRAN_CALL_C
#define FORTRAN_CALL_C ADD_
#endif

#ifdef CallFromC
#define INTFACE C_CALL
#else
#define INTFACE FORTRAN_CALL_C
#endif

/*
 *  Uncomment these macro definitions, and substitute the topology of your
 *  choice to vary the default topology (TOP = ' ') for broadcast and combines.
#define DefBSTop '1'
#define DefCombTop '1'
 */

/*
 * Uncomment this line if your MPI_Send provides a locally-blocking send
 */
//#define SndIsLocBlk

/*
 * Comment out the following line if your MPI does a data copy on every
 * non-contiguous send
 */
#define MpiBuffGood

/*
 * If your MPI cannot form data types of zero length, uncomment the
 * following definition
 */
/* #define ZeroByteTypeBug */

/*
 *  These macros set the timing and debug levels for the BLACS.  The fastest
 *  code is produced by setting both values to 0.  Higher levels provide
 *  more timing/debug information at the cost of performance.  Present levels
 *  of debug are:
 *  0 : No debug information
 *  1 : Mainly parameter checking.
 *
 *  Present levels of timing are:
 *  0 : No timings taken
 */
#ifndef BlacsDebugLvl
#define BlacsDebugLvl 0
#endif
#ifndef BlacsTimingLvl
#define BlacsTimingLvl 0
#endif

#endif
