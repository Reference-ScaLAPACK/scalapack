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

/* If somebody is building with BigMPI
 * give it priority compared to other
 * solutions since otherwise they would
 * not have bothered to configure BLACS
 * with BigMPI
 */
#ifdef BIGMPI
#define MpiInt MPI_Count
inline int _MPI_Isend(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, MPI_Request *request) {
	return MPIX_Isend_x(buf, count, datatype, dest, tag, comm, request);
}
inline int _MPI_Irecv(void *buf, MPI_Count count, MPI_Datatype datatype,
		int source, int tag, MPI_Comm comm, MPI_Request *request) {
	return MPIX_Irecv_x (buf, count, datatype, source, tag, comm, request);
}
inline int _MPI_Op_create(MPI_User_function_c *user_fn, int commute, MPI_Op *op) {
	return BigMPI_Op_create (user_fn, commute, op);
}
inline int _MPI_Type_create_struct(MPI_Count count,
                             const MPI_Count array_of_blocklengths[],
                             const MPI_Count array_of_displacements[],
                             const MPI_Datatype array_of_types[],
                             MPI_Datatype *newtype) {
	// TODO not sure what to do here
}
inline int _MPI_Type_indexed(MPI_Count count, const MPI_Count array_of_blocklengths[],
                       const MPI_Count array_of_displacements[],
                       MPI_Datatype oldtype, MPI_Datatype *newtype) {
	// TODO not sure what to do here
}
#else // no BIGMPI
/* If there is no BigMPI, but there is
 * an MPI implementation v4.x use
 * that for large tranfers. Caveat:
 * there are some incomplete v4.x
 * implementations around which will fail.
 * Hopefully they will be fixed
 * before this goes to production
 */
#if MPI_VERSION==4
#define MpiInt MPI_Count
inline int _MPI_Isend(const void *buf, MPI_Count count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, MPI_Request *request) {
	return MPI_Isend_c(buf, count, datatype, dest, tag, comm, request);
}
inline int _MPI_Irecv(void *buf, MPI_Count count, MPI_Datatype datatype,
		int source, int tag, MPI_Comm comm, MPI_Request *request) {
	return MPI_Irecv_c (buf, count, datatype, source, tag, comm, request);
}
inline int _MPI_Op_create(MPI_User_function_c *user_fn, int commute, MPI_Op *op) {
	return MPI_Op_create_c (user_fn, commute, op);
}
inline int _MPI_Type_create_struct(MPI_Count count,
                             const MPI_Count array_of_blocklengths[],
                             const MPI_Count array_of_displacements[],
                             const MPI_Datatype array_of_types[],
                             MPI_Datatype *newtype) {
	return MPI_Type_create_struct_c(count, array_of_blocklengths,
                             array_of_displacements, array_of_types,
                             newtype);
}
inline int _MPI_Type_indexed(MPI_Count count, const MPI_Count array_of_blocklengths[],
                       const MPI_Count array_of_displacements[],
                       MPI_Datatype oldtype, MPI_Datatype *newtype) {
	return MPI_Type_indexed_c(count, array_of_blocklengths, array_of_displacements,
                       oldtype, newtype);
}
/* Otherwise default to non-big
 * transfers
 */
#else // no MPI v4
#define MpiInt int
inline int _MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest,
              int tag, MPI_Comm comm, MPI_Request *request) {
	return MPI_Isend(buf, count, datatype, dest, tag, comm, request);
}
inline int _MPI_Irecv(void *buf, int count, MPI_Datatype datatype,
		int source, int tag, MPI_Comm comm, MPI_Request *request) {
	return MPI_Irecv (buf, count, datatype, source, tag, comm, request);
}
inline int _MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) {
	return MPI_Op_create (user_fn, commute, op);
}
inline int _MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[],
                           MPI_Datatype *newtype)
	return MPI_Type_create_struct(count, array_of_blocklengths,
                             array_of_displacements, array_of_types,
                             newtype);
}
int MPI_Type_indexed(int count, const int array_of_blocklengths[],
                     const int array_of_displacements[], MPI_Datatype oldtype,
                     MPI_Datatype *newtype) {
	return MPI_Type_indexed(count, array_of_blocklengths, array_of_displacements,
                       oldtype, newtype);
}
#endif // MPI version
#endif // BIGMPI

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
