

/*
 * These macros define how C routines will be called.  ADD_ assumes that
 * they will be called by fortran, which expects C routines to have an
 * underscore postfixed to the name (Suns, and the Intel expect this).
 * NOCHANGE indicates that fortran will be calling, and that it expects
 * the name called by fortran to be identical to that compiled by the C
 * (RS6K's do this).  UPCASE says it expects C routines called by fortran
 * to be in all upcase (CRAY wants this). 
 */

#define ADD_       0
#define NOCHANGE   1
#define UPCASE     2
#define C_CALL     3

#ifdef UpCase
#define F77_CALL_C UPCASE
#endif

#ifdef NoChange
#define F77_CALL_C NOCHANGE
#endif

#ifdef Add_
#define F77_CALL_C ADD_
#endif

#ifndef F77_CALL_C
#define F77_CALL_C ADD_
#endif

#if (F77_CALL_C == ADD_)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine
 * No redefinition necessary to have following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void pdgemm_(...)
 *
 * This is the default.
 */

#endif

#if (F77_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void PDGEMM(...)
 */
                                                            /* TOOLS */
#define pdlasnbt_           PDLASNBT
#define pdlachkieee_        PDLACHKIEEE
#define pdlaiectl_          PDLAIECTL
#define pdlaiectb_          PDLAIECTB

#define pslasnbt_           PSLASNBT
#define pslachkieee_        PSLACHKIEEE
#define pslaiect_           PSLAIECT

#endif

#if (F77_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran 77
 * routine call a C routine 
 * for following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call pdgemm(...)           void pdgemm(...)
 */
                                                            /* TOOLS */
#define pdlasnbt_           pdlasnbt
#define pdlachkieee_        pdlachkieee
#define pdlaiectl_          pdlaiectl
#define pdlaiectb_          pdlaiectb

#define pslasnbt_           pslasnbt
#define pslachkieee_        pslachkieee
#define pslaiect_           pslaiect
#endif
