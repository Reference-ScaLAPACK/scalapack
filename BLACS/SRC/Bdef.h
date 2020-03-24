#ifndef BDEF_H
#define BDEF_H 1

/*
 * Include the system dependant and user defined stuff
 */
#include "Bconfig.h"

/*
 * Data type defining a scope for the BLACS
 */
typedef struct bLaCsScOpE BLACSSCOPE;
struct bLaCsScOpE
{
   MPI_Comm comm;
   Int ScpId, MaxId, MinId;
   Int Np, Iam;
};
/*
 * Data type defining a context for the BLACS
 */
typedef struct bLaCsCoNtExT BLACSCONTEXT;
struct bLaCsCoNtExT
{
   BLACSSCOPE rscp, cscp, ascp, pscp; /* row, column, all, and pt2pt scopes */
   BLACSSCOPE *scp;                   /* pointer to present scope */
   Int TopsRepeat;                    /* Use only repeatable topologies? */
   Int TopsCohrnt;                    /* Use only coherent topologies? */
   Int Nb_bs, Nr_bs;           /* for bcast general tree and multiring tops */
   Int Nb_co, Nr_co;           /* for combine general tree and multiring tops */
};

/*
 *  Define the fortran data types COMPLEX*8 (SCOMPLEX)
 *  and COMPLEX*16 (DCOMPLEX).
 */
typedef struct {double r, i;} DCOMPLEX;
typedef struct {float r, i;} SCOMPLEX;

/*
 *  These variables will be defined to be MPI datatypes for complex and double
 *  complex if we are using the C interface to MPI.  If we use the fortran
 *  interface, we need to declare the contants array.  I'm too lazy to declare
 *  these guys external in every file that needs them.
 */
#ifndef GlobalVars
   extern Int *BI_COMM_WORLD;
#endif

/*
 *  Definition of buffer type for BLACS' asynchronous operations
 */
typedef struct bLaCbUfF BLACBUFF;
struct bLaCbUfF
{
   char *Buff;             /* send/recv buffer */
   Int Len;                /* length of buffer in bytes */
   Int nAops;              /* number of asynchronous operations out of buff */
   MPI_Request *Aops;   /* list of async. operations out of buff */
   MPI_Datatype dtype;  /* data type of buffer */
   Int N;                  /* number of elements of data type in buff */
   BLACBUFF *prev, *next;  /* pointer to the other BLACBUFF in queue */
};

/*
 * Pointer to the combine's vector-vector functions
 */
typedef void (*VVFUNPTR)(Int, char *, char *);
typedef void (*SDRVPTR)(BLACSCONTEXT *, Int, Int, BLACBUFF *);


#define BI_DistType                  unsigned short
#define BI_MpiDistType               MPI_UNSIGNED_SHORT

#define BUFFALIGN    8      /* force all buffers to 8 byte alignment */
#define BANYNODE     MPI_ANY_SOURCE
#define PT2PTID      9976   /* TAG used for point to point */
#define NOTINCONTEXT -1  /* Indicates node called gridmap, but not in grid */
#define MAXNCTXT     10      /* initial guess at max # of contexts */
#define MAXNSYSCTXT  10   /* initial guess at max # of system context */
#define AOPDONE      MPI_REQUEST_NULL
#define BUFWAIT      120      /* Length of time to wait for emergency buff */

/*
 * Error codes
 */
#define NORV 1          /* No receiver (only 1 proc in scoped op) */
#define NPOW2 2         /* Number of procs is not a power of 2 */
#define BADSCP 3        /* Scope not row, column or all */

/*
 * Data types
 */
#define INTEGER   3
#define SINGLE    4
#define DOUBLE    6
#define COMPLEX8  5
#define COMPLEX16 7

#define FULLCON 0      /* top is fully connected */

/*
 * Routine types
 */
#define RT_SD    1
#define RT_RV    2
#define RT_BS    3
#define RT_BR    4
#define RT_COMB  5

/*
 * Legal WHAT values for BLACS_SET
 */
#define SGET_SYSCONTXT    0
#define SGET_MSGIDS       1
#define SGET_DEBUGLVL     2
#define SGET_BLACSCONTXT 10
#define SGET_NR_BS       11
#define SGET_NB_BS       12
#define SGET_NR_CO       13
#define SGET_NB_CO       14
#define SGET_TOPSREPEAT  15
#define SGET_TOPSCOHRNT  16

/*
 * These are prototypes for error and warning functions -- I don't want
 * to prototype them in each routine.
 */
void BI_BlacsWarn(Int ConTxt, Int line, char *file, char *form, ...);
void BI_BlacsErr(Int ConTxt, Int line, char *file, char *form, ...);
Int BI_ContxtNum(BLACSCONTEXT *ctxt);

/*
 * If we've got an ANSI standard C compiler, we can use void pointers...
 */
#define BVOID void


/*
 * ========================================================================
 *     PREPROCESSOR MACRO FUNCTIONS USED FOR OPTIMIZATION & CONVENIENCE
 * ========================================================================
 */

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

/*
 * Slightly modified gridinfo substitute
 */
#define Mgridinfo(ctxt, Ng0, nprow0, npcol0, myrow0, mycol0)\
{\
   (Ng0) = (ctxt)->ascp.Np;\
   (nprow0) = (ctxt)->cscp.Np;\
   (npcol0) = (ctxt)->rscp.Np;\
   (myrow0) = (ctxt)->cscp.Iam;\
   (mycol0) = (ctxt)->rscp.Iam;\
}

/*
 * These routines return coordinates based on nodes number, or node number
 * based on coordinates.  Those routines with v after the M return virtual
 * nodes numbers (i.e., in respect to the grid, not physical node numbers)
 * based on grid coordinates, or grid coordinates based on virtual node numbers.
 */
#define Mpcoord(ctxt, node, prow, pcol)\
{\
   (prow) = (node) / (ctxt)->rscp.Np;\
   (pcol) = (node) % (ctxt)->rscp.Np;\
}
#define Mvpcoord(ctxt, node, prow, pcol) \
        Mpcoord((ctxt), (node), (prow), (pcol));

#define Mkpnum(ctxt, prow, pcol)  ( (prow)*(ctxt)->rscp.Np+(pcol) )
#define Mvkpnum(ctxt, prow, pcol) ( (prow)*(ctxt)->rscp.Np+(pcol) )

/*
 * This macro returns scoped message ID's.
 */
#define Mscopeid(ctxt) (ctxt)->scp->ScpId; \
   if (++(ctxt)->scp->ScpId == (ctxt)->scp->MaxId) \
      (ctxt)->scp->ScpId = (ctxt)->scp->MinId;

/*
 *  Get context, and check for validity if debug level is high
 */
#if (BlacsDebugLvl > 0)
#define MGetConTxt(Context, ctxtptr)\
{\
   extern BLACSCONTEXT **BI_MyContxts;\
   extern Int BI_MaxNCtxt;\
   if ( ((Context) >= BI_MaxNCtxt) || ((Context) < 0) )\
      BI_BlacsErr(-1, __LINE__, __FILE__, "Invalid context handle: %d",\
                  (Context));\
   else if (BI_MyContxts[(Context)] == NULL)\
      BI_BlacsErr(-1, __LINE__, __FILE__, "Invalid context, handle=%d",\
                  (Context));\
   (ctxtptr) = BI_MyContxts[(Context)];\
}
#else
#define MGetConTxt(Context, ctxtptr)\
{\
   extern BLACSCONTEXT **BI_MyContxts;\
   (ctxtptr) = BI_MyContxts[(Context)];\
}
#endif
/*
 * This macro handles MPI errors
 */
#if(BlacsDebugLvl > 0)
#define Mmpierror(ierr, rout, ctxt, line, file) \
{ \
   if ( (ierr) != BI_MPI_SUCCESS )\
      BI_BlacsErr(BI_ContxtNum((ctxt)), (line), (file), \
                  "MPI error %d on call to %s", (ierr), (rout)); \
}
#else
#define Mmpierror(ierr, rout, ctxt, line, file)
#endif
/*
 * A small macro useful for debugging
 */
#define ErrPrint \
{ \
   extern Int BI_Iam; \
   fprintf(stderr, "%d: line %d of file %s\n", BI_Iam, __LINE__, __FILE__); \
}

/*
 * These macros allow for the funky function declarations and character handling
 * needed on the CRAY to have a C routine callable from fortran
 */
#define F_VOID_FUNC void
#define F_INT_FUNC  Int
#define F_DOUBLE_FUNC double

#if (INTFACE == C_CALL)

#define F2C_CharTrans(c) *(c)

#else

#ifdef CRAY
#define F2C_CharTrans(c) *( _fcdtocp((c)) )
#define F_CHAR      _fcd
#else
#define F2C_CharTrans(c) *(c)
#define F_CHAR      char *
#endif

#endif

/*
 *  These macros allow for accessing values and addresses of parameters, which
 *  will be pointers if we're using fortran, and values if we're using C.
 */
#if (INTFACE == C_CALL)
#define Mpval(para) (para)
#define Mpaddress(para) (&(para))
#define Mwalltime Cdwalltime00
#else
#define Mpval(para) (*(para))
#define Mpaddress(para) (para)
#define Mwalltime dwalltime00_
#endif

/*
 * Real and complex absolute values
 */
#define Rabs(x) ( (x) < 0 ? (x) * -1 : (x) )
#define Cabs(z) ( (((z).i) < 0 ? ((z).i) * -1 : ((z).i)) + (((z).r) < 0 ? ((z).r) * -1 : ((z).r)) )

/*
 * Figures the length of packed trapezoidal matrix
 */
#define trsize(diag, m, n, bytes, length)\
{\
   if ( (diag) == 'u' ) (length) = 1;\
   else (length) = 0;\
   if ( (m) > (n) )\
      (length) = ( (n) * ( (m) - (n) ) + ( (n)*(n) ) - ( (n)*(n) )/2 +\
                   (n)/2 - (n) * (length) ) * (bytes);\
   else\
      (length) = ( (m) * ( (n) - (m) ) + ( (m)*(m) ) - ( (m)*(m) )/2 +\
                   (m)/2 - (m) * (length) ) * (bytes);\
}

/*
 * These macros call the correct packing/unpacking routines
 */
#define BI_cmvcopy(m, n, A, lda, buff) \
        BI_smvcopy(2*(m), (n), (float *) (A), 2*(lda), (float *) (buff))
#define BI_cvmcopy(m, n, A, lda, buff) \
        BI_svmcopy(2*(m), (n), (float *) (A), 2*(lda), (float *) (buff))
#define BI_zmvcopy(m, n, A, lda, buff) \
        BI_dmvcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))
#define BI_zvmcopy(m, n, A, lda, buff) \
        BI_dvmcopy(2*(m), (n), (double *) (A), 2*(lda), (double *) (buff))

/*
 * This macro avoids freeing types when the zero-byte workaround was applied
 */
#ifdef ZeroByteTypeBug
#define BI_MPI_TYPE_FREE(t) (*(t) != MPI_BYTE ? MPI_Type_free(t) : 0)
#else
#define BI_MPI_TYPE_FREE(t) MPI_Type_free(t)
#endif

#if (FORTRAN_CALL_C == NOCHANGE)
/*
 * These defines set up the naming scheme required to have a fortran
 * routine call a C routine (which is what the BLACS are written in)
 * for the following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgebs2d(...)          void dgebs2d(...)
 */

/*
 * Support routines
 */
#define blacs_pinfo_                   blacs_pinfo
#define blacs_setup_                   blacs_setup
#define setpvmtids_                    setpvmtids
#define blacs_set_                     blacs_set
#define blacs_get_                     blacs_get
#define blacs_gridinit_                blacs_gridinit
#define blacs_gridmap_                 blacs_gridmap
#define ksendid_                       ksendid
#define krecvid_                       krecvid
#define kbsid_                         kbsid
#define kbrid_                         kbrid
#define blacs_freebuff_                blacs_freebuff
#define blacs_gridexit_                blacs_gridexit
#define blacs_abort_                   blacs_abort
#define blacs_exit_                    blacs_exit
#define blacs_gridinfo_                blacs_gridinfo
#define blacs_pnum_                    blacs_pnum
#define blacs_pcoord_                  blacs_pcoord
#define dcputime00_                    dcputime00
#define dwalltime00_                   dwalltime00
#define blacs_barrier_                 blacs_barrier

/*
 * Main, type dependent, routines
 */
#define igesd2d_   igesd2d
#define igerv2d_   igerv2d
#define igebs2d_   igebs2d
#define igebr2d_   igebr2d
#define itrsd2d_   itrsd2d
#define itrrv2d_   itrrv2d
#define itrbs2d_   itrbs2d
#define itrbr2d_   itrbr2d
#define igsum2d_   igsum2d
#define igamx2d_   igamx2d
#define igamn2d_   igamn2d
#define sgesd2d_   sgesd2d
#define sgerv2d_   sgerv2d
#define sgebs2d_   sgebs2d
#define sgebr2d_   sgebr2d
#define strsd2d_   strsd2d
#define strrv2d_   strrv2d
#define strbs2d_   strbs2d
#define strbr2d_   strbr2d
#define sgsum2d_   sgsum2d
#define sgamx2d_   sgamx2d
#define sgamn2d_   sgamn2d
#define dgesd2d_   dgesd2d
#define dgerv2d_   dgerv2d
#define dgebs2d_   dgebs2d
#define dgebr2d_   dgebr2d
#define dtrsd2d_   dtrsd2d
#define dtrrv2d_   dtrrv2d
#define dtrbs2d_   dtrbs2d
#define dtrbr2d_   dtrbr2d
#define dgsum2d_   dgsum2d
#define dgamx2d_   dgamx2d
#define dgamn2d_   dgamn2d
#define cgesd2d_   cgesd2d
#define cgerv2d_   cgerv2d
#define cgebs2d_   cgebs2d
#define cgebr2d_   cgebr2d
#define ctrsd2d_   ctrsd2d
#define ctrrv2d_   ctrrv2d
#define ctrbs2d_   ctrbs2d
#define ctrbr2d_   ctrbr2d
#define cgsum2d_   cgsum2d
#define cgamx2d_   cgamx2d
#define cgamn2d_   cgamn2d
#define zgesd2d_   zgesd2d
#define zgerv2d_   zgerv2d
#define zgebs2d_   zgebs2d
#define zgebr2d_   zgebr2d
#define ztrsd2d_   ztrsd2d
#define ztrrv2d_   ztrrv2d
#define ztrbs2d_   ztrbs2d
#define ztrbr2d_   ztrbr2d
#define zgsum2d_   zgsum2d
#define zgamx2d_   zgamx2d
#define zgamn2d_   zgamn2d

#elif (FORTRAN_CALL_C == UPCASE)
/*
 * These defines set up the naming scheme required to have a fortran
 * routine call a C routine (which is what the BLACS are written in)
 * for the following Fortran to C interface:
 *           FORTRAN CALL               C DECLARATION
 *           call dgebs2d(...)          void DGEBS2D(...)
 */
/*
 * Support routines
 */
#define blacs_pinfo_                   BLACS_PINFO
#define blacs_setup_                   BLACS_SETUP
#define setpvmtids_                    SETPVMTIDS
#define blacs_set_                     BLACS_SET
#define blacs_get_                     BLACS_GET
#define blacs_gridinit_                BLACS_GRIDINIT
#define blacs_gridmap_                 BLACS_GRIDMAP
#define ksendid_                       KSENDID
#define krecvid_                       KRECVID
#define kbsid_                         KBSID
#define kbrid_                         KBRID
#define blacs_freebuff_                BLACS_FREEBUFF
#define blacs_gridexit_                BLACS_GRIDEXIT
#define blacs_abort_                   BLACS_ABORT
#define blacs_exit_                    BLACS_EXIT
#define blacs_gridinfo_                BLACS_GRIDINFO
#define blacs_pnum_                    BLACS_PNUM
#define blacs_pcoord_                  BLACS_PCOORD
#define dcputime00_                    DCPUTIME00
#define dwalltime00_                   DWALLTIME00
#define blacs_barrier_                 BLACS_BARRIER

/*
 * Main, type dependent, routines
 */
#define igesd2d_   IGESD2D
#define igerv2d_   IGERV2D
#define igebs2d_   IGEBS2D
#define igebr2d_   IGEBR2D
#define itrsd2d_   ITRSD2D
#define itrrv2d_   ITRRV2D
#define itrbs2d_   ITRBS2D
#define itrbr2d_   ITRBR2D
#define igsum2d_   IGSUM2D
#define igamx2d_   IGAMX2D
#define igamn2d_   IGAMN2D
#define sgesd2d_   SGESD2D
#define sgerv2d_   SGERV2D
#define sgebs2d_   SGEBS2D
#define sgebr2d_   SGEBR2D
#define strsd2d_   STRSD2D
#define strrv2d_   STRRV2D
#define strbs2d_   STRBS2D
#define strbr2d_   STRBR2D
#define sgsum2d_   SGSUM2D
#define sgamx2d_   SGAMX2D
#define sgamn2d_   SGAMN2D
#define dgesd2d_   DGESD2D
#define dgerv2d_   DGERV2D
#define dgebs2d_   DGEBS2D
#define dgebr2d_   DGEBR2D
#define dtrsd2d_   DTRSD2D
#define dtrrv2d_   DTRRV2D
#define dtrbs2d_   DTRBS2D
#define dtrbr2d_   DTRBR2D
#define dgsum2d_   DGSUM2D
#define dgamx2d_   DGAMX2D
#define dgamn2d_   DGAMN2D
#define cgesd2d_   CGESD2D
#define cgerv2d_   CGERV2D
#define cgebs2d_   CGEBS2D
#define cgebr2d_   CGEBR2D
#define ctrsd2d_   CTRSD2D
#define ctrrv2d_   CTRRV2D
#define ctrbs2d_   CTRBS2D
#define ctrbr2d_   CTRBR2D
#define cgsum2d_   CGSUM2D
#define cgamx2d_   CGAMX2D
#define cgamn2d_   CGAMN2D
#define zgesd2d_   ZGESD2D
#define zgerv2d_   ZGERV2D
#define zgebs2d_   ZGEBS2D
#define zgebr2d_   ZGEBR2D
#define ztrsd2d_   ZTRSD2D
#define ztrrv2d_   ZTRRV2D
#define ztrbs2d_   ZTRBS2D
#define ztrbr2d_   ZTRBR2D
#define zgsum2d_   ZGSUM2D
#define zgamx2d_   ZGAMX2D
#define zgamn2d_   ZGAMN2D

#elif (FORTRAN_CALL_C == FCISF2C)
/*
 * These defines set up the naming scheme required to have a fortran
 * routine call a C routine (which is what the BLACS are written in)
 * for systems where the fortran "compiler" is actually f2c (a fortran
 * to C conversion utility).
 */
/*
 * Initialization routines
 */
#define blacs_pinfo_    blacs_pinfo__
#define blacs_setup_    blacs_setup__
#define blacs_set_      blacs_set__
#define blacs_get_      blacs_get__
#define blacs_gridinit_ blacs_gridinit__
#define blacs_gridmap_  blacs_gridmap__
/*
 * Destruction routines
 */
#define blacs_freebuff_ blacs_freebuff__
#define blacs_gridexit_ blacs_gridexit__
#define blacs_abort_    blacs_abort__
#define blacs_exit_     blacs_exit__
/*
 * Informational & misc.
 */
#define blacs_gridinfo_ blacs_gridinfo__
#define blacs_pnum_     blacs_pnum__
#define blacs_pcoord_   blacs_pcoord__
#define blacs_barrier_  blacs_barrier__

#endif


#endif
