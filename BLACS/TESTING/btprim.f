      SUBROUTINE BTSETUP( MEM, MEMLEN, CMEM, CMEMLEN, OUTNUM,
     $                    TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX,
     $                    IAM, NNODES )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      LOGICAL TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX
      INTEGER MEMLEN, CMEMLEN, OUTNUM, IAM, NNODES
*     ..
*     .. Array Arguments ..
      INTEGER MEM(MEMLEN)
      CHARACTER*1 CMEM(CMEMLEN)
*     ..
*
*  Purpose
*  =======
*  BTSETUP:  Sets up communicator and initiliazes MPI if needed.
*
*  ====================================================================
*
*     ..
*     .. Local Scalars
      LOGICAL INIT
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*     ..
*     .. Executable Statements ..
*
      IERR = 0
      CALL MPI_INITIALIZED(INIT, IERR)
      IF (.NOT.INIT) CALL MPI_INIT(IERR)
      IF (IERR.NE.0) CALL BTMPIERR("mpi_init", IERR)
      CALL MPI_COMM_DUP(MPI_COMM_WORLD, BTCOMM, IERR)
      IF (IERR.NE.0) CALL BTMPIERR("MPI_COMM_DUP", IERR)
*
      RETURN
      END
      INTEGER FUNCTION IBTMYPROC()
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*  Purpose
*  =======
*  IBTMYPROC: returns a process number between 0 .. NPROCS-1.  On
*  systems not natively in this numbering scheme, translates to it.
*
*  ====================================================================
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Local Scalars ..
      INTEGER RANK
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*     ..
*     .. Executable Statements ..
*
      CALL MPI_COMM_RANK(BTCOMM, RANK, IERR)
      IF (IERR.NE.0) CALL BTMPIERR("MPI_COMM_RANK", IERR)
      IBTMYPROC = RANK
      RETURN
*
*     End of IBTMYPROC
*
      END
*
      INTEGER FUNCTION IBTNPROCS()
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*  Purpose
*  =======
*  IBTNPROCS: returns the number of processes in the machine.
*
*  ====================================================================
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Local Scalars ..
      INTEGER NPROC
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*     ..
*     .. Executable Statements ..
*
      CALL MPI_COMM_SIZE(BTCOMM, NPROC, IERR)
      IF (IERR.NE.0) CALL BTMPIERR("MPI_COMM_SIZE", IERR)
      IBTNPROCS = NPROC
*
      RETURN
*
*     End of IBTNPROCS
*
      END
*
      SUBROUTINE BTSEND(DTYPE, N, BUFF, DEST, MSGID)
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      INTEGER N, DTYPE, DEST, MSGID
*     ..
*     .. Array Arguments ..
      REAL BUFF(*)
*     ..
*
*     PURPOSE
*     =======
*     BTSEND: Communication primitive used to send messages independent
*     of the BLACS.  May safely be either locally or globally blocking.
*
*     Arguments
*     =========
*     DTYPE    (input) INTEGER
*              Indicates what data type BUFF is (same as PVM):
*                1  =  RAW BYTES
*                3  =  INTEGER
*                4  =  SINGLE PRECISION REAL
*                6  =  DOUBLE PRECISION REAL
*                5  =  SINGLE PRECISION COMPLEX
*                7  =  DOUBLE PRECISION COMPLEX
*
*     N        (input) INTEGER
*              The number of elements of type DTYPE in BUFF.
*
*     BUFF     (input) accepted as INTEGER array
*              The array to be communicated.  Its true data type is
*              indicated by DTYPE.
*
*     DEST      (input) INTEGER
*               The destination of the message.
*
*     MSGID     (input) INTEGER
*               The message ID (AKA message tag or type).
*
* =====================================================================
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. Local Scalars ..
      INTEGER I, IAM, MPIDTYPE
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*
      IF( DTYPE .EQ. 1 ) THEN
         MPIDTYPE = MPI_BYTE
      ELSE IF( DTYPE .EQ. 3 ) THEN
         MPIDTYPE = MPI_INTEGER
      ELSE IF( DTYPE .EQ. 4 ) THEN
         MPIDTYPE = MPI_REAL
      ELSE IF( DTYPE .EQ. 5 ) THEN
         MPIDTYPE = MPI_COMPLEX
      ELSE IF( DTYPE .EQ. 6 ) THEN
         MPIDTYPE = MPI_DOUBLE_PRECISION
      ELSE IF( DTYPE .EQ. 7 ) THEN
         MPIDTYPE = MPI_DOUBLE_COMPLEX
      END IF
*
*     Send the message
*
      IF( DEST .EQ. -1 ) THEN
         IAM = IBTMYPROC()
         DO 10 I = 0, IBTNPROCS()-1
            IF( I .NE. IAM ) THEN
               CALL MPI_SEND(BUFF, N, MPIDTYPE, I, 0, BTCOMM, IERR)
               IF (IERR.NE.0) CALL BTMPIERR("MPI_SEND", IERR)
            END IF
   10    CONTINUE
      ELSE
         CALL MPI_SEND(BUFF, N, MPIDTYPE, DEST, 0, BTCOMM, IERR)
         IF (IERR.NE.0) CALL BTMPIERR("MPI_SEND", IERR)
      END IF
*
      RETURN
*
*     End BTSEND
*
      END
*
      SUBROUTINE BTRECV(DTYPE, N, BUFF, SRC, MSGID)
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER N, DTYPE, SRC, MSGID
*     ..
*     .. Array Arguments ..
      REAL BUFF(*)
*     ..
*
*     PURPOSE
*     =======
*     BTRECV: Globally blocking receive.
*
*     Arguments
*     =========
*     DTYPE    (input) INTEGER
*              Indicates what data type BUFF is:
*                1  =  RAW BYTES
*                3  =  INTEGER
*                4  =  SINGLE PRECISION REAL
*                6  =  DOUBLE PRECISION REAL
*                5  =  SINGLE PRECISION COMPLEX
*                7  =  DOUBLE PRECISION COMPLEX
*
*     N        (input) INTEGER
*              The number of elements of type DTYPE in BUFF.
*
*     BUFF     (output) INTEGER
*              The buffer to receive into.
*
*     SRC      (input) INTEGER
*              The source of the message.
*
*     MSGID    (input) INTEGER
*              The message ID.
*
* =====================================================================
*     ..
*     .. Local Scalars ..
      INTEGER MPIDTYPE
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Local Arrays ..
      INTEGER STAT(MPI_STATUS_SIZE)
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*
      IF( DTYPE .EQ. 1 ) THEN
         MPIDTYPE = MPI_BYTE
      ELSE IF( DTYPE .EQ. 3 ) THEN
         MPIDTYPE = MPI_INTEGER
      ELSE IF( DTYPE .EQ. 4 ) THEN
         MPIDTYPE = MPI_REAL
      ELSE IF( DTYPE .EQ. 5 ) THEN
         MPIDTYPE = MPI_COMPLEX
      ELSE IF( DTYPE .EQ. 6 ) THEN
         MPIDTYPE = MPI_DOUBLE_PRECISION
      ELSE IF( DTYPE .EQ. 7 ) THEN
         MPIDTYPE = MPI_DOUBLE_COMPLEX
      END IF
*
      CALL MPI_RECV( BUFF, N, MPIDTYPE, SRC, 0, BTCOMM, STAT, IERR )
      IF (IERR.NE.0) CALL BTMPIERR("MPI_RECV", IERR)
*
      RETURN
*
*     End of BTRECV
*
      END
*
      INTEGER FUNCTION IBTSIZEOF(TYPE)
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 TYPE
*     ..
*
*  Purpose
*  =======
*  IBTSIZEOF: Returns the size, in bytes, of the 5 data types.
*  If your platform has a different size for DOUBLE PRECISION, you must
*  change the parameter statement in BLACSTEST as well.
*
*  Arguments
*  =========
*  TYPE     (input) CHARACTER*1
*           The data type who's size is to be determined:
*           'I' : INTEGER
*           'S' : SINGLE PRECISION REAL
*           'D' : DOUBLE PRECISION REAL
*           'C' : SINGLE PRECISION COMPLEX
*           'Z' : DOUBLE PRECISION COMPLEX
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  LSAME
      EXTERNAL LSAME
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*     ..
*     .. Local Scalars ..
      INTEGER LENGTH
      LOGICAL INIT
      DATA INIT /.FALSE./
*     ..
*     .. Executable Statements ..
*
*
*     Initialize MPI, if necessary
*
      IF (.NOT.INIT) THEN
         CALL MPI_INITIALIZED(INIT, IERR)
         IF (.NOT.INIT) CALL MPI_INIT(IERR)
         IF (IERR.NE.0) CALL BTMPIERR("mpi_init", IERR)
         INIT = .TRUE.
      END IF
*
      IF( LSAME(TYPE, 'I') ) THEN
         CALL MPI_TYPE_SIZE( MPI_INTEGER, LENGTH, IERR )
         IF (IERR.NE.0) CALL BTMPIERR("MPI_TYPE_SIZE", IERR)
      ELSE IF( LSAME(TYPE, 'S') ) THEN
         CALL MPI_TYPE_SIZE( MPI_REAL, LENGTH, IERR )
         IF (IERR.NE.0) CALL BTMPIERR("MPI_TYPE_SIZE", IERR)
      ELSE IF( LSAME(TYPE, 'D') ) THEN
         CALL MPI_TYPE_SIZE( MPI_DOUBLE_PRECISION, LENGTH, IERR )
         IF (IERR.NE.0) CALL BTMPIERR("MPI_TYPE_SIZE", IERR)
      ELSE IF( LSAME(TYPE, 'C') ) THEN
         CALL MPI_TYPE_SIZE( MPI_COMPLEX, LENGTH, IERR )
         IF (IERR.NE.0) CALL BTMPIERR("MPI_TYPE_SIZE", IERR)
      ELSE IF( LSAME(TYPE, 'Z') ) THEN
         CALL MPI_TYPE_SIZE( MPI_DOUBLE_COMPLEX, LENGTH, IERR )
         IF (IERR.NE.0) CALL BTMPIERR("MPI_TYPE_SIZE", IERR)
      END IF
      IBTSIZEOF = LENGTH
*
      RETURN
      END
      SUBROUTINE BTMPIERR(ROUT, IERR0)
      CHARACTER*(*) ROUT
      INTEGER IERR0
*     ..
*     .. Include Files ..
      INCLUDE 'mpif.h'
*     ..
*     .. Common Blocks ..
      COMMON /BTMPI/ BTCOMM, IERR
      INTEGER BTCOMM, IERR
*
      WRITE(*,1000) ROUT, IERR
      CALL MPI_ABORT(BTCOMM, IERR0, IERR)
*
 1000 FORMAT('Error #',I20,' from routine ',A)
      RETURN
      END
