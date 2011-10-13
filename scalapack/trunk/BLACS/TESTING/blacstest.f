      PROGRAM BLACSTEST
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*  Purpose
*  =======
*  This is the driver for the BLACS test suite.
*
*  Arguments
*  =========
*  None.  Input is done via the data files indicated below.
*
*  Input Files
*  ===========
*  The following input files must reside in the current working
*  directory:
*
*  bt.dat   -- input parameters for the test run as a whole
*  sdrv.dat -- input parameters for point-to-point testing
*  bsbr.dat -- input parameters for broadcast testing
*  comb.dat -- input parameters for combine testing
*
*  Output Files
*  ============
*  Test results are generated and sent to output file as
*  specified by the user in bt.dat.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER CMEMSIZ, MEMELTS
      PARAMETER( MEMELTS = 250000 )
      PARAMETER( CMEMSIZ = 10000 )
*     ..
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER IBTMSGID, IBTSIZEOF
      REAL SBTEPS
      DOUBLE PRECISION DBTEPS
      EXTERNAL ALLPASS, IBTMSGID, SBTEPS, DBTEPS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_PINFO, BTSETUP, RDBTIN
*     ..
*     .. Local Scalars ..
      INTEGER I, IAM, NNODES, VERB, OUTNUM, MEMLEN, NPREC, ISIZE, DSIZE
      LOGICAL TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX
*     ..
*     .. Local Arrays ..
      CHARACTER*1 CMEM(CMEMSIZ), PREC(9)
      INTEGER IPREC(9), ITMP(2)
      DOUBLE PRECISION MEM(MEMELTS)
*     ..
*     .. Executable Statements ..
*
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
*
*     Get initial process information, and initialize message IDs
*
      CALL BLACS_PINFO( IAM, NNODES )
      ITMP(1) = IBTMSGID()
*
*     Call BLACS_GRIDINIT so BLACS set up some system stuff:  should
*     make it possible for the user to print, read input files, etc.
*
      IF( NNODES .GT. 0 ) THEN
         CALL BLACS_GET( 0, 0, ITMP )
         CALL BLACS_GRIDINIT(ITMP, 'c', 1, NNODES)
         CALL BLACS_GRIDEXIT(ITMP)
      END IF
*
*     Read in what tests to do
*
      IF( IAM .EQ. 0 )
     $   CALL RDBTIN( TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX, NPREC,
     $               PREC, VERB, OUTNUM )
*
      MEMLEN = (MEMELTS * DSIZE) / ISIZE
*
*     Get process info for communication, and create virtual machine
*     if necessary
*
      CALL BTSETUP( MEM, MEMLEN, CMEM, CMEMSIZ, OUTNUM, TESTSDRV,
     $              TESTBSBR, TESTCOMB, TESTAUX, IAM, NNODES )
*
*     Send out RDBTIN information
*
      IF( IAM .EQ. 0 ) THEN
*
*        Store test info in back of precision array
*
         ITMP(1) = NPREC
         ITMP(2) = VERB
         CALL BTSEND( 3, 2, ITMP, -1, IBTMSGID() )
         DO 10 I = 1, 9
            IPREC(I) = 0
   10    CONTINUE
         DO 20 I = 1, NPREC
            IF( PREC(I) .EQ. 'I' ) THEN
               IPREC(I) = 1
            ELSE IF( PREC(I) .EQ. 'S' ) THEN
               IPREC(I) = 2
            ELSE IF( PREC(I) .EQ. 'D' ) THEN
               IPREC(I) = 3
            ELSE IF( PREC(I) .EQ. 'C' ) THEN
               IPREC(I) = 4
            ELSE IF( PREC(I) .EQ. 'Z' ) THEN
               IPREC(I) = 5
            END IF
   20    CONTINUE
         IF( TESTSDRV ) IPREC(6) = 1
         IF( TESTBSBR ) IPREC(7) = 1
         IF( TESTCOMB ) IPREC(8) = 1
         IF( TESTAUX )  IPREC(9) = 1
         CALL BTSEND( 3, 9, IPREC, -1, IBTMSGID()+1 )
      ELSE
         CALL BTRECV( 3, 2, ITMP, 0, IBTMSGID() )
         NPREC = ITMP(1)
         VERB = ITMP(2)
         CALL BTRECV( 3, 9, IPREC, 0, IBTMSGID()+1 )
         DO 30 I = 1, NPREC
            IF( IPREC(I) .EQ. 1 ) THEN
               PREC(I) = 'I'
            ELSE IF( IPREC(I) .EQ. 2 ) THEN
               PREC(I) = 'S'
            ELSE IF( IPREC(I) .EQ. 3 ) THEN
               PREC(I) = 'D'
            ELSE IF( IPREC(I) .EQ. 4 ) THEN
               PREC(I) = 'C'
            ELSE IF( IPREC(I) .EQ. 5 ) THEN
               PREC(I) = 'Z'
            END IF
   30    CONTINUE
         TESTSDRV = ( IPREC(6) .EQ. 1 )
         TESTBSBR = ( IPREC(7) .EQ. 1 )
         TESTCOMB = ( IPREC(8) .EQ. 1 )
         TESTAUX  = ( IPREC(9) .EQ. 1 )
      ENDIF
*
      IF( TESTSDRV .OR. TESTBSBR .OR. TESTCOMB .OR. TESTAUX ) THEN
*
*        Find maximal machine epsilon for single and double precision
*
         ITMP(1) = INT( SBTEPS() )
         ITMP(1) = INT( DBTEPS() )
*
         CALL RUNTESTS( MEM, MEMLEN, CMEM, CMEMSIZ, PREC, NPREC, OUTNUM,
     $                  VERB, TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX )
*
      END IF
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM,*) ' '
         WRITE(OUTNUM,1000)
         WRITE(OUTNUM,1000)
         IF( ALLPASS(.TRUE.) ) THEN
            WRITE(OUTNUM,2000) 'NO'
         ELSE
            WRITE(OUTNUM,2000) '  '
         END IF
         WRITE(OUTNUM,1000)
         WRITE(OUTNUM,1000)
         IF( OUTNUM.NE.0 .AND. OUTNUM.NE.6 ) CLOSE(OUTNUM)
      ENDIF
*
      CALL BLACS_EXIT(0)
 1000 FORMAT('=======================================')
 2000 FORMAT('THERE WERE ',A2,' FAILURES IN THIS TEST RUN')
      STOP
*
*     End BLACSTESTER
*
      END
*
      SUBROUTINE RUNTESTS( MEM, MEMLEN, CMEM, CMEMLEN, PREC, NPREC,
     $                     OUTNUM, VERB, TESTSDRV, TESTBSBR, TESTCOMB,
     $                     TESTAUX )
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, CMEMLEN, NPREC, OUTNUM, VERB, IAM, NNODES
      LOGICAL TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX
*     ..
*     .. Array Arguments ..
      CHARACTER*1 CMEM(CMEMLEN), PREC(NPREC)
      INTEGER MEM(MEMLEN)
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS, IBTMYPROC, IBTMSGID, IBTSIZEOF, SAFEINDEX
      EXTERNAL IBTNPROCS, IBTMYPROC, IBTMSGID, IBTSIZEOF, SAFEINDEX
*     ..
*     .. External Subroutines ..
      EXTERNAL CSDRVTEST, DSDRVTEST, ISDRVTEST, SSDRVTEST, ZSDRVTEST
      EXTERNAL CBSBRTEST, DBSBRTEST, IBSBRTEST, SBSBRTEST, ZBSBRTEST
      EXTERNAL ISUMTEST, SSUMTEST, DSUMTEST, CSUMTEST, ZSUMTEST
      EXTERNAL IAMXTEST, SAMXTEST, DAMXTEST, CAMXTEST, ZAMXTEST
      EXTERNAL IAMNTEST, SAMNTEST, DAMNTEST, CAMNTEST, ZAMNTEST
      EXTERNAL AUXTEST, BTSEND, BTRECV, BTINFO
*     ..
*     .. Local Scalars ..
      INTEGER NSCOPE, NOP, NTOP, NSHAPE, NMAT, NSRC, NDEST, NGRID
      INTEGER TREP, TCOH, OPPTR, SCOPEPTR, TOPPTR, UPLOPTR, DIAGPTR
      INTEGER MPTR, NPTR, LDSPTR, LDDPTR, LDIPTR
      INTEGER RSRCPTR, CSRCPTR, RDESTPTR, CDESTPTR, PPTR, QPTR
      INTEGER ISEEDPTR, RAPTR, CAPTR, CTXTPTR, WORKPTR, WORKLEN
      INTEGER MEMUSED, CMEMUSED, I, J, K
      INTEGER ISIZE, SSIZE, DSIZE, CSIZE, ZSIZE
*     ..
*     .. Local Arrays ..
      INTEGER ITMP(4)
*     ..
*     .. Executable Statements ..
*
      IAM = IBTMYPROC()
      NNODES = IBTNPROCS()
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
      DSIZE = IBTSIZEOF('D')
      CSIZE = IBTSIZEOF('C')
      ZSIZE = IBTSIZEOF('Z')
*
      IF( IAM.EQ.0 ) THEN
         CALL BLACS_GET( 0, 2, I )
         WRITE(OUTNUM,3000)
         WRITE(OUTNUM,3000)
         WRITE(OUTNUM,2000) I
         WRITE(OUTNUM,3000)
         WRITE(OUTNUM,3000)
      END IF
*
      IF( TESTAUX ) THEN
*
*        Each process will make sure that BLACS_PINFO returns
*        the same value as BLACS_SETUP, and send a packet
*        to node 0 saying whether it was.
*
         CALL BLACS_PINFO( ITMP(1), ITMP(3) )
         CALL BLACS_SETUP( ITMP(2), ITMP(4) )
         IF( IAM .EQ. 0 ) THEN
            DO 35 I = 0, NNODES-1
               IF( I .NE. 0 )
     $            CALL BTRECV( 3, 4, ITMP, I, IBTMSGID()+2 )
               IF( ITMP(1) .NE. ITMP(2) )
     $              WRITE( OUTNUM, 1000 ) ITMP(1), ITMP(2)
               IF( (ITMP(3).NE.ITMP(4)) .OR. (ITMP(3).NE.NNODES) )
     $              WRITE( OUTNUM, 1000 ) ITMP(3), ITMP(4), NNODES
   35       CONTINUE
         ELSE
            CALL BTSEND( 3, 4, ITMP, 0, IBTMSGID()+2 )
         ENDIF
      ENDIF
*
*     Run point-to-point tests as appropriate
*
      IF( TESTSDRV ) THEN
*
*        Get test info
*
         CALL BTINFO( 'SDRV', MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM,
     $                CMEMLEN, OUTNUM, NOP, NSCOPE, TREP, TCOH, NTOP,
     $                NSHAPE, NMAT, NSRC, NGRID, OPPTR, SCOPEPTR,
     $                TOPPTR, UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR,
     $                LDDPTR, LDIPTR, RSRCPTR, CSRCPTR, RDESTPTR,
     $                CDESTPTR, PPTR, QPTR )
*
*        iseedptr used as tests passed/failed array, so it must
*        be of size NTESTS -- It's not used unless VERB < 2
*
         CTXTPTR = MEMUSED + 1
         ISEEDPTR = CTXTPTR + NGRID
         MEMUSED = ISEEDPTR - 1
         IF( VERB .LT. 2 )
     $      MEMUSED = MEMUSED + NSHAPE * NMAT * NSRC * NGRID
*
         CALL MAKEGRIDS( MEM(CTXTPTR), OUTNUM, NGRID, MEM(PPTR),
     $                   MEM(QPTR) )
*
*        Call individual tests as appropriate.
*
         DO 10 I = 1, NPREC
            IF( PREC(I) .EQ. 'I' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, ISIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ISIZE
               CALL ISDRVTEST(OUTNUM, VERB, NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        MEM(RDESTPTR), MEM(CDESTPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'S' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, SSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / SSIZE
               CALL SSDRVTEST(OUTNUM, VERB, NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        MEM(RDESTPTR), MEM(CDESTPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'D' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, DSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / DSIZE
               CALL DSDRVTEST(OUTNUM, VERB, NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        MEM(RDESTPTR), MEM(CDESTPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'C' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, CSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / CSIZE
               CALL CSDRVTEST(OUTNUM, VERB, NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        MEM(RDESTPTR), MEM(CDESTPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'Z' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, ZSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ZSIZE
               CALL ZSDRVTEST(OUTNUM, VERB, NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        MEM(RDESTPTR), MEM(CDESTPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
            END IF
   10    CONTINUE
         CALL FREEGRIDS( NGRID, MEM(CTXTPTR) )
      END IF
*
      IF( TESTBSBR ) THEN
*
*        Get test info
*
         CALL BTINFO( 'BSBR', MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM,
     $                CMEMLEN, OUTNUM, NOP, NSCOPE, TREP, TCOH, NTOP,
     $                NSHAPE, NMAT, NSRC, NGRID, OPPTR, SCOPEPTR,
     $                TOPPTR, UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR,
     $                LDDPTR, LDIPTR, RSRCPTR, CSRCPTR, RDESTPTR,
     $                CDESTPTR, PPTR, QPTR )
*
*        iseedptr used as tests passed/failed array, so it must
*        be of size NTESTS -- It's not used unless VERB < 2
*
         CTXTPTR = MEMUSED + 1
         ISEEDPTR = CTXTPTR + NGRID
         MEMUSED = ISEEDPTR - 1
         IF( VERB .LT. 2 )
     $      MEMUSED = MEMUSED + NSCOPE*NTOP*NSHAPE*NMAT*NSRC*NGRID
*
         CALL MAKEGRIDS( MEM(CTXTPTR), OUTNUM, NGRID, MEM(PPTR),
     $                   MEM(QPTR) )
*
*        Call individual tests as appropriate.
*
         DO 20 I = 1, NPREC
            IF( PREC(I) .EQ. 'I' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, ISIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ISIZE
               CALL IBSBRTEST(OUTNUM, VERB, NSCOPE, CMEM(SCOPEPTR),
     $                        NTOP, CMEM(TOPPTR), NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'S' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, SSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / SSIZE
               CALL SBSBRTEST(OUTNUM, VERB, NSCOPE, CMEM(SCOPEPTR),
     $                        NTOP, CMEM(TOPPTR), NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'D' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, DSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / DSIZE
               CALL DBSBRTEST(OUTNUM, VERB, NSCOPE, CMEM(SCOPEPTR),
     $                        NTOP, CMEM(TOPPTR), NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'C' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, CSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / CSIZE
               CALL CBSBRTEST(OUTNUM, VERB, NSCOPE, CMEM(SCOPEPTR),
     $                        NTOP, CMEM(TOPPTR), NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            ELSE IF( PREC(I) .EQ. 'Z' ) THEN
*
               WORKPTR = SAFEINDEX(MEMUSED + 1, ISIZE, ZSIZE)
               WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ZSIZE
               CALL ZBSBRTEST(OUTNUM, VERB, NSCOPE, CMEM(SCOPEPTR),
     $                        NTOP, CMEM(TOPPTR), NSHAPE, CMEM(UPLOPTR),
     $                        CMEM(DIAGPTR), NMAT, MEM(MPTR),
     $                        MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR),
     $                        NSRC, MEM(RSRCPTR), MEM(CSRCPTR),
     $                        NGRID, MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                        MEM(ISEEDPTR), MEM(WORKPTR), WORKLEN)
*
            END IF
*
   20    CONTINUE
         CALL FREEGRIDS( NGRID, MEM(CTXTPTR) )
      END IF
      IF( TESTCOMB ) THEN
*
*        Get test info
*
         CALL BTINFO( 'COMB', MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM,
     $                CMEMLEN, OUTNUM, NOP, NSCOPE, TREP, TCOH, NTOP,
     $                NSHAPE, NMAT, NDEST, NGRID, OPPTR, SCOPEPTR,
     $                TOPPTR, UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR,
     $                LDDPTR, LDIPTR, RSRCPTR, CSRCPTR, RDESTPTR,
     $                CDESTPTR, PPTR, QPTR )
         CTXTPTR = MEMUSED + 1
         MEMUSED  = CTXTPTR + NGRID - 1
*
*        Find space required by RA and CA arrays
*
         K = 0
         DO 40 J = 0, NOP-1
            IF( CMEM(OPPTR+J).EQ.'>' .OR. CMEM(OPPTR+J).EQ.'<' ) THEN
               DO 30 I = 0, NMAT
*
*                 NOTE: here we assume ipre+ipost = 4*M
*
                  K = MAX0( K, 4*MEM(MPTR+I) )
                  IF ( MEM(LDIPTR+I) .NE. -1 )
     $               K = MAX0( K, MEM(NPTR+I)*MEM(LDIPTR+I) +
     $                            4*MEM(MPTR+I) )
   30          CONTINUE
            END IF
   40    CONTINUE
         RAPTR = MEMUSED + 1
         CAPTR = RAPTR + K
*
*        iseed array also used as tests passed/failed array, so it must
*        be of size MAX( 4*NNODES, NTESTS )
*
         ISEEDPTR = CAPTR + K
         I = 0
         IF( VERB.LT.2 ) I = NSCOPE * NTOP * NMAT * NDEST * NGRID
         MEMUSED = ISEEDPTR + MAX( 4*NNODES, I )
*
         CALL MAKEGRIDS( MEM(CTXTPTR), OUTNUM, NGRID, MEM(PPTR),
     $                   MEM(QPTR) )
*
*        Call individual tests as appropriate.
*
         DO 60 I = 1, NPREC
            DO 50 J = 0, NOP-1
               IF( PREC(I) .EQ. 'I' ) THEN
                  WORKPTR = SAFEINDEX(MEMUSED, ISIZE, ISIZE)
                  WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ISIZE
                  IF( CMEM(OPPTR+J) .EQ. '+' ) THEN
                     CALL ISUMTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR), NDEST,
     $                             MEM(RDESTPTR), MEM(CDESTPTR), NGRID,
     $                             MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                             MEM(ISEEDPTR), MEM(WORKPTR),
     $                             WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '>' ) THEN
                     CALL IAMXTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '<' ) THEN
                     CALL IAMNTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  END IF
               ELSE IF( PREC(I) .EQ. 'S' ) THEN
                  WORKPTR = SAFEINDEX(MEMUSED, ISIZE, SSIZE)
                  WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / SSIZE
                  IF( CMEM(OPPTR+J) .EQ. '+' ) THEN
                     CALL SSUMTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR), NDEST,
     $                             MEM(RDESTPTR), MEM(CDESTPTR), NGRID,
     $                             MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                             MEM(ISEEDPTR), MEM(WORKPTR),
     $                             WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '>' ) THEN
                     CALL SAMXTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '<' ) THEN
                     CALL SAMNTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  END IF
               ELSE IF( PREC(I) .EQ. 'C' ) THEN
                  WORKPTR = SAFEINDEX(MEMUSED, ISIZE, CSIZE)
                  WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / CSIZE
                  IF( CMEM(OPPTR+J) .EQ. '+' ) THEN
                     CALL CSUMTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR), NDEST,
     $                             MEM(RDESTPTR), MEM(CDESTPTR), NGRID,
     $                             MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                             MEM(ISEEDPTR), MEM(WORKPTR),
     $                             WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '>' ) THEN
                     CALL CAMXTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '<' ) THEN
                     CALL CAMNTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  END IF
               ELSE IF( PREC(I) .EQ. 'Z' ) THEN
                  WORKPTR = SAFEINDEX(MEMUSED, ISIZE, ZSIZE)
                  WORKLEN = ( DSIZE * (MEMLEN - WORKPTR + 1) ) / ZSIZE
                  IF( CMEM(OPPTR+J) .EQ. '+' ) THEN
                     CALL ZSUMTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR), NDEST,
     $                             MEM(RDESTPTR), MEM(CDESTPTR), NGRID,
     $                             MEM(CTXTPTR), MEM(PPTR), MEM(QPTR),
     $                             MEM(ISEEDPTR), MEM(WORKPTR),
     $                             WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '>' ) THEN
                     CALL ZAMXTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  ELSE IF( CMEM(OPPTR+J) .EQ. '<' ) THEN
                     CALL ZAMNTEST(OUTNUM, VERB, TREP, TCOH, NSCOPE,
     $                             CMEM(SCOPEPTR), NTOP, CMEM(TOPPTR),
     $                             NMAT, MEM(MPTR), MEM(NPTR),
     $                             MEM(LDSPTR), MEM(LDDPTR),
     $                             MEM(LDIPTR), NDEST, MEM(RDESTPTR),
     $                             MEM(CDESTPTR), NGRID, MEM(CTXTPTR),
     $                             MEM(PPTR), MEM(QPTR), MEM(ISEEDPTR),
     $                             MEM(RAPTR), MEM(CAPTR), K,
     $                             MEM(WORKPTR), WORKLEN)
                  END IF
               END IF
   50       CONTINUE
   60    CONTINUE
         CALL FREEGRIDS( NGRID, MEM(CTXTPTR) )
      END IF
*
      IF( TESTAUX ) THEN
         CALL AUXTEST( OUTNUM, MEM, MEMLEN )
      END IF
*
 1000 FORMAT('AUXILIARY ERROR - IAM MISMATCH: BLACS_PINFO RETURNED',I4,
     $       /,' BLACS_SETUP RETURNED',I4,'.')
 1500 FORMAT('AUXILIARY ERROR - NPROC MISMATCH: BLACS_PINFO RETURNED',
     $       I4,/,' BLACS_SETUP RETURNED',I4,', TESTER THINKS',I4,'.')
 2000 FORMAT('BEGINNING BLACS TESTING, BLACS DEBUG LEVEL =',I2)
 3000 FORMAT('==============================================')
      RETURN
*
*     End of RUNTESTS
*
      END
*
      SUBROUTINE MAKEGRIDS( CONTEXTS, OUTNUM, NGRIDS, P, Q )
      INTEGER NGRIDS, OUTNUM
      INTEGER CONTEXTS(NGRIDS), P(NGRIDS), Q(NGRIDS)
      INTEGER  IBTMYPROC
      EXTERNAL IBTMYPROC
      INTEGER NPROW, NPCOL, MYROW, MYCOL, I
*
      DO 10 I = 1, NGRIDS
         CALL BLACS_GET( 0, 0, CONTEXTS(I) )
         CALL BLACS_GRIDINIT( CONTEXTS(I), 'r', P(I), Q(I) )
   10 CONTINUE
*
      DO 20 I = 1, NGRIDS
         CALL BLACS_GRIDINFO( CONTEXTS(I), NPROW, NPCOL, MYROW, MYCOL )
         IF( NPROW .GT. 0 ) THEN
            IF( NPROW.NE.P(I) .OR. NPCOL.NE.Q(I) ) THEN
               IF( IBTMYPROC() .NE. 0 ) OUTNUM = 6
               WRITE(OUTNUM,1000) I
               IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
               CALL BLACS_ABORT( CONTEXTS(I), -1 )
            END IF
         END IF
   20 CONTINUE
*
 1000 FORMAT('Grid creation error trying to create grid #',I3)
      RETURN
      END
*
      SUBROUTINE FREEGRIDS( NGRIDS, CONTEXTS )
      INTEGER NGRIDS
      INTEGER CONTEXTS(NGRIDS)
      INTEGER I, NPROW, NPCOL, MYROW, MYCOL
*
      DO 10 I = 1, NGRIDS
         CALL BLACS_GRIDINFO( CONTEXTS(I), NPROW, NPCOL, MYROW, MYCOL )
         IF( MYROW.LT.NPROW .AND. MYCOL.LT.NPCOL )
     $      CALL BLACS_GRIDEXIT( CONTEXTS(I) )
   10 CONTINUE
      RETURN
      END
*
      SUBROUTINE AUXTEST( OUTNUM, MEM, MEMLEN )
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, MEMLEN
*     ..
*     .. Array Arguments ..
      INTEGER MEM(MEMLEN)
*     ..
*     .. External Functions ..
      LOGICAL  ALLPASS
      INTEGER  IBTMYPROC, IBTMSGID, BLACS_PNUM
      DOUBLE PRECISION DWALLTIME00
      EXTERNAL ALLPASS, IBTMYPROC, IBTMSGID, BLACS_PNUM
      EXTERNAL DWALLTIME00
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_PINFO, BLACS_GRIDINIT, BLACS_GRIDMAP
      EXTERNAL BLACS_FREEBUFF, BLACS_GRIDEXIT, BLACS_ABORT
      EXTERNAL BLACS_GRIDINFO, BLACS_PCOORD, BLACS_BARRIER
      EXTERNAL BLACS_SET
*     ..
*     .. Local Scalars ..
      LOGICAL AUXPASSED, PASSED, IPRINT
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, CTXT, CTXT2, LDA
      INTEGER I, J, K
      DOUBLE PRECISION DTIME, DEPS
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION START(2), STST(2), KEEP(2)
*     ..
*     .. Executable Statements ..
*
      IPRINT = ( IBTMYPROC() .EQ. 0 )
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) '  '
         WRITE(OUTNUM,1000)
         WRITE(OUTNUM,*) '  '
      END IF
      CALL BLACS_PINFO( I, NPROCS )
      IF( NPROCS .LT. 2 ) THEN
         IF( IPRINT )
     $      WRITE(OUTNUM,*) 'NOT ENOUGH PROCESSES TO PERFORM AUXTESTS'
         RETURN
      END IF
*
*     Make sure BLACS_PNUM and BLACS_PCOORD are inverses of each other
*
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) ' '
         WRITE(OUTNUM,*) 'RUNNING BLACS_PNUM/BLACS_PCOORD TEST'
      END IF
      PASSED = .TRUE.
      NPROCS = NPROCS - MOD(NPROCS,2)
      CALL BLACS_GET( 0, 0, CTXT )
      CALL BLACS_GRIDINIT( CTXT, 'r', 1, NPROCS )
      CALL BLACS_GRIDINFO( CTXT, NPROW, NPCOL, MYROW, MYCOL )
      IF( MYROW.GE.NPROW .OR. MYCOL.GE.NPCOL ) GOTO 100
      DO 10 I = 1, NPROCS
         K = BLACS_PNUM( CTXT, 0, I-1 )
         CALL BLACS_PCOORD( CTXT,  BLACS_PNUM( CTXT, 0, I-1 ), J, K )
         IF( PASSED ) PASSED = ( J.EQ.0 .AND. K.EQ.I-1 )
   10 CONTINUE
      K = 1
      IF( PASSED ) K = 0
      CALL IGSUM2D( CTXT, 'a', ' ', 1, 1, K, 1, -1, 0 )
      PASSED = ( K .EQ. 0 )
      AUXPASSED = PASSED
      IF( IPRINT ) THEN
         IF( PASSED ) THEN
            WRITE(OUTNUM,*) 'PASSED  BLACS_PNUM/BLACS_PCOORD TEST'
         ELSE
            WRITE(OUTNUM,*) 'FAILED  BLACS_PNUM/BLACS_PCOORD TEST'
         END IF
         WRITE(OUTNUM,*) '  '
      END IF
*
*     Test to see if DGSUM2D is repeatable when repeatability flag is set
*     Skip test if DGSUM2D is repeatable when repeatability flag is not set
*     NOTE: do not change the EPS calculation loop; it is figured in this
*           strange way so that it ports across platforms
*
      IF( IPRINT ) WRITE(OUTNUM,*) 'RUNNING REPEATABLE SUM TEST'
      J = 0
   12 CONTINUE
      PASSED = .TRUE.
      START(1) = 1.0D0
   15 CONTINUE
         DEPS = START(1)
         START(1) = START(1) / 2.0D0
         STST(1) = 1.0D0 + START(1)
      IF (STST(1) .NE. 1.0D0) GOTO 15
*
      START(1) = DEPS / DBLE(NPCOL-1)
      IF (MYCOL .EQ. 3) START(1) = 1.0D0
      START(2) = 7.00005D0 * NPCOL
      STST(1) = START(1)
      STST(2) = START(2)
      CALL BLACS_SET(CTXT, 15, J)
      CALL DGSUM2D(CTXT, 'a', 'f', 2, 1, STST, 2, -1, 0)
      KEEP(1) = STST(1)
      KEEP(2) = STST(2)
      DO 30 I = 1, 3
*
*        Have a different guy waste time so he enters combine last
*
         IF (MYCOL .EQ. I) THEN
             DTIME = DWALLTIME00()
   20        CONTINUE
             IF (DWALLTIME00() - DTIME .LT. 2.0D0) GOTO 20
         END IF
         STST(1) = START(1)
         STST(2) = START(2)
         CALL DGSUM2D(CTXT, 'a', 'f', 2, 1, STST, 2, -1, 0)
         IF ( (KEEP(1).NE.STST(1)) .OR. (KEEP(2).NE.STST(2)) )
     $      PASSED = .FALSE.
   30 CONTINUE
      K = 1
      IF (PASSED) K = 0
      CALL IGSUM2D( CTXT, 'a', ' ', 1, 1, K, 1, -1, 0 )
      PASSED = (K .EQ. 0)
      IF (J .EQ. 0) THEN
         IF (.NOT.PASSED) THEN
            J = 1
            GOTO 12
         ELSE IF( IPRINT ) THEN
            WRITE(OUTNUM,*) 'SKIPPED REPEATABLE SUM TEST'
            WRITE(OUTNUM,*) ' '
         END IF
      END IF
*
      IF (J .EQ. 1) THEN
         AUXPASSED = AUXPASSED .AND. PASSED
         IF( IPRINT ) THEN
            IF( PASSED ) THEN
               WRITE(OUTNUM,*) 'PASSED  REPEATABLE SUM TEST'
            ELSE
               WRITE(OUTNUM,*) 'FAILED  REPEATABLE SUM TEST'
            END IF
            WRITE(OUTNUM,*) ' '
         END IF
      END IF
*
*     Test BLACS_GRIDMAP: force a column major ordering, starting at an
*     arbitrary processor
*
      PASSED = .TRUE.
      IF( IPRINT ) WRITE(OUTNUM,*) 'RUNNING BLACS_GRIDMAP TEST'
      NPROW = 2
      NPCOL = NPROCS / NPROW
      DO 40 I = 0, NPROCS-1
         MEM(I+1) = BLACS_PNUM( CTXT, 0, MOD(I+NPCOL, NPROCS) )
   40 CONTINUE
      CALL BLACS_GET( CTXT, 10, CTXT2 )
      CALL BLACS_GRIDMAP( CTXT2, MEM, NPROW, NPROW, NPCOL )
      CALL BLACS_GRIDINFO( CTXT2, NPROW, NPCOL, MYROW, MYCOL )
      PASSED = ( NPROW.EQ.2 .AND. NPCOL.EQ.NPROCS/2 )
*
*     Fan in pids for final check: Note we assume SD/RV working
*
      IF( PASSED ) THEN
         K = BLACS_PNUM( CTXT2, MYROW, MYCOL )
         IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
            DO 60 J = 0, NPCOL-1
               DO 50 I = 0, NPROW-1
                  IF( I.NE.0 .OR. J.NE.0 )
     $               CALL IGERV2D( CTXT2, 1, 1, K, 1, I, J )
                  IF ( PASSED )
     $               PASSED = ( K .EQ. BLACS_PNUM(CTXT2, I, J) )
   50          CONTINUE
   60       CONTINUE
         ELSE
            CALL IGESD2D( CTXT2, 1, 1, K, 1, 0, 0 )
         END IF
      END IF
      K = 1
      IF ( PASSED ) K = 0
      CALL IGSUM2D( CTXT, 'a', ' ', 1, 1, K, 1, -1, 0 )
      PASSED = ( K .EQ. 0 )
      AUXPASSED = AUXPASSED .AND. PASSED
      IF( IPRINT ) THEN
         IF( PASSED ) THEN
            WRITE(OUTNUM,*) 'PASSED  BLACS_GRIDMAP TEST'
         ELSE
            WRITE(OUTNUM,*) 'FAILED  BLACS_GRIDMAP TEST'
         END IF
         WRITE(OUTNUM,*) ' '
      END IF
*
      IF( IPRINT ) WRITE(OUTNUM,*) 'CALL BLACS_FREEBUFF'
      CALL BLACS_FREEBUFF( CTXT, 0 )
      CALL BLACS_FREEBUFF( CTXT, 1 )
      J = 0
      CALL IGSUM2D( CTXT2, 'All', ' ', 1, 1, J, 1, -1, MYCOL )
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) 'DONE BLACS_FREEBUFF'
         WRITE(OUTNUM,*) ' '
      END IF
*
*     Make sure barriers don't interfere with each other
*
      IF( IPRINT ) WRITE(OUTNUM,*) 'CALL BARRIER'
      CALL BLACS_BARRIER(CTXT2, 'A')
      CALL BLACS_BARRIER(CTXT2, 'R')
      CALL BLACS_BARRIER(CTXT2, 'C')
      CALL BLACS_BARRIER(CTXT2, 'R')
      CALL BLACS_BARRIER(CTXT2, 'A')
      CALL BLACS_BARRIER(CTXT2, 'C')
      CALL BLACS_BARRIER(CTXT2, 'C')
      CALL BLACS_BARRIER(CTXT2, 'R')
      CALL BLACS_BARRIER(CTXT2, 'A')
      J = 0
      CALL IGSUM2D( CTXT2, 'All', ' ', 1, 1, J, 1, -1, MYCOL )
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) 'DONE BARRIER'
         WRITE(OUTNUM,*) ' '
      END IF
*
*     Ensure contiguous sends are locally-blocking
*
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) 'The following tests will hang if your BLACS'//
     $                   ' are not locally blocking:'
         WRITE(OUTNUM,*) 'RUNNING LOCALLY-BLOCKING CONTIGUOUS SEND TEST'
      END IF
      K = MIN( MEMLEN, 50000 )
*
*     Initialize send buffer
*
      DO 70 J = 1, K
         MEM(J) = 1
   70 CONTINUE
*
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL IGESD2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
         CALL IGESD2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
         CALL IGESD2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, NPROW-1, NPCOL-1 )
      ELSE IF( MYROW.EQ.NPROW-1 .AND. MYCOL.EQ.NPCOL-1 ) THEN
         CALL IGESD2D( CTXT2, K, 1, MEM, K, 0, 0 )
         CALL IGESD2D( CTXT2, K, 1, MEM, K, 0, 0 )
         CALL IGESD2D( CTXT2, K, 1, MEM, K, 0, 0 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, 0, 0 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, 0, 0 )
         CALL IGERV2D( CTXT2, K, 1, MEM, K, 0, 0 )
      END IF
      J = 0
      CALL IGSUM2D( CTXT2, 'All', ' ', 1, 1, J, 1, -1, MYCOL )
      IF( IPRINT )
     $   WRITE(OUTNUM,*) 'PASSED  LOCALLY-BLOCKING CONTIGUOUS SEND TEST'
*
*     Ensure non-contiguous sends are locally-blocking
*
      J = 4
      LDA = K / J
      I = MAX( 2, LDA / 4 )
      IF( IPRINT )
     $   WRITE(OUTNUM,*) 'RUNNING LOCALLY-BLOCKING NON-CONTIGUOUS '//
     $                   'SEND TEST'
      IF( MYROW.EQ.0 .AND. MYCOL.EQ.0 ) THEN
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, NPROW-1, NPCOL-1 )
      ELSE IF( MYROW.EQ.NPROW-1 .AND. MYCOL.EQ.NPCOL-1 ) THEN
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, 0, 0 )
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, 0, 0 )
         CALL IGESD2D( CTXT2, I, J, MEM, LDA, 0, 0 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, 0, 0 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, 0, 0 )
         CALL IGERV2D( CTXT2, I, J, MEM, LDA, 0, 0 )
      END IF
      CALL IGSUM2D( CTXT2, 'All', ' ', 1, 1, J, 1, -1, MYCOL )
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*)'PASSED  LOCALLY-BLOCKING NON-CONTIGUOUS '//
     $                  'SEND TEST'
         WRITE(OUTNUM,*) '  '
      END IF
*
*     Note that we already tested the message ID setting/getting in
*     first call to IBTMSGID()
*
      IF( IPRINT ) WRITE(OUTNUM,*) 'RUNNING BLACS_SET/BLACS_GET TESTS'
      J = 0
      CALL BLACS_SET( CTXT2, 11, 3 )
      CALL BLACS_SET( CTXT2, 12, 2 )
      CALL BLACS_GET( CTXT2, 12, I )
      CALL BLACS_GET( CTXT2, 11, K )
      IF( K.NE.3 ) J = J + 1
      IF( I.NE.2 ) J = J + 1
      CALL BLACS_SET( CTXT2, 13, 3 )
      CALL BLACS_SET( CTXT2, 14, 2 )
      CALL BLACS_GET( CTXT2, 14, I )
      CALL BLACS_GET( CTXT2, 13, K )
      IF( K.NE.3 ) J = J + 1
      IF( I.NE.2 ) J = J + 1
*
*     See if anyone had error, and print result
*
      CALL IGSUM2D( CTXT2, 'All', ' ', 1, 1, J, 1, -1, MYCOL )
      PASSED = (J .EQ. 0)
      AUXPASSED = AUXPASSED .AND. PASSED
      IF( IPRINT ) THEN
         IF( PASSED ) THEN
            WRITE(OUTNUM,*) 'PASSED  BLACS_SET/BLACS_GET TESTS'
         ELSE
            WRITE(OUTNUM,*) 'FAILED  BLACS_SET/BLACS_GET TESTS'
         END IF
         WRITE(OUTNUM,*) ' '
      END IF
*
      IF( IPRINT ) WRITE(OUTNUM,*) 'CALL BLACS_GRIDEXIT'
      CALL BLACS_GRIDEXIT(CTXT)
      CALL BLACS_GRIDEXIT(CTXT2)
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) 'DONE BLACS_GRIDEXIT'
         WRITE(OUTNUM,*) '  '
      END IF
*
  100 CONTINUE
*
      PASSED = ALLPASS(AUXPASSED)
      IF( IPRINT ) THEN
         WRITE(OUTNUM,*) 'The final auxiliary test is for BLACS_ABORT.'
         WRITE(OUTNUM,*) 'Immediately after this message, all '//
     $                   'processes should be killed.'
         WRITE(OUTNUM,*) 'If processes survive the call, your BLACS_'//
     $                   'ABORT is incorrect.'
      END IF
      CALL BLACS_PINFO( I, NPROCS )
      CALL BLACS_GET( 0, 0, CTXT )
      CALL BLACS_GRIDINIT( CTXT, 'r', 1, NPROCS )
      CALL BLACS_BARRIER(CTXT, 'A')
      CALL BLACS_GRIDINFO( CTXT, NPROW, NPCOL, MYROW, MYCOL )
*
*     Test BLACS_ABORT
*
      IF( MYROW.EQ.NPROW/2 .AND. MYCOL.EQ.NPCOL/2 ) THEN
         CALL BLACS_ABORT( CTXT, -1 )
*
*     Other procs try to cause a hang: should be killed by BLACS_ABORT
*
      ELSE
         I = 1
110      CONTINUE
            I = I + 3
            I = I - 2
            I = I - 1
         IF( I.EQ.1 ) GOTO 110
      end if
*
 1000 FORMAT('AUXILIARY TESTS: BEGIN.')
      RETURN
      END
*
      SUBROUTINE BTTRANSCHAR(TRANSTO, N, CMEM, IMEM)
      CHARACTER TRANSTO
      INTEGER N
      CHARACTER*1 CMEM(N)
      INTEGER IMEM(N)
      INTEGER I
*
      IF( TRANSTO .EQ. 'I' ) THEN
         DO 10 I = 1, N
            IMEM(I) = ICHAR( CMEM(I) )
   10    CONTINUE
      ELSE
         DO 20 I = 1, N
            CMEM(I) = CHAR( IMEM(I) )
   20    CONTINUE
      END IF
      RETURN
      END
*
      SUBROUTINE BTINFO( TEST, MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM,
     $                   CMEMLEN, OUTNUM, NOP, NSCOPE, TREP, TCOH, NTOP,
     $                   NSHAPE, NMAT, NSRC, NGRID, OPPTR, SCOPEPTR,
     $                   TOPPTR, UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR,
     $                   LDDPTR, LDIPTR, RSRCPTR, CSRCPTR, RDESTPTR,
     $                   CDESTPTR, PPTR, QPTR )
*
*     .. Scalar Arguments ..
      CHARACTER*1 TEST
      INTEGER CDESTPTR, CMEMLEN, CMEMUSED, CSRCPTR, DIAGPTR, LDDPTR,
     $        LDIPTR, LDSPTR, MEMLEN, MEMUSED, MPTR, NGRID, NMAT, NOP,
     $        NPTR, NSCOPE, NSHAPE, NSRC, NTOP, OPPTR, OUTNUM, PPTR,
     $        QPTR, RDESTPTR, RSRCPTR, SCOPEPTR, TCOH, TOPPTR, TREP,
     $        UPLOPTR
*     ..
*     .. Array Arguments ..
      CHARACTER*1 CMEM(CMEMLEN)
      INTEGER MEM(MEMLEN)
*     ..
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTMSGID, IBTSIZEOF
      EXTERNAL IBTMYPROC, IBTMSGID, IBTSIZEOF
*     ..
*     .. Local Scalars ..
      INTEGER IAM, ISIZE, DSIZE
*     ..
*     .. Local Arrays ..
      INTEGER ITMP(2)
*     ..
*     .. Executable Statements ..
*
      IAM = IBTMYPROC()
      IF( IAM .EQ. 0 ) THEN
         IF( TEST .EQ. 'S' ) THEN
            CALL RDSDRV( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
         ELSE IF( TEST .EQ. 'B' ) THEN
            CALL RDBSBR( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
         ELSE
            CALL RDCOMB( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
         END IF
         ITMP(1) = MEMUSED
         ITMP(2) = CMEMUSED
         CALL BTSEND( 3, 2, ITMP, -1, IBTMSGID()+3 )
         IF( MEMLEN .GE. MEMUSED + CMEMUSED ) THEN
            CALL BTTRANSCHAR( 'I', CMEMUSED, CMEM, MEM(MEMUSED+1) )
         ELSE
            ISIZE = IBTSIZEOF('I')
            DSIZE = IBTSIZEOF('D')
            WRITE(OUTNUM,1000) ( (MEMUSED+CMEMUSED)*ISIZE + DSIZE-1 )
     $                         / DSIZE
            CALL BLACS_ABORT(-1, -1)
         END IF
         CALL BTSEND( 3, MEMUSED+CMEMUSED, MEM, -1, IBTMSGID()+4 )
      ELSE
         CALL BTRECV( 3, 2, ITMP, 0, IBTMSGID()+3 )
         MEMUSED = ITMP(1)
         CMEMUSED = ITMP(2)
         IF( MEMLEN .GE. MEMUSED + CMEMUSED ) THEN
            CALL BTRECV( 3, MEMUSED+CMEMUSED, MEM, 0, IBTMSGID()+4 )
            CALL BTTRANSCHAR( 'C', CMEMUSED, CMEM, MEM(MEMUSED+1) )
         ELSE
            ISIZE = IBTSIZEOF('I')
            DSIZE = IBTSIZEOF('D')
            WRITE(OUTNUM,1000) ( (MEMUSED+CMEMUSED)*ISIZE + DSIZE-1 )
     $                         / DSIZE
            CALL BLACS_ABORT(-1, -1)
         END IF
      END IF
      CALL BTUNPACK( TEST, MEM, MEMUSED, NOP, NSCOPE, TREP, TCOH, NTOP,
     $               NSHAPE, NMAT, NSRC, NGRID, OPPTR, SCOPEPTR, TOPPTR,
     $               UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR, LDDPTR,
     $               LDIPTR, RSRCPTR, CSRCPTR, RDESTPTR, CDESTPTR, PPTR,
     $               QPTR)
*
 1000 FORMAT('MEM array too short to pack CMEM; increase to at least',
     $       I7)
*
      RETURN
*
*     End BTINFO
*
      END
*
      SUBROUTINE RDBTIN( TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX, NPREC,
     $                   PREC, VERB, OUTNUM )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL TESTSDRV, TESTBSBR, TESTCOMB, TESTAUX
      INTEGER NPREC, OUTNUM, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 PREC(*)
*     ..
*
*  Purpose
*  =======
*  RDBTIN:  Read and process the top-level input file BT.dat.
*
*  Arguments
*  =========
*  TESTSDRV (output) LOGICAL
*           Run any point-to-point tests?
*
*  TESTBSBR (output) LOGICAL
*           Run any broadcast tests?
*
*  TESTCOMB (output) LOGICAL
*           Run any combine-operation tests (e.g. MAX)
*
*  TESTAUX  (output) LOGICAL
*           Run any auxiliary tests?
*
*  NPREC    (output) INTEGER
*           Number of different precisions to test. (up to 5, as determined
*           by the parameter PRECMAX down in the code.)
*
*  PREC     (output) CHARACTER*1 array, dimension 5
*           Prefix letter of each precision to test, from the set
*           {'C', 'D', 'I', 'S', 'Z'}
*
*  VERB     (output) INTEGER
*           Output verbosity for this test run.
*            0 = Print only "BEGIN [SDRV/BSBR/COMB]", followed by PASSED
*                or FAILED message
*            1 = Same as 0, but also prints out header explaining all tests
*                to be run.
*            2 = Prints out info before and after every individual test.
*
*  OUTNUM   (output) INTEGER
*           Unit number for output file.
*  ======================================================================
*
*
*     .. Parameters ..
      INTEGER PRECMAX, VERBMAX, IN
      PARAMETER ( PRECMAX = 5, VERBMAX = 2, IN = 11 )
*     ..
*     .. Local Scalars ..
      INTEGER I
      CHARACTER*1 CH
      LOGICAL READERROR
*     ..
*     .. Local Arrays ..
      CHARACTER*80 HEADER, OUTNAME
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. Executable Statements
*
*     Open and read the file blacstest.dat.  Expected format is
*     -----
*     'One line of free text intended as a comment for each test run'
*     integer             Unit number of output file
*     string              Name of output file (ignored if unit = 6)
*     {'T'|'F'}           Run any point to point tests?
*     {'T'|'F'}           Run any broadcast tests?
*     {'T'|'F'}           Run any combine-operator tests?
*     {'T'|'F'}           Run the auxiliary tests?
*     integer             Number of precisions to test - up to 99
*     array of CHAR*1's   Specific precisions to test
*     integer             Output verb (1-n, n=most verbose)
*     integer             Number of nodes required by largest test case
*     -----
*     Note that the comments to the right of each line are present
*     in the sample blacstest.dat file included with this
*     distribution, but they are not required.
*
*     The array of CHAR*1's is expected to have length equal to the
*     integer in the previous line - if it is shorter, problems may
*     occur later; if it is longer, the trailing elements will just
*     be ignored.  The verb is expected to be an integer
*     between 1 and n inclusive and will be set to 1 if outside
*     this range.
*
*     Only process 0 should be calling this routine
*
      READERROR = .FALSE.
      OPEN( UNIT = IN, FILE = 'bt.dat', STATUS = 'OLD' )
      READ(IN, *) HEADER
      READ(IN, *) OUTNUM
      READ(IN, *) OUTNAME
*
*     Open and prepare output file
*
      IF( OUTNUM.NE.6 .AND. OUTNUM.NE.0 )
     $  OPEN( UNIT = OUTNUM, FILE = OUTNAME, STATUS = 'UNKNOWN' )
      WRITE(OUTNUM, *) HEADER
*
*     Determine which tests to run
*
      READ(IN, *) CH
      IF( LSAME(CH, 'T') ) THEN
         TESTSDRV = .TRUE.
      ELSE IF( LSAME(CH, 'F') ) THEN
         TESTSDRV = .FALSE.
      ELSE
         WRITE(OUTNUM, 1000) 'SDRV', CH
         READERROR = .TRUE.
      END IF
*
      READ(IN, *) CH
      IF( LSAME(CH, 'T') ) THEN
         TESTBSBR = .TRUE.
      ELSE IF(LSAME( CH, 'F') ) THEN
         TESTBSBR = .FALSE.
      ELSE
         WRITE(OUTNUM, 1000) 'BSBR', CH
         READERROR = .TRUE.
      END IF
*
      READ(IN, *) CH
      IF( LSAME(CH, 'T') ) THEN
         TESTCOMB = .TRUE.
      ELSE IF( LSAME(CH, 'F') ) THEN
         TESTCOMB = .FALSE.
      ELSE
         WRITE(OUTNUM, 1000) 'COMB', CH
         READERROR = .TRUE.
      END IF
*
      READ(IN, *) CH
      IF( LSAME(CH, 'T') ) THEN
         TESTAUX = .TRUE.
      ELSE IF( LSAME(CH, 'F') ) THEN
         TESTAUX = .FALSE.
      ELSE
         WRITE(OUTNUM, 1000) 'AUX ', CH
         READERROR = .TRUE.
      END IF
*
*     Get # of precisions, and precisions to test
*
      READ(IN, *) NPREC
      IF( NPREC .LT. 0 ) THEN
         NPREC = 0
      ELSE IF( NPREC. GT. PRECMAX ) THEN
         WRITE(OUTNUM, 2000) NPREC, PRECMAX, PRECMAX
         NPREC = PRECMAX
      END IF
*
      READ(IN, *) ( PREC(I), I = 1, NPREC )
      DO 100 I = 1, NPREC
         IF( LSAME(PREC(I), 'C') ) THEN
            PREC(I) = 'C'
         ELSE IF( LSAME(PREC(I), 'D') ) THEN
            PREC(I) = 'D'
         ELSE IF( LSAME(PREC(I), 'I') ) THEN
            PREC(I) = 'I'
         ELSE IF( LSAME(PREC(I), 'S') ) THEN
            PREC(I) = 'S'
         ELSE IF( LSAME(PREC(I), 'Z') ) THEN
            PREC(I) = 'Z'
         ELSE
            WRITE(OUTNUM, 3000) PREC(I)
            READERROR = .TRUE.
         END IF
  100 CONTINUE
*
      READ(IN, *) VERB
*
      IF( VERB .GT. VERBMAX ) THEN
         WRITE(OUTNUM, 4000) VERB, VERBMAX, VERBMAX
         VERB = VERBMAX
      ELSE IF( VERB .LT. 0 ) THEN
         WRITE(OUTNUM, 5000) VERB
         VERB = 0
      END IF
*
*     Abort if there was a fatal error
*
      IF( READERROR ) THEN
         WRITE(OUTNUM, 6000)
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE( OUTNUM )
         STOP
      END IF
*
 1000 FORMAT( 'INVALID CHARACTER FOR ',A4,' TESTS ''', A1,
     $        ''' (EXPECTED T/F)' )
 2000 FORMAT( 'NUMBER OF PRECISIONS ', I6, ' GREATER THAN ', I6,
     $        ' - SETTING TO ', I6, '.')
 3000 FORMAT( 'UNRECOGNIZABLE PRECISION ENTRY ''', A1,
     $        ''' - EXPECTED ''C'', ''D'', ''I'', ''S'', OR ''Z''.')
 4000 FORMAT( 'VERBOSITY ', I4, ' GREATER THAN ',I4,
     $        ' - SETTING TO ',I4,'.')
 5000 FORMAT( 'VERBOSITY ', I4, ' LESS THAN 0 - SETTING TO 0' )
 6000 FORMAT( 'FATAL INPUT FILE ERROR - ABORTING RUN.' )
*
      RETURN
*
*     End of RDBTIN
*
      END
*
      INTEGER FUNCTION IBTMSGID()
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     PURPOSE
*     =======
*     IBTMSGID : returns a ID for tester communication.
*
      INTEGER MINID
      INTEGER ITMP(2)
      SAVE MINID
      DATA MINID /-1/
*
*     On first call, reserve 1st 1000 IDs for tester use
*
      IF (MINID .EQ. -1) THEN
         CALL BLACS_GET( -1, 1, ITMP )
         MINID = ITMP(1)
         ITMP(1) = ITMP(1) + 1000
         CALL BLACS_SET( -1, 1, ITMP )
      END IF
*
*     return the minimum allowable ID
*
      IBTMSGID = MINID
*
      RETURN
      END
*
      SUBROUTINE BTUNPACK(TEST, MEM, MEMLEN, NOP, NSCOPE, TREP, TCOH,
     $                    NTOP, NSHAPE, NMAT, NSRC, NGRID, OPPTR,
     $                    SCOPEPTR, TOPPTR, UPLOPTR, DIAGPTR, MPTR,
     $                    NPTR, LDSPTR, LDDPTR, LDIPTR, RSRCPTR,
     $                    CSRCPTR, RDESTPTR, CDESTPTR, PPTR, QPTR)
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 TEST
      INTEGER CDESTPTR, CSRCPTR, DIAGPTR, LDDPTR, LDIPTR, LDSPTR,
     $        MEMLEN, MPTR, NGRID, NMAT, NOP, NPTR, NSCOPE, NSHAPE,
     $        NSRC, NTOP, OPPTR, PPTR, QPTR, RDESTPTR, RSRCPTR,
     $        SCOPEPTR, TCOH, TOPPTR, TREP, UPLOPTR
*     ..
*     .. Array Arguments ..
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  BTUNPACK: Figure pointers into MEM where the various input values
*  are stored.
*
*  Arguments
*  =========
*  TEST     (input) CHARACTER*1
*           The test we're unpacking for:
*            = 'S' : SDRV test
*            = 'B' : BSBR test
*            = 'C' : Combine test
*
*  MEM      (input) INTEGER array of dimension MEMLEN
*           Memory containing values and number of items.
*
*  MEMLEN   (input/output) INTEGER
*           The number of elements that are used in MEM.
*
*  .
*  .
*  .
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER NDEST, NLDI
*     ..
*     .. Executable Statements ..
*
*     Test is SDRV
*
      IF( TEST .EQ. 'S' ) THEN
         NOP    = 0
         NSHAPE = MEM(MEMLEN-3)
         NSCOPE = 0
         TREP   = 0
         TCOH   = 0
         NTOP   = 0
         NMAT   = MEM(MEMLEN-2)
         NLDI   = 0
         NSRC   = MEM(MEMLEN-1)
         NDEST  = NSRC
         NGRID  = MEM(MEMLEN)
         MEMLEN = MEMLEN - 3
*
*     Test is BSBR
*
      ELSE IF ( TEST .EQ. 'B' ) THEN
         NOP    = 0
         NSCOPE = MEM(MEMLEN-5)
         TREP   = 0
         TCOH   = 0
         NTOP   = MEM(MEMLEN-4)
         NSHAPE = MEM(MEMLEN-3)
         NMAT   = MEM(MEMLEN-2)
         NLDI   = 0
         NSRC   = MEM(MEMLEN-1)
         NDEST  = 0
         NGRID  = MEM(MEMLEN)
         MEMLEN = MEMLEN - 5
*
*     Test is COMB
*
      ELSE
         NOP    = MEM(MEMLEN-7)
         NSCOPE = MEM(MEMLEN-6)
         TREP   = MEM(MEMLEN-5)
         TCOH   = MEM(MEMLEN-4)
         NTOP   = MEM(MEMLEN-3)
         NSHAPE = 0
         NMAT   = MEM(MEMLEN-2)
         NLDI   = NMAT
         NSRC   = 0
         NDEST  = MEM(MEMLEN-1)
         NGRID  = MEM(MEMLEN)
         MEMLEN = MEMLEN - 6
      END IF
      OPPTR = 1
      SCOPEPTR = OPPTR + NOP
      TOPPTR = SCOPEPTR + NSCOPE
      UPLOPTR = TOPPTR + NTOP
      DIAGPTR = UPLOPTR + NSHAPE
      MPTR = 1
      NPTR = MPTR + NMAT
      LDSPTR = NPTR + NMAT
      LDDPTR = LDSPTR + NMAT
      LDIPTR = LDDPTR + NMAT
      RSRCPTR = LDIPTR + NLDI
      CSRCPTR = RSRCPTR + NSRC
      RDESTPTR = CSRCPTR + NSRC
      CDESTPTR = RDESTPTR + NDEST
      PPTR = CDESTPTR + NDEST
      QPTR = PPTR + NGRID
      IF( NSRC .EQ. 0 ) NSRC = NDEST
*
      RETURN
*
*     End of BTUNPACK
*
      END
*
      INTEGER FUNCTION SAFEINDEX(INDX, SIZE1, SIZE2)
*
*     .. Scalar Arguments ..
      INTEGER INDX, SIZE1, SIZE2
*     ..
*
*  If you have an array with elements of SIZE1 bytes, of which you
*  have used INDX-1 elements, returns the index necessary to keep it
*  on a SIZE2 boundary (assuming it was SIZE2 aligned in the first place).
*
*     .. Local scalars ..
      INTEGER I
*     ..
*     .. Executable Statements ..
*
*     Take into account that Fortran starts arrays at 1, not 0
*
      I = INDX - 1
   10 CONTINUE
      IF( MOD(I*SIZE1, SIZE2) .EQ. 0 ) GOTO 20
         I = I + 1
      GOTO 10
   20 CONTINUE
*
      SAFEINDEX = I + 1
*
      RETURN
      END
*
*
      SUBROUTINE RDSDRV( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMUSED, MEMLEN, CMEMUSED, CMEMLEN, OUTNUM
*     ..
*     .. Array Arguments ..
      CHARACTER*1 CMEM(CMEMLEN)
      INTEGER MEM(MEMLEN)
*     ..
*
*     Purpose
*     =======
*     RDSDRV:  Read and process the input file SDRV.dat.
*
*     Arguments
*     =========
*     MEMUSED  (output) INTEGER
*              Number of elements in MEM that this subroutine ends up using.
*
*     MEM      (output) INTEGER array of dimension memlen
*              On output, holds information read in from sdrv.dat.
*
*     MEMLEN   (input) INTEGER
*              Number of elements of MEM that this subroutine
*              may safely write into.
*
*     CMEMUSED (output) INTEGER
*              Number of elements in CMEM that this subroutine ends up using.
*
*     CMEM     (output) CHARACTER*1 array of dimension cmemlen
*              On output, holds the values for UPLO and DIAG.
*
*     CMEMLEN  (input) INTEGER
*              Number of elements of CMEM that this subroutine
*              may safely write into.
*
*     OUTNUM   (input) INTEGER
*              Unit number of the output file.
*
*     =================================================================
*
*     .. Parameters ..
      INTEGER SDIN
      PARAMETER( SDIN = 12 )
*     ..
*     .. External Functions ..
      LOGICAL  LSAME
      EXTERNAL LSAME
*     ..
*     .. Local Scalars ..
      INTEGER NSHAPE, NMAT, NSRC, NGRID, I, J
      INTEGER UPLOPTR, DIAGPTR, MPTR, NPTR, LDSPTR, LDDPTR, RSRCPTR
      INTEGER CSRCPTR, RDESTPTR, CDESTPTR, PPTR, QPTR
*     ..
*     .. Executable Statements
*
*     Open and read the file sdrv.dat.  The expected format is
*     below.
*
*------
*integer                         number of shapes of the matrix
*array of CHAR*1's               UPLO
*array of CHAR*1's               DIAG: unit diagonal or not?
*integer                         number of nmat
*array of integers               M: number of rows in matrix
*array of integers               N: number of columns in matrix
*integer                         LDA: leading dimension on source proc
*integer                         LDA: leading dimension on dest proc
*integer                         number of source/dest pairs
*array of integers               RSRC: process row of message source
*array of integers               CSRC: process column of msg. src.
*array of integers               RDEST: process row of msg. dest.
*array of integers               CDEST: process column of msg. dest.
*integer                         Number of grids
*array of integers               NPROW: number of rows in process grid
*array of integers               NPCOL: number of col's in proc. grid
*------
*  note: UPLO stands for 'upper or lower trapezoidal or general
*        rectangular.'
*  note: the text descriptions as shown above are present in
*             the sample sdrv.dat included with this distribution,
*             but are not required.
*
*     Read input file
*
      MEMUSED = 1
      CMEMUSED = 1
      OPEN(UNIT = SDIN, FILE = 'sdrv.dat', STATUS = 'OLD')
*
*     Read in number of shapes, and values of UPLO and DIAG
*
      READ(SDIN, *) NSHAPE
      UPLOPTR = CMEMUSED
      DIAGPTR = UPLOPTR + NSHAPE
      CMEMUSED = DIAGPTR + NSHAPE
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NSHAPE, 'MATRIX SHAPES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSHAPE .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'MATRIX SHAPE.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
*
*     Read in, upcase, and fatal error if UPLO/DIAG not recognized
*
      READ(SDIN, *) ( CMEM(UPLOPTR+I), I = 0, NSHAPE-1 )
      DO 30 I = 0, NSHAPE-1
         IF( LSAME(CMEM(UPLOPTR+I), 'G') ) THEN
            CMEM(UPLOPTR+I) = 'G'
         ELSE IF( LSAME(CMEM(UPLOPTR+I), 'U') ) THEN
            CMEM(UPLOPTR+I) = 'U'
         ELSE IF( LSAME(CMEM(UPLOPTR+I), 'L') ) THEN
            CMEM(UPLOPTR+I) = 'L'
         ELSE
            WRITE(OUTNUM, 3000) 'UPLO ', CMEM(UPLOPTR+I)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   30 CONTINUE
*
      READ(SDIN, *) ( CMEM(DIAGPTR+I), I = 0, NSHAPE-1 )
      DO 40 I = 0, NSHAPE-1
         IF( CMEM(UPLOPTR+I) .NE. 'G' ) THEN
            IF( LSAME(CMEM(DIAGPTR+I), 'U') ) THEN
               CMEM( DIAGPTR+I ) = 'U'
            ELSE IF( LSAME(CMEM(DIAGPTR+I), 'N') ) THEN
               CMEM(DIAGPTR+I) = 'N'
            ELSE
               WRITE(OUTNUM, 3000) 'DIAG ', CMEM(DIAGPTR+I)
               IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
               STOP
            END IF
         END IF
   40 CONTINUE
*
*     Read in number of matrices, and values for M, N, LDASRC, and LDADEST
*
      READ(SDIN, *) NMAT
      MPTR = MEMUSED
      NPTR = MPTR + NMAT
      LDSPTR = NPTR + NMAT
      LDDPTR = LDSPTR + NMAT
      MEMUSED = LDDPTR + NMAT
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'MATRICES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NMAT .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'MATRIX.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM( MPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( NPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDSPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDDPTR+I ), I = 0, NMAT-1 )
*
*     Make sure matrix values are legal
*
      CALL CHKMATDAT( OUTNUM, 'SDRV.dat', .FALSE., NMAT, MEM(MPTR),
     $                MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR), MEM(LDDPTR) )
*
*     Read in number of src/dest pairs, and values of src/dest
*
      READ(SDIN, *) NSRC
      RSRCPTR  = MEMUSED
      CSRCPTR  = RSRCPTR  + NSRC
      RDESTPTR = CSRCPTR  + NSRC
      CDESTPTR = RDESTPTR + NSRC
      MEMUSED  = CDESTPTR + NSRC
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'SRC/DEST.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSRC .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'SRC/DEST.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM(RSRCPTR+I), I = 0, NSRC-1 )
      READ(SDIN, *) ( MEM(CSRCPTR+I), I = 0, NSRC-1 )
      READ(SDIN, *) ( MEM(RDESTPTR+I), I = 0, NSRC-1 )
      READ(SDIN, *) ( MEM(CDESTPTR+I), I = 0, NSRC-1 )
*
*     Read in number of grids pairs, and values of P (process rows) and
*     Q (process columns)
*
      READ(SDIN, *) NGRID
      PPTR = MEMUSED
      QPTR = PPTR + NGRID
      MEMUSED = QPTR + NGRID
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NGRID, 'PROCESS GRIDS.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NGRID .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'PROCESS GRID'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE( OUTNUM )
         STOP
      END IF
*
      READ(SDIN, *) ( MEM(PPTR+I), I = 0, NGRID-1 )
      READ(SDIN, *) ( MEM(QPTR+I), I = 0, NGRID-1 )
      IF( SDIN .NE. 6 .AND. SDIN .NE. 0 ) CLOSE( SDIN )
*
*     Fatal error if we've got an illegal grid
*
      DO 70 J = 0, NGRID-1
         IF( MEM(PPTR+J).LT.1 .OR. MEM(QPTR+J).LT.1 ) THEN
            WRITE(OUTNUM, 4000) MEM(PPTR+J), MEM(QPTR+J)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   70 CONTINUE
*
*     Prepare output variables
*
      MEM(MEMUSED)   = NSHAPE
      MEM(MEMUSED+1) = NMAT
      MEM(MEMUSED+2) = NSRC
      MEM(MEMUSED+3) = NGRID
      MEMUSED = MEMUSED + 3
      CMEMUSED = CMEMUSED - 1
*
 1000 FORMAT('Mem too short (',I4,') to handle',I4,' ',A20)
 2000 FORMAT('Must have at least one ',A20)
 3000 FORMAT('UNRECOGNIZABLE ',A5,' ''', A1, '''.')
 4000 FORMAT('Illegal process grid: {',I3,',',I3,'}.')
*
      RETURN
*
*     End of RDSDRV.
*
      END
*
      SUBROUTINE CHKMATDAT( NOUT, INFILE, TSTFLAG, NMAT, M0, N0,
     $                      LDAS0, LDAD0, LDI0 )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL TSTFLAG
      INTEGER NOUT, NMAT
*     ..
*     .. Array Arguments ..
      CHARACTER*8 INFILE
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
*     ..
*   Purpose
*  =======
*  CHKMATDAT: Checks that matrix data is correct.
*
*  Arguments
*  =========
*  NOUT    (input) INTEGER
*          The device number to write output to.
*
*  INFILE  (input) CHARACTER*8
*          The name of the input file where matrix values came from.
*
*  TSTFLAG (input) LOGICAL
*          Whether to test RCFLAG (LDI) values or not.
*
*  NMAT    (input) INTEGER
*          The number of matrices to be tested.
*
*  M0      (input) INTEGER array of dimension (NMAT)
*          Values of M to be tested.
*
*  M0      (input) INTEGER array of dimension (NMAT)
*          Values of M to be tested.
*
*  N0      (input) INTEGER array of dimension (NMAT)
*          Values of N to be tested.
*
*  LDAS0   (input) INTEGER array of dimension (NMAT)
*          Values of LDAS (leading dimension of A on source process)
*          to be tested.
*
*  LDAD0   (input) INTEGER array of dimension (NMAT)
*          Values of LDAD (leading dimension of A on destination
*          process) to be tested.
*
*  ====================================================================
*
*     .. Local Scalars ..
      LOGICAL MATOK
      INTEGER I
*     ..
*     .. Executable Statements ..
      MATOK = .TRUE.
      DO 10 I = 1, NMAT
         IF( M0(I) .LT. 0 ) THEN
            WRITE(NOUT,1000) INFILE, 'M', M0(I)
            MATOK = .FALSE.
         ELSE IF( N0(I) .LT. 0 ) THEN
            WRITE(NOUT,1000) INFILE, 'N', N0(I)
            MATOK = .FALSE.
         ELSE IF( LDAS0(I) .LT. M0(I) ) THEN
            WRITE(NOUT,2000) INFILE, 'LDASRC', LDAS0(I), M0(I)
            MATOK = .FALSE.
         ELSE IF( LDAD0(I) .LT. M0(I) ) THEN
            WRITE(NOUT,2000) INFILE, 'LDADST', LDAD0(I), M0(I)
            MATOK = .FALSE.
         ELSE IF( TSTFLAG ) THEN
            IF( (LDI0(I).LT.M0(I)) .AND. (LDI0(I).NE.-1) ) THEN
               WRITE(NOUT,2000) INFILE, 'RCFLAG', LDI0(I), M0(I)
               MATOK = .FALSE.
            END IF
         END IF
   10 CONTINUE
*
      IF( .NOT.MATOK ) THEN
         IF( NOUT .NE. 6 .AND. NOUT .NE. 0 ) CLOSE(NOUT)
         CALL BLACS_ABORT(-1, 1)
      END IF
*
 1000 FORMAT(A8,' INPUT ERROR: Illegal ',A1,'; value=',I6,'.')
 2000 FORMAT(A8,' INPUT ERROR: Illegal ',A6,'; value=',I6,', but M=',I6)
*
      RETURN
      END
*
      LOGICAL FUNCTION ALLPASS( THISTEST )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL THISTEST
*     ..
*  Purpose
*  =======
*  ALLPASS: Returns whether all tests have passed so far.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL PASSHIST
*     ..
*     .. Save Statement ..
      SAVE PASSHIST
*     ..
*     .. Data Statements ..
      DATA PASSHIST /.TRUE./
*     ..
*     .. Executable Statements ..
      PASSHIST = (PASSHIST .AND. THISTEST)
      ALLPASS = PASSHIST
*
      RETURN
      END
*
      SUBROUTINE RDBSBR( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMUSED, MEMLEN, CMEMUSED, CMEMLEN, OUTNUM
*     ..
*     .. Array Arguments ..
      CHARACTER*1 CMEM(CMEMLEN)
      INTEGER MEM(MEMLEN)
*     ..
*
*     Purpose
*     =======
*     RDBSBR:  Read and process the input file BSBR.dat.
*
*     Arguments
*     =========
*     MEMUSED  (output) INTEGER
*              Number of elements in MEM that this subroutine ends up using.
*
*     MEM      (output) INTEGER array of dimension memlen
*              On output, holds information read in from sdrv.dat.
*
*     MEMLEN   (input) INTEGER
*              Number of elements of MEM that this subroutine
*              may safely write into.
*
*     CMEMUSED (output) INTEGER
*              Number of elements in CMEM that this subroutine ends up using.
*
*     CMEM     (output) CHARACTER*1 array of dimension cmemlen
*              On output, holds the values for UPLO and DIAG.
*
*     CMEMLEN  (input) INTEGER
*              Number of elements of CMEM that this subroutine
*              may safely write into.
*
*     OUTNUM   (input) INTEGER
*              Unit number of the output file.
*
*     =================================================================
*
*     .. Parameters ..
      INTEGER SDIN
      PARAMETER( SDIN = 12 )
*     ..
*     .. External Functions ..
      LOGICAL  LSAME
      EXTERNAL LSAME
*     ..
*     .. Local Scalars ..
      INTEGER NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID, I, J
      INTEGER SCOPEPTR, TOPPTR, UPLOPTR, DIAGPTR, MPTR, NPTR
      INTEGER LDSPTR, LDDPTR, RSRCPTR, CSRCPTR, PPTR, QPTR
*     ..
*     .. Executable Statements
*
*     Open and read the file bsbr.dat.  The expected format is
*     below.
*
*------
*integer                         Number of scopes
*array of CHAR*1's               Values for Scopes
*integer                         Number of topologies
*array of CHAR*1's               Values for TOP
*integer                         number of shapes of the matrix
*array of CHAR*1's               UPLO
*array of CHAR*1's               DIAG: unit diagonal or not?
*integer                         number of nmat
*array of integers               M: number of rows in matrix
*array of integers               N: number of columns in matrix
*integer                         LDA: leading dimension on source proc
*integer                         LDA: leading dimension on dest proc
*integer                         number of source/dest pairs
*array of integers               RSRC: process row of message source
*array of integers               CSRC: process column of msg. src.
*integer                         Number of grids
*array of integers               NPROW: number of rows in process grid
*array of integers               NPCOL: number of col's in proc. grid
*------
*  note: UPLO stands for 'upper or lower trapezoidal or general
*        rectangular.'
*  note: the text descriptions as shown above are present in
*             the sample bsbr.dat included with this distribution,
*             but are not required.
*
*     Read input file
*
      MEMUSED = 1
      CMEMUSED = 1
      OPEN(UNIT = SDIN, FILE = 'bsbr.dat', STATUS = 'OLD')
*
*     Read in scopes and topologies
*
      READ(SDIN, *) NSCOPE
      SCOPEPTR = CMEMUSED
      CMEMUSED = SCOPEPTR + NSCOPE
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NSCOPE, 'SCOPES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSCOPE .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'SCOPE.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
*
      READ(SDIN, *) ( CMEM(SCOPEPTR+I), I = 0, NSCOPE-1 )
      DO 20 I = 0, NSCOPE-1
         IF( LSAME(CMEM(SCOPEPTR+I), 'R') ) THEN
            CMEM(SCOPEPTR+I) = 'R'
         ELSE IF( LSAME(CMEM(SCOPEPTR+I), 'C') ) THEN
            CMEM(SCOPEPTR+I) = 'C'
         ELSE IF( LSAME(CMEM(SCOPEPTR+I), 'A') ) THEN
            CMEM(SCOPEPTR+I) = 'A'
         ELSE
            WRITE(OUTNUM, 3000) 'SCOPE', CMEM(SCOPEPTR+I)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   20 CONTINUE
*
      READ(SDIN, *) NTOP
      TOPPTR = CMEMUSED
      CMEMUSED = TOPPTR + NTOP
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NTOP, 'TOPOLOGIES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NTOP .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'TOPOLOGY.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( CMEM(TOPPTR+I), I = 0, NTOP-1 )
*
*
*     Read in number of shapes, and values of UPLO and DIAG
*
      READ(SDIN, *) NSHAPE
      UPLOPTR = CMEMUSED
      DIAGPTR = UPLOPTR + NSHAPE
      CMEMUSED = DIAGPTR + NSHAPE
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NSHAPE, 'MATRIX SHAPES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSHAPE .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'MATRIX SHAPE.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
*
*     Read in, upcase, and fatal error if UPLO/DIAG not recognized
*
      READ(SDIN, *) ( CMEM(UPLOPTR+I), I = 0, NSHAPE-1 )
      DO 30 I = 0, NSHAPE-1
         IF( LSAME(CMEM(UPLOPTR+I), 'G') ) THEN
            CMEM(UPLOPTR+I) = 'G'
         ELSE IF( LSAME(CMEM(UPLOPTR+I), 'U') ) THEN
            CMEM(UPLOPTR+I) = 'U'
         ELSE IF( LSAME(CMEM(UPLOPTR+I), 'L') ) THEN
            CMEM(UPLOPTR+I) = 'L'
         ELSE
            WRITE(OUTNUM, 3000) 'UPLO ', CMEM(UPLOPTR+I)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   30 CONTINUE
*
      READ(SDIN, *) ( CMEM(DIAGPTR+I), I = 0, NSHAPE-1 )
      DO 40 I = 0, NSHAPE-1
         IF( CMEM(UPLOPTR+I) .NE. 'G' ) THEN
            IF( LSAME(CMEM(DIAGPTR+I), 'U') ) THEN
               CMEM( DIAGPTR+I ) = 'U'
            ELSE IF( LSAME(CMEM(DIAGPTR+I), 'N') ) THEN
               CMEM(DIAGPTR+I) = 'N'
            ELSE
               WRITE(OUTNUM, 3000) 'DIAG ', CMEM(DIAGPTR+I)
               IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
               STOP
            END IF
         END IF
   40 CONTINUE
*
*     Read in number of matrices, and values for M, N, LDASRC, and LDADEST
*
      READ(SDIN, *) NMAT
      MPTR = MEMUSED
      NPTR = MPTR + NMAT
      LDSPTR = NPTR + NMAT
      LDDPTR = LDSPTR + NMAT
      MEMUSED = LDDPTR + NMAT
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'MATRICES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NMAT .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'MATRIX.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM( MPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( NPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDSPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDDPTR+I ), I = 0, NMAT-1 )
*
*     Make sure matrix values are legal
*
      CALL CHKMATDAT( OUTNUM, 'BSBR.dat', .FALSE., NMAT, MEM(MPTR),
     $                MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR), MEM(LDDPTR) )
*
*     Read in number of src pairs, and values of src
*
      READ(SDIN, *) NSRC
      RSRCPTR  = MEMUSED
      CSRCPTR  = RSRCPTR  + NSRC
      MEMUSED  = CSRCPTR + NSRC
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'SRC.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSRC .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'SRC.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM(RSRCPTR+I), I = 0, NSRC-1 )
      READ(SDIN, *) ( MEM(CSRCPTR+I), I = 0, NSRC-1 )
*
*     Read in number of grids pairs, and values of P (process rows) and
*     Q (process columns)
*
      READ(SDIN, *) NGRID
      PPTR = MEMUSED
      QPTR = PPTR + NGRID
      MEMUSED = QPTR + NGRID
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NGRID, 'PROCESS GRIDS.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NGRID .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'PROCESS GRID'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE( OUTNUM )
         STOP
      END IF
*
      READ(SDIN, *) ( MEM(PPTR+I), I = 0, NGRID-1 )
      READ(SDIN, *) ( MEM(QPTR+I), I = 0, NGRID-1 )
      IF( SDIN .NE. 6 .AND. SDIN .NE. 0 ) CLOSE( SDIN )
*
*     Fatal error if we've got an illegal grid
*
      DO 70 J = 0, NGRID-1
         IF( MEM(PPTR+J).LT.1 .OR. MEM(QPTR+J).LT.1 ) THEN
            WRITE(OUTNUM, 4000) MEM(PPTR+J), MEM(QPTR+J)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   70 CONTINUE
*
*     Prepare output variables
*
      MEM(MEMUSED)   = NSCOPE
      MEM(MEMUSED+1) = NTOP
      MEM(MEMUSED+2) = NSHAPE
      MEM(MEMUSED+3) = NMAT
      MEM(MEMUSED+4) = NSRC
      MEM(MEMUSED+5) = NGRID
      MEMUSED = MEMUSED + 5
      CMEMUSED = CMEMUSED - 1
*
 1000 FORMAT('Mem too short (',I4,') to handle',I4,' ',A20)
 2000 FORMAT('Must have at least one ',A20)
 3000 FORMAT('UNRECOGNIZABLE ',A5,' ''', A1, '''.')
 4000 FORMAT('Illegal process grid: {',I3,',',I3,'}.')
*
      RETURN
*
*     End of RDBSBR.
*
      END
*
*
      SUBROUTINE ISDRVTEST( OUTNUM, VERB, NSHAPE, UPLO0, DIAG0,
     $                      NMAT, M0, N0, LDAS0, LDAD0, NSRC, RSRC0,
     $                      CSRC0, RDEST0, CDEST0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSHAPE, NMAT, NSRC, NGRID, MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), RDEST0(NSRC), CDEST0(NSRC)
      INTEGER CONTEXT0(NGRID), P0(NGRID), Q0(NGRID), TFAIL(*)
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ITESTSDRV:  Test integer send/recv
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) INTEGER array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL ITRSD2D, IGESD2D, ITRRV2D, IGERV2D
      EXTERNAL IINITMAT, ICHKMAT, ICHKPAD, IBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 UPLO, DIAG
      LOGICAL TESTOK
      INTEGER IAM, I, K, IGR, ISH, IMA, ISO, MYROW, MYCOL, IPRE, IPOST
      INTEGER M, N, NPROW, NPCOL, RSRC, CSRC, RDEST, CDEST
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE
      INTEGER SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -1
      RCHECKVAL = -2
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      ISIZE = IBTSIZEOF('I')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( ISIZE * (MEMLEN-I) ) / ( ISIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SDRV tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         DO 80 ISH = 1, NSHAPE
            UPLO = UPLO0(ISH)
            DIAG = DIAG0(ISH)
*
            DO 70 IMA = 1, NMAT
               M = M0(IMA)
               N = N0(IMA)
               LDASRC = LDAS0(IMA)
               LDADST = LDAD0(IMA)
*
               DO 60 ISO = 1, NSRC
                  TESTNUM = TESTNUM + 1
                  RSRC = RSRC0(ISO)
                  CSRC = CSRC0(ISO)
                  IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
                  RDEST = RDEST0(ISO)
                  CDEST = CDEST0(ISO)
                  IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     IF( IAM .EQ. 0 ) THEN
                        WRITE(OUTNUM, 7000) TESTNUM, 'RUNNING',
     $                                      UPLO, DIAG, M, N,
     $                                      LDASRC, LDADST, RSRC, CSRC,
     $                                      RDEST, CDEST, NPROW, NPCOL
                     END IF
                  END IF
*
                  TESTOK = .TRUE.
                  IPRE  = 2 * M
                  IPOST = IPRE
                  APTR = IPRE + 1
*
*                 source process generates matrix and sends it
*
                  IF( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
                     CALL IINITMAT( UPLO, DIAG, M, N, MEM, LDASRC,
     $                              IPRE, IPOST, SCHECKVAL, TESTNUM,
     $                              MYROW, MYCOL )
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                         CALL ITRSD2D( CONTEXT, UPLO, DIAG, M, N,
     $                                 MEM(APTR), LDASRC, RDEST, CDEST )
                     ELSE
                         CALL IGESD2D( CONTEXT, M, N, MEM(APTR),
     $                                 LDASRC, RDEST, CDEST )
                     END IF
                  END IF
*
                  IF( MYROW .EQ. RDEST .AND. MYCOL .EQ. CDEST ) THEN
*
*                    Pad entire matrix area
*
                     DO 50 K = 1, IPRE+IPOST+LDADST*N
                        MEM(K) = RCHECKVAL
   50                CONTINUE
*
*                    Receive matrix
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                        CALL ITRRV2D( CONTEXT, UPLO, DIAG, M, N,
     $                                MEM(APTR), LDADST, RSRC, CSRC )
                     ELSE
                        CALL IGERV2D( CONTEXT, M, N, MEM(APTR),
     $                                LDADST, RSRC, CSRC )
                     END IF
*
*                    Check for errors in matrix or padding
*
                     I = NERR
                     CALL ICHKMAT( UPLO, DIAG, M, N, MEM(APTR), LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, TESTNUM, MAXERR,
     $                        NERR, MEM(ERRIPTR), MEM(ERRDPTR) )
*
                     CALL ICHKPAD( UPLO, DIAG, M, N, MEM, LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST,
     $                        RCHECKVAL, TESTNUM, MAXERR, NERR,
     $                        MEM(ERRIPTR), MEM(ERRDPTR) )
                     TESTOK = I .EQ. NERR
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     I = NERR
                     CALL IBTCHECKIN( 0, OUTNUM, MAXERR, NERR,
     $                                MEM(ERRIPTR), MEM(ERRDPTR),
     $                                TFAIL )
                     IF( IAM .EQ. 0 ) THEN
                        IF( TESTOK .AND. I.EQ.NERR ) THEN
                           WRITE(OUTNUM, 7000) TESTNUM, 'PASSED ',
     $                           UPLO, DIAG, M, N, LDASRC, LDADST,
     $                           RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ELSE
                           NFAIL = NFAIL + 1
                           WRITE(OUTNUM, 7000) TESTNUM, 'FAILED ',
     $                          UPLO, DIAG, M, N, LDASRC, LDADST,
     $                          RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ENDIF
                     END IF
*
*                    Once we've printed out errors, can re-use buf space
*
                     NERR = 0
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL IBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('INTEGER SDRV TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS UPLO DIA     M     N  LDAS  LDAD RSRC ',
     $       'CSRC RDEST CDEST    P    Q')
 6000 FORMAT(' ----- ------- ---- --- ----- ----- ----- ----- ---- ',
     $       '---- ----- ----- ---- ----')
 7000 FORMAT(I6,1X,A7,4X,A1,3X,A1,4I6,2I5,2I6,2I5)
 8000 FORMAT('INTEGER SDRV TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('INTEGER SDRV TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ISDRVTEST.
*
      END
*
*
      SUBROUTINE SSDRVTEST( OUTNUM, VERB, NSHAPE, UPLO0, DIAG0,
     $                      NMAT, M0, N0, LDAS0, LDAD0, NSRC, RSRC0,
     $                      CSRC0, RDEST0, CDEST0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSHAPE, NMAT, NSRC, NGRID, MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), RDEST0(NSRC), CDEST0(NSRC)
      INTEGER CONTEXT0(NGRID), P0(NGRID), Q0(NGRID), TFAIL(*)
      REAL MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  STESTSDRV:  Test real send/recv
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) REAL array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL STRSD2D, SGESD2D, STRRV2D, SGERV2D
      EXTERNAL SINITMAT, SCHKMAT, SCHKPAD, SBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 UPLO, DIAG
      LOGICAL TESTOK
      INTEGER IAM, I, K, IGR, ISH, IMA, ISO, MYROW, MYCOL, IPRE, IPOST
      INTEGER M, N, NPROW, NPCOL, RSRC, CSRC, RDEST, CDEST
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, SSIZE
      REAL SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -0.01E0
      RCHECKVAL = -0.02E0
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( SSIZE * (MEMLEN-I) ) / ( SSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SDRV tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         DO 80 ISH = 1, NSHAPE
            UPLO = UPLO0(ISH)
            DIAG = DIAG0(ISH)
*
            DO 70 IMA = 1, NMAT
               M = M0(IMA)
               N = N0(IMA)
               LDASRC = LDAS0(IMA)
               LDADST = LDAD0(IMA)
*
               DO 60 ISO = 1, NSRC
                  TESTNUM = TESTNUM + 1
                  RSRC = RSRC0(ISO)
                  CSRC = CSRC0(ISO)
                  IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
                  RDEST = RDEST0(ISO)
                  CDEST = CDEST0(ISO)
                  IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     IF( IAM .EQ. 0 ) THEN
                        WRITE(OUTNUM, 7000) TESTNUM, 'RUNNING',
     $                                      UPLO, DIAG, M, N,
     $                                      LDASRC, LDADST, RSRC, CSRC,
     $                                      RDEST, CDEST, NPROW, NPCOL
                     END IF
                  END IF
*
                  TESTOK = .TRUE.
                  IPRE  = 2 * M
                  IPOST = IPRE
                  APTR = IPRE + 1
*
*                 source process generates matrix and sends it
*
                  IF( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
                     CALL SINITMAT( UPLO, DIAG, M, N, MEM, LDASRC,
     $                              IPRE, IPOST, SCHECKVAL, TESTNUM,
     $                              MYROW, MYCOL )
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                         CALL STRSD2D( CONTEXT, UPLO, DIAG, M, N,
     $                                 MEM(APTR), LDASRC, RDEST, CDEST )
                     ELSE
                         CALL SGESD2D( CONTEXT, M, N, MEM(APTR),
     $                                 LDASRC, RDEST, CDEST )
                     END IF
                  END IF
*
                  IF( MYROW .EQ. RDEST .AND. MYCOL .EQ. CDEST ) THEN
*
*                    Pad entire matrix area
*
                     DO 50 K = 1, IPRE+IPOST+LDADST*N
                        MEM(K) = RCHECKVAL
   50                CONTINUE
*
*                    Receive matrix
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                        CALL STRRV2D( CONTEXT, UPLO, DIAG, M, N,
     $                                MEM(APTR), LDADST, RSRC, CSRC )
                     ELSE
                        CALL SGERV2D( CONTEXT, M, N, MEM(APTR),
     $                                LDADST, RSRC, CSRC )
                     END IF
*
*                    Check for errors in matrix or padding
*
                     I = NERR
                     CALL SCHKMAT( UPLO, DIAG, M, N, MEM(APTR), LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, TESTNUM, MAXERR,
     $                        NERR, MEM(ERRIPTR), MEM(ERRDPTR) )
*
                     CALL SCHKPAD( UPLO, DIAG, M, N, MEM, LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST,
     $                        RCHECKVAL, TESTNUM, MAXERR, NERR,
     $                        MEM(ERRIPTR), MEM(ERRDPTR) )
                     TESTOK = I .EQ. NERR
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     I = NERR
                     CALL SBTCHECKIN( 0, OUTNUM, MAXERR, NERR,
     $                                MEM(ERRIPTR), MEM(ERRDPTR),
     $                                TFAIL )
                     IF( IAM .EQ. 0 ) THEN
                        IF( TESTOK .AND. I.EQ.NERR ) THEN
                           WRITE(OUTNUM, 7000) TESTNUM, 'PASSED ',
     $                           UPLO, DIAG, M, N, LDASRC, LDADST,
     $                           RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ELSE
                           NFAIL = NFAIL + 1
                           WRITE(OUTNUM, 7000) TESTNUM, 'FAILED ',
     $                          UPLO, DIAG, M, N, LDASRC, LDADST,
     $                          RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ENDIF
                     END IF
*
*                    Once we've printed out errors, can re-use buf space
*
                     NERR = 0
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL SBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('REAL SDRV TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS UPLO DIA     M     N  LDAS  LDAD RSRC ',
     $       'CSRC RDEST CDEST    P    Q')
 6000 FORMAT(' ----- ------- ---- --- ----- ----- ----- ----- ---- ',
     $       '---- ----- ----- ---- ----')
 7000 FORMAT(I6,1X,A7,4X,A1,3X,A1,4I6,2I5,2I6,2I5)
 8000 FORMAT('REAL SDRV TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('REAL SDRV TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of SSDRVTEST.
*
      END
*
*
      SUBROUTINE DSDRVTEST( OUTNUM, VERB, NSHAPE, UPLO0, DIAG0,
     $                      NMAT, M0, N0, LDAS0, LDAD0, NSRC, RSRC0,
     $                      CSRC0, RDEST0, CDEST0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSHAPE, NMAT, NSRC, NGRID, MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), RDEST0(NSRC), CDEST0(NSRC)
      INTEGER CONTEXT0(NGRID), P0(NGRID), Q0(NGRID), TFAIL(*)
      DOUBLE PRECISION MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  DTESTSDRV:  Test double precision send/recv
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) DOUBLE PRECISION array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL DTRSD2D, DGESD2D, DTRRV2D, DGERV2D
      EXTERNAL DINITMAT, DCHKMAT, DCHKPAD, DBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 UPLO, DIAG
      LOGICAL TESTOK
      INTEGER IAM, I, K, IGR, ISH, IMA, ISO, MYROW, MYCOL, IPRE, IPOST
      INTEGER M, N, NPROW, NPCOL, RSRC, CSRC, RDEST, CDEST
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, DSIZE
      DOUBLE PRECISION SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -0.01D0
      RCHECKVAL = -0.02D0
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( DSIZE * (MEMLEN-I) ) / ( DSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SDRV tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         DO 80 ISH = 1, NSHAPE
            UPLO = UPLO0(ISH)
            DIAG = DIAG0(ISH)
*
            DO 70 IMA = 1, NMAT
               M = M0(IMA)
               N = N0(IMA)
               LDASRC = LDAS0(IMA)
               LDADST = LDAD0(IMA)
*
               DO 60 ISO = 1, NSRC
                  TESTNUM = TESTNUM + 1
                  RSRC = RSRC0(ISO)
                  CSRC = CSRC0(ISO)
                  IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
                  RDEST = RDEST0(ISO)
                  CDEST = CDEST0(ISO)
                  IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     IF( IAM .EQ. 0 ) THEN
                        WRITE(OUTNUM, 7000) TESTNUM, 'RUNNING',
     $                                      UPLO, DIAG, M, N,
     $                                      LDASRC, LDADST, RSRC, CSRC,
     $                                      RDEST, CDEST, NPROW, NPCOL
                     END IF
                  END IF
*
                  TESTOK = .TRUE.
                  IPRE  = 2 * M
                  IPOST = IPRE
                  APTR = IPRE + 1
*
*                 source process generates matrix and sends it
*
                  IF( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
                     CALL DINITMAT( UPLO, DIAG, M, N, MEM, LDASRC,
     $                              IPRE, IPOST, SCHECKVAL, TESTNUM,
     $                              MYROW, MYCOL )
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                         CALL DTRSD2D( CONTEXT, UPLO, DIAG, M, N,
     $                                 MEM(APTR), LDASRC, RDEST, CDEST )
                     ELSE
                         CALL DGESD2D( CONTEXT, M, N, MEM(APTR),
     $                                 LDASRC, RDEST, CDEST )
                     END IF
                  END IF
*
                  IF( MYROW .EQ. RDEST .AND. MYCOL .EQ. CDEST ) THEN
*
*                    Pad entire matrix area
*
                     DO 50 K = 1, IPRE+IPOST+LDADST*N
                        MEM(K) = RCHECKVAL
   50                CONTINUE
*
*                    Receive matrix
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                        CALL DTRRV2D( CONTEXT, UPLO, DIAG, M, N,
     $                                MEM(APTR), LDADST, RSRC, CSRC )
                     ELSE
                        CALL DGERV2D( CONTEXT, M, N, MEM(APTR),
     $                                LDADST, RSRC, CSRC )
                     END IF
*
*                    Check for errors in matrix or padding
*
                     I = NERR
                     CALL DCHKMAT( UPLO, DIAG, M, N, MEM(APTR), LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, TESTNUM, MAXERR,
     $                        NERR, MEM(ERRIPTR), MEM(ERRDPTR) )
*
                     CALL DCHKPAD( UPLO, DIAG, M, N, MEM, LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST,
     $                        RCHECKVAL, TESTNUM, MAXERR, NERR,
     $                        MEM(ERRIPTR), MEM(ERRDPTR) )
                     TESTOK = I .EQ. NERR
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     I = NERR
                     CALL DBTCHECKIN( 0, OUTNUM, MAXERR, NERR,
     $                                MEM(ERRIPTR), MEM(ERRDPTR),
     $                                TFAIL )
                     IF( IAM .EQ. 0 ) THEN
                        IF( TESTOK .AND. I.EQ.NERR ) THEN
                           WRITE(OUTNUM, 7000) TESTNUM, 'PASSED ',
     $                           UPLO, DIAG, M, N, LDASRC, LDADST,
     $                           RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ELSE
                           NFAIL = NFAIL + 1
                           WRITE(OUTNUM, 7000) TESTNUM, 'FAILED ',
     $                          UPLO, DIAG, M, N, LDASRC, LDADST,
     $                          RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ENDIF
                     END IF
*
*                    Once we've printed out errors, can re-use buf space
*
                     NERR = 0
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL DBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE PRECISION SDRV TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS UPLO DIA     M     N  LDAS  LDAD RSRC ',
     $       'CSRC RDEST CDEST    P    Q')
 6000 FORMAT(' ----- ------- ---- --- ----- ----- ----- ----- ---- ',
     $       '---- ----- ----- ---- ----')
 7000 FORMAT(I6,1X,A7,4X,A1,3X,A1,4I6,2I5,2I6,2I5)
 8000 FORMAT('DOUBLE PRECISION SDRV TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('DOUBLE PRECISION SDRV TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of DSDRVTEST.
*
      END
*
*
      SUBROUTINE CSDRVTEST( OUTNUM, VERB, NSHAPE, UPLO0, DIAG0,
     $                      NMAT, M0, N0, LDAS0, LDAD0, NSRC, RSRC0,
     $                      CSRC0, RDEST0, CDEST0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSHAPE, NMAT, NSRC, NGRID, MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), RDEST0(NSRC), CDEST0(NSRC)
      INTEGER CONTEXT0(NGRID), P0(NGRID), Q0(NGRID), TFAIL(*)
      COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  CTESTSDRV:  Test complex send/recv
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL CTRSD2D, CGESD2D, CTRRV2D, CGERV2D
      EXTERNAL CINITMAT, CCHKMAT, CCHKPAD, CBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 UPLO, DIAG
      LOGICAL TESTOK
      INTEGER IAM, I, K, IGR, ISH, IMA, ISO, MYROW, MYCOL, IPRE, IPOST
      INTEGER M, N, NPROW, NPCOL, RSRC, CSRC, RDEST, CDEST
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, CSIZE
      COMPLEX SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = CMPLX( -0.01, -0.01 )
      RCHECKVAL = CMPLX( -0.02, -0.02 )
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      CSIZE = IBTSIZEOF('C')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( CSIZE * (MEMLEN-I) ) / ( CSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SDRV tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         DO 80 ISH = 1, NSHAPE
            UPLO = UPLO0(ISH)
            DIAG = DIAG0(ISH)
*
            DO 70 IMA = 1, NMAT
               M = M0(IMA)
               N = N0(IMA)
               LDASRC = LDAS0(IMA)
               LDADST = LDAD0(IMA)
*
               DO 60 ISO = 1, NSRC
                  TESTNUM = TESTNUM + 1
                  RSRC = RSRC0(ISO)
                  CSRC = CSRC0(ISO)
                  IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
                  RDEST = RDEST0(ISO)
                  CDEST = CDEST0(ISO)
                  IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     IF( IAM .EQ. 0 ) THEN
                        WRITE(OUTNUM, 7000) TESTNUM, 'RUNNING',
     $                                      UPLO, DIAG, M, N,
     $                                      LDASRC, LDADST, RSRC, CSRC,
     $                                      RDEST, CDEST, NPROW, NPCOL
                     END IF
                  END IF
*
                  TESTOK = .TRUE.
                  IPRE  = 2 * M
                  IPOST = IPRE
                  APTR = IPRE + 1
*
*                 source process generates matrix and sends it
*
                  IF( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
                     CALL CINITMAT( UPLO, DIAG, M, N, MEM, LDASRC,
     $                              IPRE, IPOST, SCHECKVAL, TESTNUM,
     $                              MYROW, MYCOL )
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                         CALL CTRSD2D( CONTEXT, UPLO, DIAG, M, N,
     $                                 MEM(APTR), LDASRC, RDEST, CDEST )
                     ELSE
                         CALL CGESD2D( CONTEXT, M, N, MEM(APTR),
     $                                 LDASRC, RDEST, CDEST )
                     END IF
                  END IF
*
                  IF( MYROW .EQ. RDEST .AND. MYCOL .EQ. CDEST ) THEN
*
*                    Pad entire matrix area
*
                     DO 50 K = 1, IPRE+IPOST+LDADST*N
                        MEM(K) = RCHECKVAL
   50                CONTINUE
*
*                    Receive matrix
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                        CALL CTRRV2D( CONTEXT, UPLO, DIAG, M, N,
     $                                MEM(APTR), LDADST, RSRC, CSRC )
                     ELSE
                        CALL CGERV2D( CONTEXT, M, N, MEM(APTR),
     $                                LDADST, RSRC, CSRC )
                     END IF
*
*                    Check for errors in matrix or padding
*
                     I = NERR
                     CALL CCHKMAT( UPLO, DIAG, M, N, MEM(APTR), LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, TESTNUM, MAXERR,
     $                        NERR, MEM(ERRIPTR), MEM(ERRDPTR) )
*
                     CALL CCHKPAD( UPLO, DIAG, M, N, MEM, LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST,
     $                        RCHECKVAL, TESTNUM, MAXERR, NERR,
     $                        MEM(ERRIPTR), MEM(ERRDPTR) )
                     TESTOK = I .EQ. NERR
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     I = NERR
                     CALL CBTCHECKIN( 0, OUTNUM, MAXERR, NERR,
     $                                MEM(ERRIPTR), MEM(ERRDPTR),
     $                                TFAIL )
                     IF( IAM .EQ. 0 ) THEN
                        IF( TESTOK .AND. I.EQ.NERR ) THEN
                           WRITE(OUTNUM, 7000) TESTNUM, 'PASSED ',
     $                           UPLO, DIAG, M, N, LDASRC, LDADST,
     $                           RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ELSE
                           NFAIL = NFAIL + 1
                           WRITE(OUTNUM, 7000) TESTNUM, 'FAILED ',
     $                          UPLO, DIAG, M, N, LDASRC, LDADST,
     $                          RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ENDIF
                     END IF
*
*                    Once we've printed out errors, can re-use buf space
*
                     NERR = 0
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL CBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('COMPLEX SDRV TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS UPLO DIA     M     N  LDAS  LDAD RSRC ',
     $       'CSRC RDEST CDEST    P    Q')
 6000 FORMAT(' ----- ------- ---- --- ----- ----- ----- ----- ---- ',
     $       '---- ----- ----- ---- ----')
 7000 FORMAT(I6,1X,A7,4X,A1,3X,A1,4I6,2I5,2I6,2I5)
 8000 FORMAT('COMPLEX SDRV TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('COMPLEX SDRV TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of CSDRVTEST.
*
      END
*
*
      SUBROUTINE ZSDRVTEST( OUTNUM, VERB, NSHAPE, UPLO0, DIAG0,
     $                      NMAT, M0, N0, LDAS0, LDAD0, NSRC, RSRC0,
     $                      CSRC0, RDEST0, CDEST0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSHAPE, NMAT, NSRC, NGRID, MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), RDEST0(NSRC), CDEST0(NSRC)
      INTEGER CONTEXT0(NGRID), P0(NGRID), Q0(NGRID), TFAIL(*)
      DOUBLE COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ZTESTSDRV:  Test double complex send/recv
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNSRC)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) DOUBLE COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL ALLPASS
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL ZTRSD2D, ZGESD2D, ZTRRV2D, ZGERV2D
      EXTERNAL ZINITMAT, ZCHKMAT, ZCHKPAD, ZBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 UPLO, DIAG
      LOGICAL TESTOK
      INTEGER IAM, I, K, IGR, ISH, IMA, ISO, MYROW, MYCOL, IPRE, IPOST
      INTEGER M, N, NPROW, NPCOL, RSRC, CSRC, RDEST, CDEST
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, ZSIZE
      DOUBLE COMPLEX SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = DCMPLX( -0.01D0, -0.01D0 )
      RCHECKVAL = DCMPLX( -0.02D0, -0.02D0 )
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      ZSIZE = IBTSIZEOF('Z')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( ZSIZE * (MEMLEN-I) ) / ( ZSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SDRV tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         DO 80 ISH = 1, NSHAPE
            UPLO = UPLO0(ISH)
            DIAG = DIAG0(ISH)
*
            DO 70 IMA = 1, NMAT
               M = M0(IMA)
               N = N0(IMA)
               LDASRC = LDAS0(IMA)
               LDADST = LDAD0(IMA)
*
               DO 60 ISO = 1, NSRC
                  TESTNUM = TESTNUM + 1
                  RSRC = RSRC0(ISO)
                  CSRC = CSRC0(ISO)
                  IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
                  RDEST = RDEST0(ISO)
                  CDEST = CDEST0(ISO)
                  IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                     NSKIP = NSKIP + 1
                     GOTO 60
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     IF( IAM .EQ. 0 ) THEN
                        WRITE(OUTNUM, 7000) TESTNUM, 'RUNNING',
     $                                      UPLO, DIAG, M, N,
     $                                      LDASRC, LDADST, RSRC, CSRC,
     $                                      RDEST, CDEST, NPROW, NPCOL
                     END IF
                  END IF
*
                  TESTOK = .TRUE.
                  IPRE  = 2 * M
                  IPOST = IPRE
                  APTR = IPRE + 1
*
*                 source process generates matrix and sends it
*
                  IF( MYROW .EQ. RSRC .AND. MYCOL .EQ. CSRC ) THEN
                     CALL ZINITMAT( UPLO, DIAG, M, N, MEM, LDASRC,
     $                              IPRE, IPOST, SCHECKVAL, TESTNUM,
     $                              MYROW, MYCOL )
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                         CALL ZTRSD2D( CONTEXT, UPLO, DIAG, M, N,
     $                                 MEM(APTR), LDASRC, RDEST, CDEST )
                     ELSE
                         CALL ZGESD2D( CONTEXT, M, N, MEM(APTR),
     $                                 LDASRC, RDEST, CDEST )
                     END IF
                  END IF
*
                  IF( MYROW .EQ. RDEST .AND. MYCOL .EQ. CDEST ) THEN
*
*                    Pad entire matrix area
*
                     DO 50 K = 1, IPRE+IPOST+LDADST*N
                        MEM(K) = RCHECKVAL
   50                CONTINUE
*
*                    Receive matrix
*
                     IF( UPLO .EQ. 'U' .OR. UPLO .EQ. 'L' ) THEN
                        CALL ZTRRV2D( CONTEXT, UPLO, DIAG, M, N,
     $                                MEM(APTR), LDADST, RSRC, CSRC )
                     ELSE
                        CALL ZGERV2D( CONTEXT, M, N, MEM(APTR),
     $                                LDADST, RSRC, CSRC )
                     END IF
*
*                    Check for errors in matrix or padding
*
                     I = NERR
                     CALL ZCHKMAT( UPLO, DIAG, M, N, MEM(APTR), LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, TESTNUM, MAXERR,
     $                        NERR, MEM(ERRIPTR), MEM(ERRDPTR) )
*
                     CALL ZCHKPAD( UPLO, DIAG, M, N, MEM, LDADST,
     $                        RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST,
     $                        RCHECKVAL, TESTNUM, MAXERR, NERR,
     $                        MEM(ERRIPTR), MEM(ERRDPTR) )
                     TESTOK = I .EQ. NERR
                  END IF
*
                  IF( VERB .GT. 1 ) THEN
                     I = NERR
                     CALL ZBTCHECKIN( 0, OUTNUM, MAXERR, NERR,
     $                                MEM(ERRIPTR), MEM(ERRDPTR),
     $                                TFAIL )
                     IF( IAM .EQ. 0 ) THEN
                        IF( TESTOK .AND. I.EQ.NERR ) THEN
                           WRITE(OUTNUM, 7000) TESTNUM, 'PASSED ',
     $                           UPLO, DIAG, M, N, LDASRC, LDADST,
     $                           RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ELSE
                           NFAIL = NFAIL + 1
                           WRITE(OUTNUM, 7000) TESTNUM, 'FAILED ',
     $                          UPLO, DIAG, M, N, LDASRC, LDADST,
     $                          RSRC, CSRC, RDEST, CDEST, NPROW, NPCOL
                        ENDIF
                     END IF
*
*                    Once we've printed out errors, can re-use buf space
*
                     NERR = 0
                  END IF
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL ZBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE COMPLEX SDRV TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS UPLO DIA     M     N  LDAS  LDAD RSRC ',
     $       'CSRC RDEST CDEST    P    Q')
 6000 FORMAT(' ----- ------- ---- --- ----- ----- ----- ----- ---- ',
     $       '---- ----- ----- ---- ----')
 7000 FORMAT(I6,1X,A7,4X,A1,3X,A1,4I6,2I5,2I6,2I5)
 8000 FORMAT('DOUBLE COMPLEX SDRV TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('DOUBLE COMPLEX SDRV TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ZSDRVTEST.
*
      END
*
*
      SUBROUTINE IBSBRTEST( OUTNUM, VERB, NSCOPE, SCOPE0, NTOP, TOP0,
     $                      NSHAPE, UPLO0, DIAG0, NMAT, M0, N0, LDAS0,
     $                      LDAD0, NSRC, RSRC0, CSRC0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID
      INTEGER MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), TFAIL(*)
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ITESTBSBR:  Test integer broadcast
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) INTEGER array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL ITRBS2D, IGEBS2D, ITRBR2D, IGEBR2D
      EXTERNAL IINITMAT, ICHKMAT, ICHKPAD, IBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP, UPLO, DIAG
      LOGICAL TESTOK, INGRID
      INTEGER IAM, I, K, J, IGR, ISH, IMA, ISO, ISC, ITO
      INTEGER M, N, NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
      INTEGER ISTART, ISTOP, IPRE, IPOST, SETWHAT
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE
      INTEGER SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -1
      RCHECKVAL = -2
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      ISIZE = IBTSIZEOF('I')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( ISIZE * (MEMLEN-I) ) / ( ISIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run BSBR tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         INGRID = ( NPROW .GT. 0 )
*
         DO 100 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 90 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multipath ('M') or general tree ('T'),
*              need to loop over calls to BLACS_SET
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 11
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 12
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 80 ISH = 1, NSHAPE
                  UPLO = UPLO0(ISH)
                  DIAG = DIAG0(ISH)
*
                  DO 70 IMA = 1, NMAT
                     M = M0(IMA)
                     N = N0(IMA)
                     LDASRC = LDAS0(IMA)
                     LDADST = LDAD0(IMA)
*
                     DO 60 ISO = 1, NSRC
                        TESTNUM = TESTNUM + 1
                        RSRC = RSRC0(ISO)
                        CSRC = CSRC0(ISO)
                        IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                           NSKIP = NSKIP + 1
                           GOTO 60
                        END IF
                        IF( VERB .GT. 1 ) THEN
                           IF( IAM .EQ. 0 ) THEN
                              WRITE(OUTNUM, 7000)
     $                        TESTNUM, 'RUNNING',SCOPE, TOP, UPLO, DIAG,
     $                        M, N, LDASRC, LDADST, RSRC, CSRC,
     $                        NPROW, NPCOL
                           END IF
                        END IF
*
                        TESTOK = .TRUE.
                        IPRE  = 2 * M
                        IPOST = IPRE
                        APTR = IPRE + 1
*
*                       If I am in scope
*
                        IF( (MYROW.EQ.RSRC .AND. SCOPE.EQ.'R') .OR.
     $                       (MYCOL.EQ.CSRC .AND. SCOPE.EQ.'C') .OR.
     $                       (SCOPE .EQ. 'A') ) THEN
*
*                          source process generates matrix and sends it
*
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL IINITMAT(UPLO, DIAG, M, N, MEM,
     $                                      LDASRC, IPRE, IPOST,
     $                                      SCHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              DO 20 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 20
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                     CALL ITRBS2D(CONTEXT, SCOPE, TOP,
     $                                            UPLO, DIAG, M, N,
     $                                            MEM(APTR), LDASRC )
                                 ELSE
                                     CALL IGEBS2D(CONTEXT, SCOPE, TOP,
     $                                            M, N, MEM(APTR),
     $                                            LDASRC )
                                 END IF
   20                         CONTINUE
*
*                          Destination processes
*
                           ELSE IF( INGRID ) THEN
                              DO 40 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 40
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*                                Pad entire matrix area
*
                                 DO 30 K = 1, IPRE+IPOST+LDADST*N
                                    MEM(K) = RCHECKVAL
   30                            CONTINUE
*
*                                Receive matrix
*
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                    CALL ITRBR2D(CONTEXT, SCOPE, TOP,
     $                                           UPLO, DIAG, M, N,
     $                                           MEM(APTR), LDADST,
     $                                           RSRC, CSRC)
                                 ELSE
                                    CALL IGEBR2D(CONTEXT, SCOPE, TOP,
     $                                           M, N, MEM(APTR),
     $                                           LDADST, RSRC, CSRC)
                                 END IF
*
*                                Check for errors in matrix or padding
*
                                 I = NERR
                                 CALL ICHKMAT(UPLO, DIAG, M, N,
     $                                   MEM(APTR), LDADST, RSRC, CSRC,
     $                                   MYROW, MYCOL, TESTNUM, MAXERR,
     $                                   NERR, MEM(ERRIPTR),
     $                                   MEM(ERRDPTR))
*
                                 CALL ICHKPAD(UPLO, DIAG, M, N, MEM,
     $                                   LDADST, RSRC, CSRC, MYROW,
     $                                   MYCOL, IPRE, IPOST, RCHECKVAL,
     $                                   TESTNUM, MAXERR, NERR,
     $                                   MEM(ERRIPTR), MEM(ERRDPTR))
   40                         CONTINUE
                              TESTOK = ( I .EQ. NERR )
                           END IF
                        END IF
*
                        IF( VERB .GT. 1 ) THEN
                           I = NERR
                           CALL IBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                                     MEM(ERRIPTR), MEM(ERRDPTR),
     $                                     TFAIL)
                           IF( IAM .EQ. 0 ) THEN
                              TESTOK = ( TESTOK .AND. (I.EQ.NERR) )
                              IF( TESTOK ) THEN
                                 WRITE(OUTNUM,7000)TESTNUM,'PASSED ',
     $                                 SCOPE, TOP, UPLO, DIAG, M, N,
     $                                 LDASRC, LDADST, RSRC, CSRC,
     $                                 NPROW, NPCOL
                              ELSE
                                 NFAIL = NFAIL + 1
                                 WRITE(OUTNUM,7000)TESTNUM,'FAILED ',
     $                                SCOPE, TOP, UPLO, DIAG, M, N,
     $                                LDASRC, LDADST, RSRC, CSRC,
     $                                NPROW, NPCOL
                              END IF
                           END IF
*
*                          Once we've printed out errors, can re-use buf space
*
                           NERR = 0
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL IBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('INTEGER BSBR TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS SCOPE TOP UPLO DIAG     M     N  LDAS ',
     $       ' LDAD RSRC CSRC    P    Q')
 6000 FORMAT(' ----- ------- ----- --- ---- ---- ----- ----- ----- ',
     $       '----- ---- ---- ---- ----')
 7000 FORMAT(I6,1X,A7,5X,A1,3X,A1,2(4X,A1), 4I6, 4I5)
 8000 FORMAT('INTEGER BSBR TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('INTEGER BSBR TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of IBSBRTEST.
*
      END
*
*
      SUBROUTINE SBSBRTEST( OUTNUM, VERB, NSCOPE, SCOPE0, NTOP, TOP0,
     $                      NSHAPE, UPLO0, DIAG0, NMAT, M0, N0, LDAS0,
     $                      LDAD0, NSRC, RSRC0, CSRC0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID
      INTEGER MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), TFAIL(*)
      REAL MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  STESTBSBR:  Test real broadcast
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) REAL array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL STRBS2D, SGEBS2D, STRBR2D, SGEBR2D
      EXTERNAL SINITMAT, SCHKMAT, SCHKPAD, SBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP, UPLO, DIAG
      LOGICAL TESTOK, INGRID
      INTEGER IAM, I, K, J, IGR, ISH, IMA, ISO, ISC, ITO
      INTEGER M, N, NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
      INTEGER ISTART, ISTOP, IPRE, IPOST, SETWHAT
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, SSIZE
      REAL SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -0.01E0
      RCHECKVAL = -0.02E0
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( SSIZE * (MEMLEN-I) ) / ( SSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run BSBR tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         INGRID = ( NPROW .GT. 0 )
*
         DO 100 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 90 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multipath ('M') or general tree ('T'),
*              need to loop over calls to BLACS_SET
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 11
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 12
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 80 ISH = 1, NSHAPE
                  UPLO = UPLO0(ISH)
                  DIAG = DIAG0(ISH)
*
                  DO 70 IMA = 1, NMAT
                     M = M0(IMA)
                     N = N0(IMA)
                     LDASRC = LDAS0(IMA)
                     LDADST = LDAD0(IMA)
*
                     DO 60 ISO = 1, NSRC
                        TESTNUM = TESTNUM + 1
                        RSRC = RSRC0(ISO)
                        CSRC = CSRC0(ISO)
                        IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                           NSKIP = NSKIP + 1
                           GOTO 60
                        END IF
                        IF( VERB .GT. 1 ) THEN
                           IF( IAM .EQ. 0 ) THEN
                              WRITE(OUTNUM, 7000)
     $                        TESTNUM, 'RUNNING',SCOPE, TOP, UPLO, DIAG,
     $                        M, N, LDASRC, LDADST, RSRC, CSRC,
     $                        NPROW, NPCOL
                           END IF
                        END IF
*
                        TESTOK = .TRUE.
                        IPRE  = 2 * M
                        IPOST = IPRE
                        APTR = IPRE + 1
*
*                       If I am in scope
*
                        IF( (MYROW.EQ.RSRC .AND. SCOPE.EQ.'R') .OR.
     $                       (MYCOL.EQ.CSRC .AND. SCOPE.EQ.'C') .OR.
     $                       (SCOPE .EQ. 'A') ) THEN
*
*                          source process generates matrix and sends it
*
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL SINITMAT(UPLO, DIAG, M, N, MEM,
     $                                      LDASRC, IPRE, IPOST,
     $                                      SCHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              DO 20 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 20
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                     CALL STRBS2D(CONTEXT, SCOPE, TOP,
     $                                            UPLO, DIAG, M, N,
     $                                            MEM(APTR), LDASRC )
                                 ELSE
                                     CALL SGEBS2D(CONTEXT, SCOPE, TOP,
     $                                            M, N, MEM(APTR),
     $                                            LDASRC )
                                 END IF
   20                         CONTINUE
*
*                          Destination processes
*
                           ELSE IF( INGRID ) THEN
                              DO 40 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 40
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*                                Pad entire matrix area
*
                                 DO 30 K = 1, IPRE+IPOST+LDADST*N
                                    MEM(K) = RCHECKVAL
   30                            CONTINUE
*
*                                Receive matrix
*
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                    CALL STRBR2D(CONTEXT, SCOPE, TOP,
     $                                           UPLO, DIAG, M, N,
     $                                           MEM(APTR), LDADST,
     $                                           RSRC, CSRC)
                                 ELSE
                                    CALL SGEBR2D(CONTEXT, SCOPE, TOP,
     $                                           M, N, MEM(APTR),
     $                                           LDADST, RSRC, CSRC)
                                 END IF
*
*                                Check for errors in matrix or padding
*
                                 I = NERR
                                 CALL SCHKMAT(UPLO, DIAG, M, N,
     $                                   MEM(APTR), LDADST, RSRC, CSRC,
     $                                   MYROW, MYCOL, TESTNUM, MAXERR,
     $                                   NERR, MEM(ERRIPTR),
     $                                   MEM(ERRDPTR))
*
                                 CALL SCHKPAD(UPLO, DIAG, M, N, MEM,
     $                                   LDADST, RSRC, CSRC, MYROW,
     $                                   MYCOL, IPRE, IPOST, RCHECKVAL,
     $                                   TESTNUM, MAXERR, NERR,
     $                                   MEM(ERRIPTR), MEM(ERRDPTR))
   40                         CONTINUE
                              TESTOK = ( I .EQ. NERR )
                           END IF
                        END IF
*
                        IF( VERB .GT. 1 ) THEN
                           I = NERR
                           CALL SBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                                     MEM(ERRIPTR), MEM(ERRDPTR),
     $                                     TFAIL)
                           IF( IAM .EQ. 0 ) THEN
                              TESTOK = ( TESTOK .AND. (I.EQ.NERR) )
                              IF( TESTOK ) THEN
                                 WRITE(OUTNUM,7000)TESTNUM,'PASSED ',
     $                                 SCOPE, TOP, UPLO, DIAG, M, N,
     $                                 LDASRC, LDADST, RSRC, CSRC,
     $                                 NPROW, NPCOL
                              ELSE
                                 NFAIL = NFAIL + 1
                                 WRITE(OUTNUM,7000)TESTNUM,'FAILED ',
     $                                SCOPE, TOP, UPLO, DIAG, M, N,
     $                                LDASRC, LDADST, RSRC, CSRC,
     $                                NPROW, NPCOL
                              END IF
                           END IF
*
*                          Once we've printed out errors, can re-use buf space
*
                           NERR = 0
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL SBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('REAL BSBR TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS SCOPE TOP UPLO DIAG     M     N  LDAS ',
     $       ' LDAD RSRC CSRC    P    Q')
 6000 FORMAT(' ----- ------- ----- --- ---- ---- ----- ----- ----- ',
     $       '----- ---- ---- ---- ----')
 7000 FORMAT(I6,1X,A7,5X,A1,3X,A1,2(4X,A1), 4I6, 4I5)
 8000 FORMAT('REAL BSBR TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('REAL BSBR TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of SBSBRTEST.
*
      END
*
*
      SUBROUTINE DBSBRTEST( OUTNUM, VERB, NSCOPE, SCOPE0, NTOP, TOP0,
     $                      NSHAPE, UPLO0, DIAG0, NMAT, M0, N0, LDAS0,
     $                      LDAD0, NSRC, RSRC0, CSRC0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID
      INTEGER MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), TFAIL(*)
      DOUBLE PRECISION MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  DTESTBSBR:  Test double precision broadcast
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) DOUBLE PRECISION array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL DTRBS2D, DGEBS2D, DTRBR2D, DGEBR2D
      EXTERNAL DINITMAT, DCHKMAT, DCHKPAD, DBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP, UPLO, DIAG
      LOGICAL TESTOK, INGRID
      INTEGER IAM, I, K, J, IGR, ISH, IMA, ISO, ISC, ITO
      INTEGER M, N, NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
      INTEGER ISTART, ISTOP, IPRE, IPOST, SETWHAT
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, DSIZE
      DOUBLE PRECISION SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = -0.01D0
      RCHECKVAL = -0.02D0
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( DSIZE * (MEMLEN-I) ) / ( DSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run BSBR tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         INGRID = ( NPROW .GT. 0 )
*
         DO 100 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 90 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multipath ('M') or general tree ('T'),
*              need to loop over calls to BLACS_SET
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 11
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 12
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 80 ISH = 1, NSHAPE
                  UPLO = UPLO0(ISH)
                  DIAG = DIAG0(ISH)
*
                  DO 70 IMA = 1, NMAT
                     M = M0(IMA)
                     N = N0(IMA)
                     LDASRC = LDAS0(IMA)
                     LDADST = LDAD0(IMA)
*
                     DO 60 ISO = 1, NSRC
                        TESTNUM = TESTNUM + 1
                        RSRC = RSRC0(ISO)
                        CSRC = CSRC0(ISO)
                        IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                           NSKIP = NSKIP + 1
                           GOTO 60
                        END IF
                        IF( VERB .GT. 1 ) THEN
                           IF( IAM .EQ. 0 ) THEN
                              WRITE(OUTNUM, 7000)
     $                        TESTNUM, 'RUNNING',SCOPE, TOP, UPLO, DIAG,
     $                        M, N, LDASRC, LDADST, RSRC, CSRC,
     $                        NPROW, NPCOL
                           END IF
                        END IF
*
                        TESTOK = .TRUE.
                        IPRE  = 2 * M
                        IPOST = IPRE
                        APTR = IPRE + 1
*
*                       If I am in scope
*
                        IF( (MYROW.EQ.RSRC .AND. SCOPE.EQ.'R') .OR.
     $                       (MYCOL.EQ.CSRC .AND. SCOPE.EQ.'C') .OR.
     $                       (SCOPE .EQ. 'A') ) THEN
*
*                          source process generates matrix and sends it
*
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL DINITMAT(UPLO, DIAG, M, N, MEM,
     $                                      LDASRC, IPRE, IPOST,
     $                                      SCHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              DO 20 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 20
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                     CALL DTRBS2D(CONTEXT, SCOPE, TOP,
     $                                            UPLO, DIAG, M, N,
     $                                            MEM(APTR), LDASRC )
                                 ELSE
                                     CALL DGEBS2D(CONTEXT, SCOPE, TOP,
     $                                            M, N, MEM(APTR),
     $                                            LDASRC )
                                 END IF
   20                         CONTINUE
*
*                          Destination processes
*
                           ELSE IF( INGRID ) THEN
                              DO 40 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 40
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*                                Pad entire matrix area
*
                                 DO 30 K = 1, IPRE+IPOST+LDADST*N
                                    MEM(K) = RCHECKVAL
   30                            CONTINUE
*
*                                Receive matrix
*
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                    CALL DTRBR2D(CONTEXT, SCOPE, TOP,
     $                                           UPLO, DIAG, M, N,
     $                                           MEM(APTR), LDADST,
     $                                           RSRC, CSRC)
                                 ELSE
                                    CALL DGEBR2D(CONTEXT, SCOPE, TOP,
     $                                           M, N, MEM(APTR),
     $                                           LDADST, RSRC, CSRC)
                                 END IF
*
*                                Check for errors in matrix or padding
*
                                 I = NERR
                                 CALL DCHKMAT(UPLO, DIAG, M, N,
     $                                   MEM(APTR), LDADST, RSRC, CSRC,
     $                                   MYROW, MYCOL, TESTNUM, MAXERR,
     $                                   NERR, MEM(ERRIPTR),
     $                                   MEM(ERRDPTR))
*
                                 CALL DCHKPAD(UPLO, DIAG, M, N, MEM,
     $                                   LDADST, RSRC, CSRC, MYROW,
     $                                   MYCOL, IPRE, IPOST, RCHECKVAL,
     $                                   TESTNUM, MAXERR, NERR,
     $                                   MEM(ERRIPTR), MEM(ERRDPTR))
   40                         CONTINUE
                              TESTOK = ( I .EQ. NERR )
                           END IF
                        END IF
*
                        IF( VERB .GT. 1 ) THEN
                           I = NERR
                           CALL DBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                                     MEM(ERRIPTR), MEM(ERRDPTR),
     $                                     TFAIL)
                           IF( IAM .EQ. 0 ) THEN
                              TESTOK = ( TESTOK .AND. (I.EQ.NERR) )
                              IF( TESTOK ) THEN
                                 WRITE(OUTNUM,7000)TESTNUM,'PASSED ',
     $                                 SCOPE, TOP, UPLO, DIAG, M, N,
     $                                 LDASRC, LDADST, RSRC, CSRC,
     $                                 NPROW, NPCOL
                              ELSE
                                 NFAIL = NFAIL + 1
                                 WRITE(OUTNUM,7000)TESTNUM,'FAILED ',
     $                                SCOPE, TOP, UPLO, DIAG, M, N,
     $                                LDASRC, LDADST, RSRC, CSRC,
     $                                NPROW, NPCOL
                              END IF
                           END IF
*
*                          Once we've printed out errors, can re-use buf space
*
                           NERR = 0
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL DBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE PRECISION BSBR TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS SCOPE TOP UPLO DIAG     M     N  LDAS ',
     $       ' LDAD RSRC CSRC    P    Q')
 6000 FORMAT(' ----- ------- ----- --- ---- ---- ----- ----- ----- ',
     $       '----- ---- ---- ---- ----')
 7000 FORMAT(I6,1X,A7,5X,A1,3X,A1,2(4X,A1), 4I6, 4I5)
 8000 FORMAT('DOUBLE PRECISION BSBR TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('DOUBLE PRECISION BSBR TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of DBSBRTEST.
*
      END
*
*
      SUBROUTINE CBSBRTEST( OUTNUM, VERB, NSCOPE, SCOPE0, NTOP, TOP0,
     $                      NSHAPE, UPLO0, DIAG0, NMAT, M0, N0, LDAS0,
     $                      LDAD0, NSRC, RSRC0, CSRC0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID
      INTEGER MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), TFAIL(*)
      COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  CTESTBSBR:  Test complex broadcast
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL CTRBS2D, CGEBS2D, CTRBR2D, CGEBR2D
      EXTERNAL CINITMAT, CCHKMAT, CCHKPAD, CBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP, UPLO, DIAG
      LOGICAL TESTOK, INGRID
      INTEGER IAM, I, K, J, IGR, ISH, IMA, ISO, ISC, ITO
      INTEGER M, N, NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
      INTEGER ISTART, ISTOP, IPRE, IPOST, SETWHAT
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, CSIZE
      COMPLEX SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = CMPLX( -0.01, -0.01 )
      RCHECKVAL = CMPLX( -0.02, -0.02 )
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      CSIZE = IBTSIZEOF('C')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( CSIZE * (MEMLEN-I) ) / ( CSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run BSBR tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         INGRID = ( NPROW .GT. 0 )
*
         DO 100 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 90 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multipath ('M') or general tree ('T'),
*              need to loop over calls to BLACS_SET
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 11
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 12
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 80 ISH = 1, NSHAPE
                  UPLO = UPLO0(ISH)
                  DIAG = DIAG0(ISH)
*
                  DO 70 IMA = 1, NMAT
                     M = M0(IMA)
                     N = N0(IMA)
                     LDASRC = LDAS0(IMA)
                     LDADST = LDAD0(IMA)
*
                     DO 60 ISO = 1, NSRC
                        TESTNUM = TESTNUM + 1
                        RSRC = RSRC0(ISO)
                        CSRC = CSRC0(ISO)
                        IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                           NSKIP = NSKIP + 1
                           GOTO 60
                        END IF
                        IF( VERB .GT. 1 ) THEN
                           IF( IAM .EQ. 0 ) THEN
                              WRITE(OUTNUM, 7000)
     $                        TESTNUM, 'RUNNING',SCOPE, TOP, UPLO, DIAG,
     $                        M, N, LDASRC, LDADST, RSRC, CSRC,
     $                        NPROW, NPCOL
                           END IF
                        END IF
*
                        TESTOK = .TRUE.
                        IPRE  = 2 * M
                        IPOST = IPRE
                        APTR = IPRE + 1
*
*                       If I am in scope
*
                        IF( (MYROW.EQ.RSRC .AND. SCOPE.EQ.'R') .OR.
     $                       (MYCOL.EQ.CSRC .AND. SCOPE.EQ.'C') .OR.
     $                       (SCOPE .EQ. 'A') ) THEN
*
*                          source process generates matrix and sends it
*
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL CINITMAT(UPLO, DIAG, M, N, MEM,
     $                                      LDASRC, IPRE, IPOST,
     $                                      SCHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              DO 20 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 20
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                     CALL CTRBS2D(CONTEXT, SCOPE, TOP,
     $                                            UPLO, DIAG, M, N,
     $                                            MEM(APTR), LDASRC )
                                 ELSE
                                     CALL CGEBS2D(CONTEXT, SCOPE, TOP,
     $                                            M, N, MEM(APTR),
     $                                            LDASRC )
                                 END IF
   20                         CONTINUE
*
*                          Destination processes
*
                           ELSE IF( INGRID ) THEN
                              DO 40 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 40
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*                                Pad entire matrix area
*
                                 DO 30 K = 1, IPRE+IPOST+LDADST*N
                                    MEM(K) = RCHECKVAL
   30                            CONTINUE
*
*                                Receive matrix
*
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                    CALL CTRBR2D(CONTEXT, SCOPE, TOP,
     $                                           UPLO, DIAG, M, N,
     $                                           MEM(APTR), LDADST,
     $                                           RSRC, CSRC)
                                 ELSE
                                    CALL CGEBR2D(CONTEXT, SCOPE, TOP,
     $                                           M, N, MEM(APTR),
     $                                           LDADST, RSRC, CSRC)
                                 END IF
*
*                                Check for errors in matrix or padding
*
                                 I = NERR
                                 CALL CCHKMAT(UPLO, DIAG, M, N,
     $                                   MEM(APTR), LDADST, RSRC, CSRC,
     $                                   MYROW, MYCOL, TESTNUM, MAXERR,
     $                                   NERR, MEM(ERRIPTR),
     $                                   MEM(ERRDPTR))
*
                                 CALL CCHKPAD(UPLO, DIAG, M, N, MEM,
     $                                   LDADST, RSRC, CSRC, MYROW,
     $                                   MYCOL, IPRE, IPOST, RCHECKVAL,
     $                                   TESTNUM, MAXERR, NERR,
     $                                   MEM(ERRIPTR), MEM(ERRDPTR))
   40                         CONTINUE
                              TESTOK = ( I .EQ. NERR )
                           END IF
                        END IF
*
                        IF( VERB .GT. 1 ) THEN
                           I = NERR
                           CALL CBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                                     MEM(ERRIPTR), MEM(ERRDPTR),
     $                                     TFAIL)
                           IF( IAM .EQ. 0 ) THEN
                              TESTOK = ( TESTOK .AND. (I.EQ.NERR) )
                              IF( TESTOK ) THEN
                                 WRITE(OUTNUM,7000)TESTNUM,'PASSED ',
     $                                 SCOPE, TOP, UPLO, DIAG, M, N,
     $                                 LDASRC, LDADST, RSRC, CSRC,
     $                                 NPROW, NPCOL
                              ELSE
                                 NFAIL = NFAIL + 1
                                 WRITE(OUTNUM,7000)TESTNUM,'FAILED ',
     $                                SCOPE, TOP, UPLO, DIAG, M, N,
     $                                LDASRC, LDADST, RSRC, CSRC,
     $                                NPROW, NPCOL
                              END IF
                           END IF
*
*                          Once we've printed out errors, can re-use buf space
*
                           NERR = 0
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL CBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('COMPLEX BSBR TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS SCOPE TOP UPLO DIAG     M     N  LDAS ',
     $       ' LDAD RSRC CSRC    P    Q')
 6000 FORMAT(' ----- ------- ----- --- ---- ---- ----- ----- ----- ',
     $       '----- ---- ---- ---- ----')
 7000 FORMAT(I6,1X,A7,5X,A1,3X,A1,2(4X,A1), 4I6, 4I5)
 8000 FORMAT('COMPLEX BSBR TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('COMPLEX BSBR TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of CBSBRTEST.
*
      END
*
*
      SUBROUTINE ZBSBRTEST( OUTNUM, VERB, NSCOPE, SCOPE0, NTOP, TOP0,
     $                      NSHAPE, UPLO0, DIAG0, NMAT, M0, N0, LDAS0,
     $                      LDAD0, NSRC, RSRC0, CSRC0, NGRID, CONTEXT0,
     $                      P0, Q0, TFAIL, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER OUTNUM, VERB, NSCOPE, NTOP, NSHAPE, NMAT, NSRC, NGRID
      INTEGER MEMLEN
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      CHARACTER*1 UPLO0(NSHAPE), DIAG0(NSHAPE)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RSRC0(NSRC), CSRC0(NSRC), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), TFAIL(*)
      DOUBLE COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ZTESTBSBR:  Test double complex broadcast
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NSHAPE   (input) INTEGER
*           The number of matrix shapes to be tested.
*
*  UPLO0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of UPLO to be tested.
*
*  DIAG0    (input) CHARACTER*1 array of dimension (NSHAPE)
*           Values of DIAG to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NSRC     (input) INTEGER
*           The number of sources to be tested.
*
*  RSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of RSRC (row coordinate of source) to be tested.
*
*  CSRC0    (input) INTEGER array of dimension (NDEST)
*           Values of CSRC (column coordinate of source) to be tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  TFAIL    (workspace) INTEGER array of dimension (NTESTS)
*           If VERB < 2, serves to indicate which tests fail.  This
*           requires workspace of NTESTS (number of tests performed).
*
*  MEM      (workspace) DOUBLE COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO
      EXTERNAL ZTRBS2D, ZGEBS2D, ZTRBR2D, ZGEBR2D
      EXTERNAL ZINITMAT, ZCHKMAT, ZCHKPAD, ZBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP, UPLO, DIAG
      LOGICAL TESTOK, INGRID
      INTEGER IAM, I, K, J, IGR, ISH, IMA, ISO, ISC, ITO
      INTEGER M, N, NPROW, NPCOL, MYROW, MYCOL, RSRC, CSRC
      INTEGER ISTART, ISTOP, IPRE, IPOST, SETWHAT
      INTEGER NERR, NSKIP, NFAIL, TESTNUM, CONTEXT, MAXERR, LDASRC
      INTEGER LDADST, ERRDPTR, APTR, ERRIPTR, ISIZE, ZSIZE
      DOUBLE COMPLEX SCHECKVAL, RCHECKVAL
*     ..
*     .. Executable Statements ..
*
      SCHECKVAL = DCMPLX( -0.01D0, -0.01D0 )
      RCHECKVAL = DCMPLX( -0.02D0, -0.02D0 )
*
      IAM = IBTMYPROC()
      ISIZE = IBTSIZEOF('I')
      ZSIZE = IBTSIZEOF('Z')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NSHAPE:', NSHAPE
            WRITE(OUTNUM, 3000) ' UPLO :', ( UPLO0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 3000) ' DIAG :', ( DIAG0(I), I = 1, NSHAPE )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NSRC  :', NSRC
            WRITE(OUTNUM, 2000) ' RSRC :',( RSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) ' CSRC :',( CSRC0(I), I = 1, NSRC )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,5000)
            WRITE(OUTNUM,6000)
         END IF
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + 4 * M0(IMA)
         IF( K .GT. I ) I = K
   10 CONTINUE
      MAXERR = ( ZSIZE * (MEMLEN-I) ) / ( ZSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run BSBR tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 110 IGR = 1, NGRID
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
*
         INGRID = ( NPROW .GT. 0 )
*
         DO 100 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 90 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multipath ('M') or general tree ('T'),
*              need to loop over calls to BLACS_SET
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 11
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 12
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 80 ISH = 1, NSHAPE
                  UPLO = UPLO0(ISH)
                  DIAG = DIAG0(ISH)
*
                  DO 70 IMA = 1, NMAT
                     M = M0(IMA)
                     N = N0(IMA)
                     LDASRC = LDAS0(IMA)
                     LDADST = LDAD0(IMA)
*
                     DO 60 ISO = 1, NSRC
                        TESTNUM = TESTNUM + 1
                        RSRC = RSRC0(ISO)
                        CSRC = CSRC0(ISO)
                        IF( RSRC.GE.P0(IGR) .OR. CSRC.GE.Q0(IGR) ) THEN
                           NSKIP = NSKIP + 1
                           GOTO 60
                        END IF
                        IF( VERB .GT. 1 ) THEN
                           IF( IAM .EQ. 0 ) THEN
                              WRITE(OUTNUM, 7000)
     $                        TESTNUM, 'RUNNING',SCOPE, TOP, UPLO, DIAG,
     $                        M, N, LDASRC, LDADST, RSRC, CSRC,
     $                        NPROW, NPCOL
                           END IF
                        END IF
*
                        TESTOK = .TRUE.
                        IPRE  = 2 * M
                        IPOST = IPRE
                        APTR = IPRE + 1
*
*                       If I am in scope
*
                        IF( (MYROW.EQ.RSRC .AND. SCOPE.EQ.'R') .OR.
     $                       (MYCOL.EQ.CSRC .AND. SCOPE.EQ.'C') .OR.
     $                       (SCOPE .EQ. 'A') ) THEN
*
*                          source process generates matrix and sends it
*
                           IF( MYROW.EQ.RSRC .AND. MYCOL.EQ.CSRC ) THEN
                              CALL ZINITMAT(UPLO, DIAG, M, N, MEM,
     $                                      LDASRC, IPRE, IPOST,
     $                                      SCHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              DO 20 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 20
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                     CALL ZTRBS2D(CONTEXT, SCOPE, TOP,
     $                                            UPLO, DIAG, M, N,
     $                                            MEM(APTR), LDASRC )
                                 ELSE
                                     CALL ZGEBS2D(CONTEXT, SCOPE, TOP,
     $                                            M, N, MEM(APTR),
     $                                            LDASRC )
                                 END IF
   20                         CONTINUE
*
*                          Destination processes
*
                           ELSE IF( INGRID ) THEN
                              DO 40 J = ISTART, ISTOP
                                 IF( J.EQ.0 ) GOTO 40
                                 IF( SETWHAT.NE.0 )
     $                              CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*                                Pad entire matrix area
*
                                 DO 30 K = 1, IPRE+IPOST+LDADST*N
                                    MEM(K) = RCHECKVAL
   30                            CONTINUE
*
*                                Receive matrix
*
                                 IF( UPLO.EQ.'U' .OR. UPLO.EQ.'L' ) THEN
                                    CALL ZTRBR2D(CONTEXT, SCOPE, TOP,
     $                                           UPLO, DIAG, M, N,
     $                                           MEM(APTR), LDADST,
     $                                           RSRC, CSRC)
                                 ELSE
                                    CALL ZGEBR2D(CONTEXT, SCOPE, TOP,
     $                                           M, N, MEM(APTR),
     $                                           LDADST, RSRC, CSRC)
                                 END IF
*
*                                Check for errors in matrix or padding
*
                                 I = NERR
                                 CALL ZCHKMAT(UPLO, DIAG, M, N,
     $                                   MEM(APTR), LDADST, RSRC, CSRC,
     $                                   MYROW, MYCOL, TESTNUM, MAXERR,
     $                                   NERR, MEM(ERRIPTR),
     $                                   MEM(ERRDPTR))
*
                                 CALL ZCHKPAD(UPLO, DIAG, M, N, MEM,
     $                                   LDADST, RSRC, CSRC, MYROW,
     $                                   MYCOL, IPRE, IPOST, RCHECKVAL,
     $                                   TESTNUM, MAXERR, NERR,
     $                                   MEM(ERRIPTR), MEM(ERRDPTR))
   40                         CONTINUE
                              TESTOK = ( I .EQ. NERR )
                           END IF
                        END IF
*
                        IF( VERB .GT. 1 ) THEN
                           I = NERR
                           CALL ZBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                                     MEM(ERRIPTR), MEM(ERRDPTR),
     $                                     TFAIL)
                           IF( IAM .EQ. 0 ) THEN
                              TESTOK = ( TESTOK .AND. (I.EQ.NERR) )
                              IF( TESTOK ) THEN
                                 WRITE(OUTNUM,7000)TESTNUM,'PASSED ',
     $                                 SCOPE, TOP, UPLO, DIAG, M, N,
     $                                 LDASRC, LDADST, RSRC, CSRC,
     $                                 NPROW, NPCOL
                              ELSE
                                 NFAIL = NFAIL + 1
                                 WRITE(OUTNUM,7000)TESTNUM,'FAILED ',
     $                                SCOPE, TOP, UPLO, DIAG, M, N,
     $                                LDASRC, LDADST, RSRC, CSRC,
     $                                NPROW, NPCOL
                              END IF
                           END IF
*
*                          Once we've printed out errors, can re-use buf space
*
                           NERR = 0
                        END IF
   60                CONTINUE
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL ZBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), TFAIL )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 8000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 9000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE COMPLEX BSBR TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 5000 FORMAT(' TEST#  STATUS SCOPE TOP UPLO DIAG     M     N  LDAS ',
     $       ' LDAD RSRC CSRC    P    Q')
 6000 FORMAT(' ----- ------- ----- --- ---- ---- ----- ----- ----- ',
     $       '----- ---- ---- ---- ----')
 7000 FORMAT(I6,1X,A7,5X,A1,3X,A1,2(4X,A1), 4I6, 4I5)
 8000 FORMAT('DOUBLE COMPLEX BSBR TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 9000 FORMAT('DOUBLE COMPLEX BSBR TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ZBSBRTEST.
*
      END
*
*
      SUBROUTINE RDCOMB( MEMUSED, MEM, MEMLEN, CMEMUSED, CMEM, CMEMLEN,
     $                   OUTNUM )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMUSED, MEMLEN, CMEMUSED, CMEMLEN, OUTNUM
*     ..
*     .. Array Arguments ..
      CHARACTER*1 CMEM(CMEMLEN)
      INTEGER MEM(MEMLEN)
*     ..
*
*     Purpose
*     =======
*     RDCOMB:  Read and process the input file COMB.dat.
*
*     Arguments
*     =========
*     MEMUSED  (output) INTEGER
*              Number of elements in MEM that this subroutine ends up using.
*
*     MEM      (output) INTEGER array of dimension memlen
*              On output, holds information read in from sdrv.dat.
*
*     MEMLEN   (input) INTEGER
*              Number of elements of MEM that this subroutine
*              may safely write into.
*
*     CMEMUSED (output) INTEGER
*              Number of elements in CMEM that this subroutine ends up using.
*
*     CMEM     (output) CHARACTER*1 array of dimension cmemlen
*              On output, holds the values for UPLO and DIAG.
*
*     CMEMLEN  (input) INTEGER
*              Number of elements of CMEM that this subroutine
*              may safely write into.
*
*     OUTNUM   (input) INTEGER
*              Unit number of the output file.
*
*     =================================================================
*
*     .. Parameters ..
      INTEGER SDIN
      PARAMETER( SDIN = 12 )
*     ..
*     .. External Functions ..
      LOGICAL  LSAME
      EXTERNAL LSAME
*     ..
*     .. Local Scalars ..
      INTEGER TOPSREPEAT, TOPSCOHRNT, NOPS, NSCOPE, NTOP, NMAT, NDEST
      INTEGER NGRID, I, J, OPPTR, SCOPEPTR, TOPPTR, MPTR, NPTR
      INTEGER LDSPTR, LDDPTR, LDIPTR, RDESTPTR, CDESTPTR, PPTR, QPTR
*     ..
*     .. Executable Statements
*
*     Open and read the file comb.dat.  The expected format is
*     below.
*
*------
*integer                         Number of operations
*array of CHAR*1's               OPs: '+', '>', '<'
*integer                         Number of scopes
*array of CHAR*1's               Values for Scopes
*HAR*1                           Repeatability flag ('R', 'N', 'B')
*HAR*1                           Coherency flag ('C', 'N', 'B')
*integer                         Number of topologies
*array of CHAR*1's               Values for TOP
*integer                         number of nmat
*array of integers               M: number of rows in matrix
*array of integers               N: number of columns in matrix
*integer                         LDA: leading dimension on source proc
*integer                         LDA: leading dimension on dest proc
*integer                         number of source/dest pairs
*array of integers               RDEST: process row of msg. dest.
*array of integers               CDEST: process column of msg. dest.
*integer                         Number of grids
*array of integers               NPROW: number of rows in process grid
*array of integers               NPCOL: number of col's in proc. grid
*------
*  note: the text descriptions as shown above are present in
*             the sample comb.dat included with this distribution,
*             but are not required.
*
*     Read input file
*
      MEMUSED = 1
      CMEMUSED = 1
      OPEN(UNIT = SDIN, FILE = 'comb.dat', STATUS = 'OLD')
*
*     Get what operations to test (+, >, <)
*
      READ(SDIN, *) NOPS
      OPPTR = CMEMUSED
      CMEMUSED = OPPTR + NOPS
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NOPS, 'OPERATIONS.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NOPS .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'OPERATIONS.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
*
      READ(SDIN, *) ( CMEM(OPPTR+I), I = 0, NOPS-1 )
      DO 10 I = 0, NOPS-1
         IF( (CMEM(OPPTR+I).NE.'+') .AND. (CMEM(OPPTR+I).NE.'>') .AND.
     $       (CMEM(OPPTR+I).NE.'<') ) THEN
            WRITE(OUTNUM,5000) CMEM(OPPTR+I)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   10 CONTINUE
*
*     Read in scopes and topologies
*
      READ(SDIN, *) NSCOPE
      SCOPEPTR = CMEMUSED
      CMEMUSED = SCOPEPTR + NSCOPE
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NSCOPE, 'SCOPES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NSCOPE .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'SCOPE.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
*
      READ(SDIN, *) ( CMEM(SCOPEPTR+I), I = 0, NSCOPE-1 )
      DO 20 I = 0, NSCOPE-1
         IF( LSAME(CMEM(SCOPEPTR+I), 'R') ) THEN
            CMEM(SCOPEPTR+I) = 'R'
         ELSE IF( LSAME(CMEM(SCOPEPTR+I), 'C') ) THEN
            CMEM(SCOPEPTR+I) = 'C'
         ELSE IF( LSAME(CMEM(SCOPEPTR+I), 'A') ) THEN
            CMEM(SCOPEPTR+I) = 'A'
         ELSE
            WRITE(OUTNUM, 3000) 'SCOPE', CMEM(SCOPEPTR+I)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   20 CONTINUE
*
      READ(SDIN, *) TOPSREPEAT
      READ(SDIN, *) TOPSCOHRNT
*
      READ(SDIN, *) NTOP
      TOPPTR = CMEMUSED
      CMEMUSED = TOPPTR + NTOP
      IF ( CMEMUSED .GT. CMEMLEN ) THEN
         WRITE(OUTNUM, 1000) CMEMLEN, NTOP, 'TOPOLOGIES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NTOP .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'TOPOLOGY.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( CMEM(TOPPTR+I), I = 0, NTOP-1 )
*
*
*     Read in number of matrices, and values for M, N, LDASRC, and LDADEST
*
      READ(SDIN, *) NMAT
      MPTR = MEMUSED
      NPTR = MPTR + NMAT
      LDSPTR = NPTR + NMAT
      LDDPTR = LDSPTR + NMAT
      LDIPTR = LDDPTR + NMAT
      MEMUSED = LDIPTR + NMAT
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'MATRICES.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NMAT .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'MATRIX.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM( MPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( NPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDSPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDDPTR+I ), I = 0, NMAT-1 )
      READ(SDIN, *) ( MEM( LDIPTR+I ), I = 0, NMAT-1 )
*
*     Make sure matrix values are legal
*
      CALL CHKMATDAT( OUTNUM, 'COMB.dat', .TRUE., NMAT, MEM(MPTR),
     $                MEM(NPTR), MEM(LDSPTR), MEM(LDDPTR), MEM(LDIPTR) )
*
*     Read in number of dest pairs, and values of dest
*
      READ(SDIN, *) NDEST
      RDESTPTR  = MEMUSED
      CDESTPTR  = RDESTPTR  + NDEST
      MEMUSED  = CDESTPTR + NDEST
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NMAT, 'DEST.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NDEST .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'DEST.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      END IF
      READ(SDIN, *) ( MEM(RDESTPTR+I), I = 0, NDEST-1 )
      READ(SDIN, *) ( MEM(CDESTPTR+I), I = 0, NDEST-1 )
*
*     Read in number of grids pairs, and values of P (process rows) and
*     Q (process columns)
*
      READ(SDIN, *) NGRID
      PPTR = MEMUSED
      QPTR = PPTR + NGRID
      MEMUSED = QPTR + NGRID
      IF( MEMUSED .GT. MEMLEN ) THEN
         WRITE(OUTNUM, 1000) MEMLEN, NGRID, 'PROCESS GRIDS.'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
         STOP
      ELSE IF( NGRID .LT. 1 ) THEN
         WRITE(OUTNUM, 2000) 'PROCESS GRID'
         IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE( OUTNUM )
         STOP
      END IF
*
      READ(SDIN, *) ( MEM(PPTR+I), I = 0, NGRID-1 )
      READ(SDIN, *) ( MEM(QPTR+I), I = 0, NGRID-1 )
      IF( SDIN .NE. 6 .AND. SDIN .NE. 0 ) CLOSE( SDIN )
*
*     Fatal error if we've got an illegal grid
*
      DO 70 J = 0, NGRID-1
         IF( MEM(PPTR+J).LT.1 .OR. MEM(QPTR+J).LT.1 ) THEN
            WRITE(OUTNUM, 4000) MEM(PPTR+J), MEM(QPTR+J)
            IF( OUTNUM .NE. 6 .AND. OUTNUM .NE. 0 ) CLOSE(OUTNUM)
            STOP
         END IF
   70 CONTINUE
*
*     Prepare output variables
*
      MEM(MEMUSED)   = NOPS
      MEM(MEMUSED+1) = NSCOPE
      MEM(MEMUSED+2) = TOPSREPEAT
      MEM(MEMUSED+3) = TOPSCOHRNT
      MEM(MEMUSED+4) = NTOP
      MEM(MEMUSED+5) = NMAT
      MEM(MEMUSED+6) = NDEST
      MEM(MEMUSED+7) = NGRID
      MEMUSED = MEMUSED + 7
      CMEMUSED = CMEMUSED - 1
*
 1000 FORMAT('Mem too short (',I4,') to handle',I4,' ',A20)
 2000 FORMAT('Must have at least one ',A20)
 3000 FORMAT('UNRECOGNIZABLE ',A5,' ''', A1, '''.')
 4000 FORMAT('Illegal process grid: {',I3,',',I3,'}.')
 5000 FORMAT('Illegal OP value ''',A1,''':, expected ''+'' (SUM),',
     $       ' ''>'' (MAX), or ''<'' (MIN).')
*
      RETURN
*
*     End of RDCOMB.
*
      END
*
*
      SUBROUTINE IBTCHECKIN( NFTESTS, OUTNUM, MAXERR, NERR, IERR,
     $                       IVAL, TFAILED )
      INTEGER NFTESTS, OUTNUM, MAXERR, NERR
      INTEGER IERR(*), TFAILED(*)
      INTEGER IVAL(*)
*
*  Purpose
*  =======
*  IBTCHECKIN: Process 0 receives error report from all processes.
*
*  Arguments
*  =========
*  NFTESTS  (input/output) INTEGER
*           if NFTESTS is <= 0 upon entry, NFTESTS is not written to.
*           Otherwise, on entry it specifies the total number of tests
*           run, and on exit it is the number of tests which failed.
*
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRIBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (workspace) INTEGER array, dimension NFTESTS
*          Workspace used to keep track of which tests failed.
*          If input of NFTESTS < 1, this array not accessed.
*
*  ===================================================================
*
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTMSGID
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID
*     ..
*     .. Local Scalars ..
      LOGICAL COUNTING
      INTEGER K, NERR2, IAM, NPROCS, NTESTS
*
*     Proc 0 collects error info from everyone
*
      IAM = IBTMYPROC()
      NPROCS = IBTNPROCS()
*
      IF( IAM .EQ. 0 ) THEN
*
*        If we are finding out how many failed tests there are, initialize
*        the total number of tests (NTESTS), and zero the test failed array
*
         COUNTING = NFTESTS .GT. 0
         IF( COUNTING ) THEN
            NTESTS = NFTESTS
            DO 10 K = 1, NTESTS
               TFAILED(K) = 0
   10       CONTINUE
         END IF
*
         CALL IPRINTERRS(OUTNUM, MAXERR, NERR, IERR, IVAL, COUNTING,
     $                   TFAILED)
*
         DO 20 K = 1, NPROCS-1
            CALL BTSEND(3, 0, K, K, IBTMSGID()+50)
            CALL BTRECV(3, 1, NERR2, K, IBTMSGID()+50)
            IF( NERR2 .GT. 0 ) THEN
               NERR = NERR + NERR2
               CALL BTRECV(3, NERR2*6, IERR, K, IBTMSGID()+51)
               CALL BTRECV(3, NERR2*2, IVAL, K, IBTMSGID()+51)
               CALL IPRINTERRS(OUTNUM, MAXERR, NERR2, IERR, IVAL,
     $                         COUNTING, TFAILED)
            END IF
   20    CONTINUE
*
*        Count up number of tests that failed
*
         IF( COUNTING ) THEN
            NFTESTS = 0
            DO 30 K = 1, NTESTS
               NFTESTS = NFTESTS + TFAILED(K)
   30       CONTINUE
         END IF
*
*     Send my error info to proc 0
*
      ELSE
         CALL BTRECV(3, 0, K, 0, IBTMSGID()+50)
         CALL BTSEND(3, 1, NERR, 0, IBTMSGID()+50)
         IF( NERR .GT. 0 ) THEN
            CALL BTSEND(3, NERR*6, IERR, 0, IBTMSGID()+51)
            CALL BTSEND(3, NERR*2, IVAL, 0, IBTMSGID()+51)
         END IF
      ENDIF
*
      RETURN
*
*     End of IBTCHECKIN
*
      END
*
      SUBROUTINE IINITMAT(UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL, TESTNUM, MYROW, MYCOL)
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST, TESTNUM, MYROW, MYCOL
      INTEGER CHECKVAL
      INTEGER MEM(*)
*
*     .. External Subroutines ..
      EXTERNAL IGENMAT, IPADMAT
*     ..
*     .. Executable Statements ..
*
      CALL IGENMAT( M, N, MEM(IPRE+1), LDA, TESTNUM, MYROW, MYCOL )
      CALL IPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST, CHECKVAL )
*
      RETURN
      END
*
      SUBROUTINE IGENMAT( M, N, A, LDA, TESTNUM, MYROW, MYCOL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA, TESTNUM, MYROW, MYCOL
*     ..
*     .. Array Arguments ..
      INTEGER A(LDA,N)
*     ..
*
*  Purpose
*  =======
*  IGENMAT: Generates an M-by-N matrix filled with random elements.
*
*  Arguments
*  =========
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (output) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  TESTNUM  (input) INTEGER
*           Unique number for this test case, used as a basis for
*           the random seeds.
*
*  ====================================================================
*
*     .. External Functions ..
      INTEGER IBTNPROCS
      INTEGER IBTRAN
      EXTERNAL IBTRAN, IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. Executable Statements ..
*
*     ISEED's four values must be positive integers less than 4096,
*     fourth one has to be odd. (see _LARND).  Use some goofy
*     functions to come up with seed values which together should
*     be unique.
*
      NPROCS = IBTNPROCS()
      SRC = MYROW * NPROCS + MYCOL
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
      DO 10 J = 1, N
         DO 10 I = 1, M
            A(I, J) = IBTRAN( ISEED )
   10 CONTINUE
*
      RETURN
*
*     End of IGENMAT.
*
      END
*
      INTEGER FUNCTION IBTRAN(ISEED)
      INTEGER ISEED(*)
*
*     .. External Functions ..
      DOUBLE PRECISION DLARND
      EXTERNAL DLARND
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION DVAL
*     ..
*     .. Executable Statements ..
*
      DVAL = 1.0D6 * DLARND(2, ISEED)
      IBTRAN = INT(DVAL)
*
      RETURN
*
*     End of Ibtran
*
      END
*
      SUBROUTINE IPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST
      INTEGER CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER MEM( * )
*     ..
*
*  Purpose
*  =======
*
*  IPADMAT: Pad Matrix.
*  This routines surrounds a matrix with a guardzone initialized to the
*  value CHECKVAL.  There are three distinct guardzones:
*  - A contiguous zone of size IPRE immediately before the start
*    of the matrix.
*  - A contiguous zone of size IPOST immedately after the end of the
*    matrix.
*  - Interstitial zones within each column of the matrix, in the
*    elements A( M+1:LDA, J ).
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM      (output) integer array, dimension (IPRE+IPOST+LDA*N)
*           The address IPRE elements ahead of the matrix A you want to
*           pad, which is then of dimension (LDA,N).
*
*  IPRE     (input) INTEGER
*           The size of the guard zone ahead of the matrix A.
*
*  IPOST    (input) INTEGER
*           The size of the guard zone behind the matrix A.
*
*  CHECKVAL (input) integer
*           The value to insert into the guard zones.
*
*  ====================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            MEM( I ) = CHECKVAL
   10    CONTINUE
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            MEM( I ) = CHECKVAL
   20    CONTINUE
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K+LDA-M-1
               MEM( I ) = CHECKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
*     If the matrix is upper or lower trapezoidal, calculate the
*     additional triangular area which needs to be padded,  Each
*     element referred to is in the Ith row and the Jth column.
*
      IF( UPLO .EQ. 'U' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 41 I = 1, M
                  DO 42 J = 1, I
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   42             CONTINUE
   41          CONTINUE
            ELSE
               DO 43 I = 2, M
                  DO 44 J = 1, I-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   44             CONTINUE
   43          CONTINUE
            END IF
         ELSE
            IF( DIAG .EQ. 'U' ) THEN
               DO 45 I = M-N+1, M
                  DO 46 J = 1, I-(M-N)
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   46             CONTINUE
   45          CONTINUE
            ELSE
               DO 47 I = M-N+2, M
                  DO 48 J = 1, I-(M-N)-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   48             CONTINUE
   47          CONTINUE
            END IF
         END IF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 49 I = 1, M
                  DO 50 J = N-M+I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   50             CONTINUE
   49          CONTINUE
            ELSE
               DO 51 I = 1, M-1
                  DO 52 J = N-M+I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   52             CONTINUE
   51          CONTINUE
            END IF
         ELSE
            IF( UPLO .EQ. 'U' ) THEN
               DO 53 I = 1, N
                  DO 54 J = I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   54             CONTINUE
   53          CONTINUE
            ELSE
               DO 55 I = 1, N-1
                  DO 56 J = I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   56             CONTINUE
   55          CONTINUE
            END IF
         END IF
      END IF
*
*     End of IPADMAT.
*
      RETURN
      END
*
      SUBROUTINE ICHKPAD( UPLO, DIAG, M, N, MEM, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, IPRE, IPOST, CHECKVAL,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST
      INTEGER TESTNUM, MAXERR, NERR
      INTEGER CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      INTEGER MEM(*), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  ICHKPAD: Check padding put in by PADMAT.
*  Checks that padding around target matrix has not been overwritten
*  by the previous point-to-point or broadcast send.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM       (input) integer array, dimension(IPRE+IPOST+LDA*N).
*            Memory location IPRE elements in front of the matrix A.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*  IPRE     (input) INTEGER
*           The size of the guard zone before the start of A.
*
*  IPOST    (input) INTEGER
*           The size of guard zone after A.
*
*  CHECKVAL (input) integer
*           The value to pad matrix with.
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRIBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL ISTRAP
      INTEGER I, J, K, IRST, IRND, ICST, ICND, SRC, DEST
      INTEGER NPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE - I + 1
                  ERRIBUF(6, NERR) = ERR_PRE
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   10    CONTINUE
      END IF
*
*     Check buffer behind A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I - J + 1
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = ERR_POST
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         DO 40 J = 1, N
            DO 30 I = M+1, LDA
               K = IPRE + (J-1)*LDA + I
               IF( MEM(K) .NE. CHECKVAL) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_GAP
                     ERRDBUF(1, NERR) = MEM(K)
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Determine limits of trapezoidal matrix
*
      ISTRAP = .FALSE.
      IF( UPLO .EQ. 'U' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 2
            IRND = M
            ICST = 1
            ICND = M - 1
         ELSEIF( M .GT. N ) THEN
            IRST = ( M-N ) + 2
            IRND = M
            ICST = 1
            ICND = N - 1
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            IRST = IRST - 1
            ICND = ICND + 1
         ENDIF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 1
            IRND = 1
            ICST = ( N-M ) + 2
            ICND = N
         ELSEIF( M .GT. N ) THEN
            IRST = 1
            IRND = 1
            ICST = 2
            ICND = N
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            ICST = ICST - 1
         ENDIF
      ENDIF
*
*     Check elements and report any errors
*
      IF( ISTRAP ) THEN
         DO 100 J = ICST, ICND
            DO 105 I = IRST, IRND
               IF( MEM( IPRE + (J-1)*LDA + I ) .NE. CHECKVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_TRI
                     ERRDBUF(1, NERR) = MEM( IPRE + (J-1)*LDA + I )
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
  105       CONTINUE
*
*           Update the limits to allow filling in padding
*
            IF( UPLO .EQ. 'U' ) THEN
               IRST = IRST + 1
            ELSE
               IRND = IRND + 1
            ENDIF
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of ICHKPAD.
*
      END
*
      SUBROUTINE ICHKMAT( UPLO, DIAG, M, N, A, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, TESTNUM, MAXERR, NERR,
     $                    ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      INTEGER A(LDA,N), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  iCHKMAT:  Check matrix to see whether there were any transmission
*            errors.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (input) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRIBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC, DEST
      LOGICAL USEIT
      INTEGER COMPVAL
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      INTEGER IBTRAN
      EXTERNAL IBTRAN, IBTNPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Initialize ISEED with the same values as used in IGENMAT.
*
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
*     Generate the elements randomly with the same method used in GENMAT.
*     Note that for trapezoidal matrices, we generate all elements in the
*     enclosing rectangle and then ignore the complementary triangle.
*
      DO 100 J = 1, N
         DO 105 I = 1, M
            COMPVAL = IBTRAN( ISEED )
*
*           Now determine whether we actually need this value.  The
*           strategy is to chop out the proper triangle based on what
*           particular kind of trapezoidal matrix we're dealing with.
*
            USEIT = .TRUE.
            IF( UPLO .EQ. 'U' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF( UPLO .EQ. 'L' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J. GE. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J .GE. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            END IF
*
*           Compare the generated value to the one that's in the
*           received matrix.  If they don't match, tack another
*           error record onto what's already there.
*
            IF( USEIT ) THEN
               IF( A(I,J) .NE. COMPVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = 5
                     ERRDBUF(1, NERR) = A(I, J)
                     ERRDBUF(2, NERR) = COMPVAL
                  END IF
               END IF
            END IF
  105    CONTINUE
  100 CONTINUE
      RETURN
*
*     End of ICHKMAT.
*
      END
*
      SUBROUTINE IPRINTERRS( OUTNUM, MAXERR, NERR,
     $                       ERRIBUF, ERRDBUF, COUNTING, TFAILED )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL COUNTING
      INTEGER OUTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), TFAILED(*)
      INTEGER ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  IPRINTERRS: Print errors that have been recorded
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRIBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (input/ourput) INTEGER array, dimension NTESTS
*          Workspace used to keep track of which tests failed.
*          This array not accessed unless COUNTING is true.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      EXTERNAL IBTMYPROC, IBTNPROCS
*     ..
*     .. Local Scalars ..
      CHARACTER*1 MAT
      LOGICAL MATISINT
      INTEGER OLDTEST, NPROCS, PROW, PCOL, I, ERRTYPE
*     ..
*     .. Executable Statements ..
*
      IF( (IBTMYPROC().NE.0) .OR. (NERR.LE.0) ) RETURN
      OLDTEST = -1
      NPROCS = IBTNPROCS()
      PROW = ERRIBUF(3,1) / NPROCS
      PCOL = MOD( ERRIBUF(3,1), NPROCS )
      IF( NERR .GT. MAXERR ) WRITE(OUTNUM,13000)
*
      DO 20 I = 1, MIN( NERR, MAXERR )
         IF( ERRIBUF(1,I) .NE. OLDTEST ) THEN
            IF( OLDTEST .NE. -1 )
     $         WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM,1000) PROW, PCOL, ERRIBUF(1,I)
            IF( COUNTING ) TFAILED( ERRIBUF(1,I) ) = 1
            OLDTEST = ERRIBUF(1, I)
         END IF
*
*        Print out error message depending on type of error
*
         ERRTYPE = ERRIBUF(6, I)
         IF( ERRTYPE .LT. -10 ) THEN
            ERRTYPE = -ERRTYPE - 10
            MAT = 'C'
            MATISINT = .TRUE.
         ELSE IF( ERRTYPE .LT. 0 ) THEN
            ERRTYPE = -ERRTYPE
            MAT = 'R'
            MATISINT = .TRUE.
         ELSE
            MATISINT = .FALSE.
         END IF
*
*        RA/CA arrays from MAX/MIN have different printing protocol
*
         IF( MATISINT ) THEN
            IF( ERRIBUF(2, I) .EQ. -1 ) THEN
               WRITE(OUTNUM,11000) ERRIBUF(4,I), ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,7000) ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,8000) ERRIBUF(4,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,9000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,10000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $                             INT( ERRDBUF(2,I) ),
     $                             INT( ERRDBUF(1,I) )
            END IF
*
*        Have memory overwrites in matrix A
*
         ELSE
            IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,2000) ERRIBUF(5,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,3000) ERRIBUF(4,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,4000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_TRI ) THEN
               WRITE(OUTNUM,5000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE
               WRITE(OUTNUM,6000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            END IF
         END IF
   20 CONTINUE
      WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
*
 1000 FORMAT('PROCESS {',I4,',',I4,'} REPORTS ERRORS IN TEST#',I6,':')
 2000 FORMAT('   Buffer overwrite ',I4,
     $       ' elements before the start of A:',/,
     $       '   Expected=',I12,
     $       '; Received=',I12)
 3000 FORMAT('   Buffer overwrite ',I4,' elements after the end of A:',
     $       /,'   Expected=',I12,
     $       '; Received=',I12)
 4000 FORMAT('   LDA-M gap overwrite at postion (',I4,',',I4,'):',/,
     $       '   Expected=',I12,
     $       '; Received=',I12)
 5000 FORMAT('   Complementory triangle overwrite at A(',I4,',',I4,
     $       '):',/,'   Expected=',I12,
     $       '; Received=',I12)
 6000 FORMAT('   Invalid element at A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,
     $       '; Received=',I12)
 7000 FORMAT('   Buffer overwrite ',I4,' elements before the start of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
 8000 FORMAT('   Buffer overwrite ',I4,' elements after the end of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
*
 9000 FORMAT('   LD',A1,'A-M gap overwrite at postion (',I4,',',I4,'):'
     $       ,/,'   Expected=',I12,'; Received=',I12)
*
10000 FORMAT('   Invalid element at ',A1,'A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,'; Received=',I12)
11000 FORMAT('   Overwrite at position (',I4,',',I4,') of non-existent '
     $       ,A1,'A array.',/,'   Expected=',I12,'; Received=',I12)
12000 FORMAT('PROCESS {',I4,',',I4,'} DONE ERROR REPORT FOR TEST#',
     $       I6,'.')
13000 FORMAT('WARNING: There were more errors than could be recorded.',
     $       /,'Increase MEMELTS to get complete listing.')
      RETURN
*
*     End IPRINTERRS
*
      END
*
*
      SUBROUTINE SBTCHECKIN( NFTESTS, OUTNUM, MAXERR, NERR, IERR,
     $                       SVAL, TFAILED )
      INTEGER NFTESTS, OUTNUM, MAXERR, NERR
      INTEGER IERR(*), TFAILED(*)
      REAL SVAL(*)
*
*  Purpose
*  =======
*  SBTCHECKIN: Process 0 receives error report from all processes.
*
*  Arguments
*  =========
*  NFTESTS  (input/output) INTEGER
*           if NFTESTS is <= 0 upon entry, NFTESTS is not written to.
*           Otherwise, on entry it specifies the total number of tests
*           run, and on exit it is the number of tests which failed.
*
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRSBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (workspace) INTEGER array, dimension NFTESTS
*          Workspace used to keep track of which tests failed.
*          If input of NFTESTS < 1, this array not accessed.
*
*  ===================================================================
*
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTMSGID
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID
*     ..
*     .. Local Scalars ..
      LOGICAL COUNTING
      INTEGER K, NERR2, IAM, NPROCS, NTESTS
*
*     Proc 0 collects error info from everyone
*
      IAM = IBTMYPROC()
      NPROCS = IBTNPROCS()
*
      IF( IAM .EQ. 0 ) THEN
*
*        If we are finding out how many failed tests there are, initialize
*        the total number of tests (NTESTS), and zero the test failed array
*
         COUNTING = NFTESTS .GT. 0
         IF( COUNTING ) THEN
            NTESTS = NFTESTS
            DO 10 K = 1, NTESTS
               TFAILED(K) = 0
   10       CONTINUE
         END IF
*
         CALL SPRINTERRS(OUTNUM, MAXERR, NERR, IERR, SVAL, COUNTING,
     $                   TFAILED)
*
         DO 20 K = 1, NPROCS-1
            CALL BTSEND(3, 0, K, K, IBTMSGID()+50)
            CALL BTRECV(3, 1, NERR2, K, IBTMSGID()+50)
            IF( NERR2 .GT. 0 ) THEN
               NERR = NERR + NERR2
               CALL BTRECV(3, NERR2*6, IERR, K, IBTMSGID()+51)
               CALL BTRECV(4, NERR2*2, SVAL, K, IBTMSGID()+51)
               CALL SPRINTERRS(OUTNUM, MAXERR, NERR2, IERR, SVAL,
     $                         COUNTING, TFAILED)
            END IF
   20    CONTINUE
*
*        Count up number of tests that failed
*
         IF( COUNTING ) THEN
            NFTESTS = 0
            DO 30 K = 1, NTESTS
               NFTESTS = NFTESTS + TFAILED(K)
   30       CONTINUE
         END IF
*
*     Send my error info to proc 0
*
      ELSE
         CALL BTRECV(3, 0, K, 0, IBTMSGID()+50)
         CALL BTSEND(3, 1, NERR, 0, IBTMSGID()+50)
         IF( NERR .GT. 0 ) THEN
            CALL BTSEND(3, NERR*6, IERR, 0, IBTMSGID()+51)
            CALL BTSEND(4, NERR*2, SVAL, 0, IBTMSGID()+51)
         END IF
      ENDIF
*
      RETURN
*
*     End of SBTCHECKIN
*
      END
*
      SUBROUTINE SINITMAT(UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL, TESTNUM, MYROW, MYCOL)
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST, TESTNUM, MYROW, MYCOL
      REAL CHECKVAL
      REAL MEM(*)
*
*     .. External Subroutines ..
      EXTERNAL SGENMAT, SPADMAT
*     ..
*     .. Executable Statements ..
*
      CALL SGENMAT( M, N, MEM(IPRE+1), LDA, TESTNUM, MYROW, MYCOL )
      CALL SPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST, CHECKVAL )
*
      RETURN
      END
*
      SUBROUTINE SGENMAT( M, N, A, LDA, TESTNUM, MYROW, MYCOL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA, TESTNUM, MYROW, MYCOL
*     ..
*     .. Array Arguments ..
      REAL A(LDA,N)
*     ..
*
*  Purpose
*  =======
*  SGENMAT: Generates an M-by-N matrix filled with random elements.
*
*  Arguments
*  =========
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (output) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  TESTNUM  (input) INTEGER
*           Unique number for this test case, used as a basis for
*           the random seeds.
*
*  ====================================================================
*
*     .. External Functions ..
      INTEGER IBTNPROCS
      REAL SBTRAN
      EXTERNAL SBTRAN, IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. Executable Statements ..
*
*     ISEED's four values must be positive integers less than 4096,
*     fourth one has to be odd. (see _LARND).  Use some goofy
*     functions to come up with seed values which together should
*     be unique.
*
      NPROCS = IBTNPROCS()
      SRC = MYROW * NPROCS + MYCOL
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
      DO 10 J = 1, N
         DO 10 I = 1, M
            A(I, J) = SBTRAN( ISEED )
   10 CONTINUE
*
      RETURN
*
*     End of SGENMAT.
*
      END
*
      REAL FUNCTION SBTRAN(ISEED)
      INTEGER ISEED(*)
*
*     .. External Functions ..
      DOUBLE PRECISION DLARND
      EXTERNAL DLARND
*     .. Executable Statements ..
*
      SBTRAN = REAL( DLARND(2, ISEED) )
*
      RETURN
*
*     End of Sbtran
*
      END
*
      SUBROUTINE SPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST
      REAL CHECKVAL
*     ..
*     .. Array Arguments ..
      REAL MEM( * )
*     ..
*
*  Purpose
*  =======
*
*  SPADMAT: Pad Matrix.
*  This routines surrounds a matrix with a guardzone initialized to the
*  value CHECKVAL.  There are three distinct guardzones:
*  - A contiguous zone of size IPRE immediately before the start
*    of the matrix.
*  - A contiguous zone of size IPOST immedately after the end of the
*    matrix.
*  - Interstitial zones within each column of the matrix, in the
*    elements A( M+1:LDA, J ).
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM      (output) real array, dimension (IPRE+IPOST+LDA*N)
*           The address IPRE elements ahead of the matrix A you want to
*           pad, which is then of dimension (LDA,N).
*
*  IPRE     (input) INTEGER
*           The size of the guard zone ahead of the matrix A.
*
*  IPOST    (input) INTEGER
*           The size of the guard zone behind the matrix A.
*
*  CHECKVAL (input) real
*           The value to insert into the guard zones.
*
*  ====================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            MEM( I ) = CHECKVAL
   10    CONTINUE
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            MEM( I ) = CHECKVAL
   20    CONTINUE
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K+LDA-M-1
               MEM( I ) = CHECKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
*     If the matrix is upper or lower trapezoidal, calculate the
*     additional triangular area which needs to be padded,  Each
*     element referred to is in the Ith row and the Jth column.
*
      IF( UPLO .EQ. 'U' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 41 I = 1, M
                  DO 42 J = 1, I
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   42             CONTINUE
   41          CONTINUE
            ELSE
               DO 43 I = 2, M
                  DO 44 J = 1, I-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   44             CONTINUE
   43          CONTINUE
            END IF
         ELSE
            IF( DIAG .EQ. 'U' ) THEN
               DO 45 I = M-N+1, M
                  DO 46 J = 1, I-(M-N)
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   46             CONTINUE
   45          CONTINUE
            ELSE
               DO 47 I = M-N+2, M
                  DO 48 J = 1, I-(M-N)-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   48             CONTINUE
   47          CONTINUE
            END IF
         END IF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 49 I = 1, M
                  DO 50 J = N-M+I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   50             CONTINUE
   49          CONTINUE
            ELSE
               DO 51 I = 1, M-1
                  DO 52 J = N-M+I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   52             CONTINUE
   51          CONTINUE
            END IF
         ELSE
            IF( UPLO .EQ. 'U' ) THEN
               DO 53 I = 1, N
                  DO 54 J = I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   54             CONTINUE
   53          CONTINUE
            ELSE
               DO 55 I = 1, N-1
                  DO 56 J = I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   56             CONTINUE
   55          CONTINUE
            END IF
         END IF
      END IF
*
*     End of SPADMAT.
*
      RETURN
      END
*
      SUBROUTINE SCHKPAD( UPLO, DIAG, M, N, MEM, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, IPRE, IPOST, CHECKVAL,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST
      INTEGER TESTNUM, MAXERR, NERR
      REAL CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      REAL MEM(*), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  SCHKPAD: Check padding put in by PADMAT.
*  Checks that padding around target matrix has not been overwritten
*  by the previous point-to-point or broadcast send.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM       (input) real array, dimension(IPRE+IPOST+LDA*N).
*            Memory location IPRE elements in front of the matrix A.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*  IPRE     (input) INTEGER
*           The size of the guard zone before the start of A.
*
*  IPOST    (input) INTEGER
*           The size of guard zone after A.
*
*  CHECKVAL (input) real
*           The value to pad matrix with.
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRSBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL ISTRAP
      INTEGER I, J, K, IRST, IRND, ICST, ICND, SRC, DEST
      INTEGER NPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE - I + 1
                  ERRIBUF(6, NERR) = ERR_PRE
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   10    CONTINUE
      END IF
*
*     Check buffer behind A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I - J + 1
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = ERR_POST
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         DO 40 J = 1, N
            DO 30 I = M+1, LDA
               K = IPRE + (J-1)*LDA + I
               IF( MEM(K) .NE. CHECKVAL) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_GAP
                     ERRDBUF(1, NERR) = MEM(K)
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Determine limits of trapezoidal matrix
*
      ISTRAP = .FALSE.
      IF( UPLO .EQ. 'U' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 2
            IRND = M
            ICST = 1
            ICND = M - 1
         ELSEIF( M .GT. N ) THEN
            IRST = ( M-N ) + 2
            IRND = M
            ICST = 1
            ICND = N - 1
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            IRST = IRST - 1
            ICND = ICND + 1
         ENDIF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 1
            IRND = 1
            ICST = ( N-M ) + 2
            ICND = N
         ELSEIF( M .GT. N ) THEN
            IRST = 1
            IRND = 1
            ICST = 2
            ICND = N
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            ICST = ICST - 1
         ENDIF
      ENDIF
*
*     Check elements and report any errors
*
      IF( ISTRAP ) THEN
         DO 100 J = ICST, ICND
            DO 105 I = IRST, IRND
               IF( MEM( IPRE + (J-1)*LDA + I ) .NE. CHECKVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_TRI
                     ERRDBUF(1, NERR) = MEM( IPRE + (J-1)*LDA + I )
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
  105       CONTINUE
*
*           Update the limits to allow filling in padding
*
            IF( UPLO .EQ. 'U' ) THEN
               IRST = IRST + 1
            ELSE
               IRND = IRND + 1
            ENDIF
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of SCHKPAD.
*
      END
*
      SUBROUTINE SCHKMAT( UPLO, DIAG, M, N, A, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, TESTNUM, MAXERR, NERR,
     $                    ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      REAL A(LDA,N), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  sCHKMAT:  Check matrix to see whether there were any transmission
*            errors.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (input) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRSBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC, DEST
      LOGICAL USEIT
      REAL COMPVAL
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      REAL SBTRAN
      EXTERNAL SBTRAN, IBTNPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Initialize ISEED with the same values as used in SGENMAT.
*
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
*     Generate the elements randomly with the same method used in GENMAT.
*     Note that for trapezoidal matrices, we generate all elements in the
*     enclosing rectangle and then ignore the complementary triangle.
*
      DO 100 J = 1, N
         DO 105 I = 1, M
            COMPVAL = SBTRAN( ISEED )
*
*           Now determine whether we actually need this value.  The
*           strategy is to chop out the proper triangle based on what
*           particular kind of trapezoidal matrix we're dealing with.
*
            USEIT = .TRUE.
            IF( UPLO .EQ. 'U' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF( UPLO .EQ. 'L' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J. GE. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J .GE. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            END IF
*
*           Compare the generated value to the one that's in the
*           received matrix.  If they don't match, tack another
*           error record onto what's already there.
*
            IF( USEIT ) THEN
               IF( A(I,J) .NE. COMPVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = 5
                     ERRDBUF(1, NERR) = A(I, J)
                     ERRDBUF(2, NERR) = COMPVAL
                  END IF
               END IF
            END IF
  105    CONTINUE
  100 CONTINUE
      RETURN
*
*     End of SCHKMAT.
*
      END
*
      SUBROUTINE SPRINTERRS( OUTNUM, MAXERR, NERR,
     $                       ERRIBUF, ERRDBUF, COUNTING, TFAILED )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL COUNTING
      INTEGER OUTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), TFAILED(*)
      REAL ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  SPRINTERRS: Print errors that have been recorded
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRSBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (input/ourput) INTEGER array, dimension NTESTS
*          Workspace used to keep track of which tests failed.
*          This array not accessed unless COUNTING is true.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      EXTERNAL IBTMYPROC, IBTNPROCS
*     ..
*     .. Local Scalars ..
      CHARACTER*1 MAT
      LOGICAL MATISINT
      INTEGER OLDTEST, NPROCS, PROW, PCOL, I, ERRTYPE
*     ..
*     .. Executable Statements ..
*
      IF( (IBTMYPROC().NE.0) .OR. (NERR.LE.0) ) RETURN
      OLDTEST = -1
      NPROCS = IBTNPROCS()
      PROW = ERRIBUF(3,1) / NPROCS
      PCOL = MOD( ERRIBUF(3,1), NPROCS )
      IF( NERR .GT. MAXERR ) WRITE(OUTNUM,13000)
*
      DO 20 I = 1, MIN( NERR, MAXERR )
         IF( ERRIBUF(1,I) .NE. OLDTEST ) THEN
            IF( OLDTEST .NE. -1 )
     $         WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM,1000) PROW, PCOL, ERRIBUF(1,I)
            IF( COUNTING ) TFAILED( ERRIBUF(1,I) ) = 1
            OLDTEST = ERRIBUF(1, I)
         END IF
*
*        Print out error message depending on type of error
*
         ERRTYPE = ERRIBUF(6, I)
         IF( ERRTYPE .LT. -10 ) THEN
            ERRTYPE = -ERRTYPE - 10
            MAT = 'C'
            MATISINT = .TRUE.
         ELSE IF( ERRTYPE .LT. 0 ) THEN
            ERRTYPE = -ERRTYPE
            MAT = 'R'
            MATISINT = .TRUE.
         ELSE
            MATISINT = .FALSE.
         END IF
*
*        RA/CA arrays from MAX/MIN have different printing protocol
*
         IF( MATISINT ) THEN
            IF( ERRIBUF(2, I) .EQ. -1 ) THEN
               WRITE(OUTNUM,11000) ERRIBUF(4,I), ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,7000) ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,8000) ERRIBUF(4,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,9000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,10000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $                             INT( ERRDBUF(2,I) ),
     $                             INT( ERRDBUF(1,I) )
            END IF
*
*        Have memory overwrites in matrix A
*
         ELSE
            IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,2000) ERRIBUF(5,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,3000) ERRIBUF(4,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,4000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_TRI ) THEN
               WRITE(OUTNUM,5000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE
               WRITE(OUTNUM,6000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            END IF
         END IF
   20 CONTINUE
      WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
*
 1000 FORMAT('PROCESS {',I4,',',I4,'} REPORTS ERRORS IN TEST#',I6,':')
 2000 FORMAT('   Buffer overwrite ',I4,
     $       ' elements before the start of A:',/,
     $       '   Expected=',G15.8,
     $       '; Received=',G15.8)
 3000 FORMAT('   Buffer overwrite ',I4,' elements after the end of A:',
     $       /,'   Expected=',G15.8,
     $       '; Received=',G15.8)
 4000 FORMAT('   LDA-M gap overwrite at postion (',I4,',',I4,'):',/,
     $       '   Expected=',G15.8,
     $       '; Received=',G15.8)
 5000 FORMAT('   Complementory triangle overwrite at A(',I4,',',I4,
     $       '):',/,'   Expected=',G15.8,
     $       '; Received=',G15.8)
 6000 FORMAT('   Invalid element at A(',I4,',',I4,'):',/,
     $       '   Expected=',G15.8,
     $       '; Received=',G15.8)
 7000 FORMAT('   Buffer overwrite ',I4,' elements before the start of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
 8000 FORMAT('   Buffer overwrite ',I4,' elements after the end of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
*
 9000 FORMAT('   LD',A1,'A-M gap overwrite at postion (',I4,',',I4,'):'
     $       ,/,'   Expected=',I12,'; Received=',I12)
*
10000 FORMAT('   Invalid element at ',A1,'A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,'; Received=',I12)
11000 FORMAT('   Overwrite at position (',I4,',',I4,') of non-existent '
     $       ,A1,'A array.',/,'   Expected=',I12,'; Received=',I12)
12000 FORMAT('PROCESS {',I4,',',I4,'} DONE ERROR REPORT FOR TEST#',
     $       I6,'.')
13000 FORMAT('WARNING: There were more errors than could be recorded.',
     $       /,'Increase MEMELTS to get complete listing.')
      RETURN
*
*     End SPRINTERRS
*
      END
*
*
      SUBROUTINE DBTCHECKIN( NFTESTS, OUTNUM, MAXERR, NERR, IERR,
     $                       DVAL, TFAILED )
      INTEGER NFTESTS, OUTNUM, MAXERR, NERR
      INTEGER IERR(*), TFAILED(*)
      DOUBLE PRECISION DVAL(*)
*
*  Purpose
*  =======
*  DBTCHECKIN: Process 0 receives error report from all processes.
*
*  Arguments
*  =========
*  NFTESTS  (input/output) INTEGER
*           if NFTESTS is <= 0 upon entry, NFTESTS is not written to.
*           Otherwise, on entry it specifies the total number of tests
*           run, and on exit it is the number of tests which failed.
*
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRDBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (workspace) INTEGER array, dimension NFTESTS
*          Workspace used to keep track of which tests failed.
*          If input of NFTESTS < 1, this array not accessed.
*
*  ===================================================================
*
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTMSGID
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID
*     ..
*     .. Local Scalars ..
      LOGICAL COUNTING
      INTEGER K, NERR2, IAM, NPROCS, NTESTS
*
*     Proc 0 collects error info from everyone
*
      IAM = IBTMYPROC()
      NPROCS = IBTNPROCS()
*
      IF( IAM .EQ. 0 ) THEN
*
*        If we are finding out how many failed tests there are, initialize
*        the total number of tests (NTESTS), and zero the test failed array
*
         COUNTING = NFTESTS .GT. 0
         IF( COUNTING ) THEN
            NTESTS = NFTESTS
            DO 10 K = 1, NTESTS
               TFAILED(K) = 0
   10       CONTINUE
         END IF
*
         CALL DPRINTERRS(OUTNUM, MAXERR, NERR, IERR, DVAL, COUNTING,
     $                   TFAILED)
*
         DO 20 K = 1, NPROCS-1
            CALL BTSEND(3, 0, K, K, IBTMSGID()+50)
            CALL BTRECV(3, 1, NERR2, K, IBTMSGID()+50)
            IF( NERR2 .GT. 0 ) THEN
               NERR = NERR + NERR2
               CALL BTRECV(3, NERR2*6, IERR, K, IBTMSGID()+51)
               CALL BTRECV(6, NERR2*2, DVAL, K, IBTMSGID()+51)
               CALL DPRINTERRS(OUTNUM, MAXERR, NERR2, IERR, DVAL,
     $                         COUNTING, TFAILED)
            END IF
   20    CONTINUE
*
*        Count up number of tests that failed
*
         IF( COUNTING ) THEN
            NFTESTS = 0
            DO 30 K = 1, NTESTS
               NFTESTS = NFTESTS + TFAILED(K)
   30       CONTINUE
         END IF
*
*     Send my error info to proc 0
*
      ELSE
         CALL BTRECV(3, 0, K, 0, IBTMSGID()+50)
         CALL BTSEND(3, 1, NERR, 0, IBTMSGID()+50)
         IF( NERR .GT. 0 ) THEN
            CALL BTSEND(3, NERR*6, IERR, 0, IBTMSGID()+51)
            CALL BTSEND(6, NERR*2, DVAL, 0, IBTMSGID()+51)
         END IF
      ENDIF
*
      RETURN
*
*     End of DBTCHECKIN
*
      END
*
      SUBROUTINE DINITMAT(UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL, TESTNUM, MYROW, MYCOL)
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST, TESTNUM, MYROW, MYCOL
      DOUBLE PRECISION CHECKVAL
      DOUBLE PRECISION MEM(*)
*
*     .. External Subroutines ..
      EXTERNAL DGENMAT, DPADMAT
*     ..
*     .. Executable Statements ..
*
      CALL DGENMAT( M, N, MEM(IPRE+1), LDA, TESTNUM, MYROW, MYCOL )
      CALL DPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST, CHECKVAL )
*
      RETURN
      END
*
      SUBROUTINE DGENMAT( M, N, A, LDA, TESTNUM, MYROW, MYCOL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA, TESTNUM, MYROW, MYCOL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,N)
*     ..
*
*  Purpose
*  =======
*  DGENMAT: Generates an M-by-N matrix filled with random elements.
*
*  Arguments
*  =========
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (output) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  TESTNUM  (input) INTEGER
*           Unique number for this test case, used as a basis for
*           the random seeds.
*
*  ====================================================================
*
*     .. External Functions ..
      INTEGER IBTNPROCS
      DOUBLE PRECISION DBTRAN
      EXTERNAL DBTRAN, IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. Executable Statements ..
*
*     ISEED's four values must be positive integers less than 4096,
*     fourth one has to be odd. (see _LARND).  Use some goofy
*     functions to come up with seed values which together should
*     be unique.
*
      NPROCS = IBTNPROCS()
      SRC = MYROW * NPROCS + MYCOL
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
      DO 10 J = 1, N
         DO 10 I = 1, M
            A(I, J) = DBTRAN( ISEED )
   10 CONTINUE
*
      RETURN
*
*     End of DGENMAT.
*
      END
*
      DOUBLE PRECISION FUNCTION DBTRAN(ISEED)
      INTEGER ISEED(*)
*
*     .. External Functions ..
      DOUBLE PRECISION DLARND
      EXTERNAL DLARND
*     .. Executable Statements ..
*
      DBTRAN = DLARND(2, ISEED)
*
      RETURN
*
*     End of Dbtran
*
      END
*
      SUBROUTINE DPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST
      DOUBLE PRECISION CHECKVAL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION MEM( * )
*     ..
*
*  Purpose
*  =======
*
*  DPADMAT: Pad Matrix.
*  This routines surrounds a matrix with a guardzone initialized to the
*  value CHECKVAL.  There are three distinct guardzones:
*  - A contiguous zone of size IPRE immediately before the start
*    of the matrix.
*  - A contiguous zone of size IPOST immedately after the end of the
*    matrix.
*  - Interstitial zones within each column of the matrix, in the
*    elements A( M+1:LDA, J ).
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM      (output) double precision array, dimension (IPRE+IPOST+LDA*N)
*           The address IPRE elements ahead of the matrix A you want to
*           pad, which is then of dimension (LDA,N).
*
*  IPRE     (input) INTEGER
*           The size of the guard zone ahead of the matrix A.
*
*  IPOST    (input) INTEGER
*           The size of the guard zone behind the matrix A.
*
*  CHECKVAL (input) double precision
*           The value to insert into the guard zones.
*
*  ====================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            MEM( I ) = CHECKVAL
   10    CONTINUE
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            MEM( I ) = CHECKVAL
   20    CONTINUE
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K+LDA-M-1
               MEM( I ) = CHECKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
*     If the matrix is upper or lower trapezoidal, calculate the
*     additional triangular area which needs to be padded,  Each
*     element referred to is in the Ith row and the Jth column.
*
      IF( UPLO .EQ. 'U' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 41 I = 1, M
                  DO 42 J = 1, I
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   42             CONTINUE
   41          CONTINUE
            ELSE
               DO 43 I = 2, M
                  DO 44 J = 1, I-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   44             CONTINUE
   43          CONTINUE
            END IF
         ELSE
            IF( DIAG .EQ. 'U' ) THEN
               DO 45 I = M-N+1, M
                  DO 46 J = 1, I-(M-N)
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   46             CONTINUE
   45          CONTINUE
            ELSE
               DO 47 I = M-N+2, M
                  DO 48 J = 1, I-(M-N)-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   48             CONTINUE
   47          CONTINUE
            END IF
         END IF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 49 I = 1, M
                  DO 50 J = N-M+I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   50             CONTINUE
   49          CONTINUE
            ELSE
               DO 51 I = 1, M-1
                  DO 52 J = N-M+I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   52             CONTINUE
   51          CONTINUE
            END IF
         ELSE
            IF( UPLO .EQ. 'U' ) THEN
               DO 53 I = 1, N
                  DO 54 J = I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   54             CONTINUE
   53          CONTINUE
            ELSE
               DO 55 I = 1, N-1
                  DO 56 J = I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   56             CONTINUE
   55          CONTINUE
            END IF
         END IF
      END IF
*
*     End of DPADMAT.
*
      RETURN
      END
*
      SUBROUTINE DCHKPAD( UPLO, DIAG, M, N, MEM, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, IPRE, IPOST, CHECKVAL,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST
      INTEGER TESTNUM, MAXERR, NERR
      DOUBLE PRECISION CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      DOUBLE PRECISION MEM(*), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  DCHKPAD: Check padding put in by PADMAT.
*  Checks that padding around target matrix has not been overwritten
*  by the previous point-to-point or broadcast send.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM       (input) double precision array, dimension(IPRE+IPOST+LDA*N).
*            Memory location IPRE elements in front of the matrix A.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*  IPRE     (input) INTEGER
*           The size of the guard zone before the start of A.
*
*  IPOST    (input) INTEGER
*           The size of guard zone after A.
*
*  CHECKVAL (input) double precision
*           The value to pad matrix with.
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRDBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL ISTRAP
      INTEGER I, J, K, IRST, IRND, ICST, ICND, SRC, DEST
      INTEGER NPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE - I + 1
                  ERRIBUF(6, NERR) = ERR_PRE
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   10    CONTINUE
      END IF
*
*     Check buffer behind A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I - J + 1
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = ERR_POST
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         DO 40 J = 1, N
            DO 30 I = M+1, LDA
               K = IPRE + (J-1)*LDA + I
               IF( MEM(K) .NE. CHECKVAL) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_GAP
                     ERRDBUF(1, NERR) = MEM(K)
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Determine limits of trapezoidal matrix
*
      ISTRAP = .FALSE.
      IF( UPLO .EQ. 'U' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 2
            IRND = M
            ICST = 1
            ICND = M - 1
         ELSEIF( M .GT. N ) THEN
            IRST = ( M-N ) + 2
            IRND = M
            ICST = 1
            ICND = N - 1
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            IRST = IRST - 1
            ICND = ICND + 1
         ENDIF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 1
            IRND = 1
            ICST = ( N-M ) + 2
            ICND = N
         ELSEIF( M .GT. N ) THEN
            IRST = 1
            IRND = 1
            ICST = 2
            ICND = N
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            ICST = ICST - 1
         ENDIF
      ENDIF
*
*     Check elements and report any errors
*
      IF( ISTRAP ) THEN
         DO 100 J = ICST, ICND
            DO 105 I = IRST, IRND
               IF( MEM( IPRE + (J-1)*LDA + I ) .NE. CHECKVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_TRI
                     ERRDBUF(1, NERR) = MEM( IPRE + (J-1)*LDA + I )
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
  105       CONTINUE
*
*           Update the limits to allow filling in padding
*
            IF( UPLO .EQ. 'U' ) THEN
               IRST = IRST + 1
            ELSE
               IRND = IRND + 1
            ENDIF
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of DCHKPAD.
*
      END
*
      SUBROUTINE DCHKMAT( UPLO, DIAG, M, N, A, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, TESTNUM, MAXERR, NERR,
     $                    ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      DOUBLE PRECISION A(LDA,N), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  dCHKMAT:  Check matrix to see whether there were any transmission
*            errors.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (input) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRDBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC, DEST
      LOGICAL USEIT
      DOUBLE PRECISION COMPVAL
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      DOUBLE PRECISION DBTRAN
      EXTERNAL DBTRAN, IBTNPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Initialize ISEED with the same values as used in DGENMAT.
*
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
*     Generate the elements randomly with the same method used in GENMAT.
*     Note that for trapezoidal matrices, we generate all elements in the
*     enclosing rectangle and then ignore the complementary triangle.
*
      DO 100 J = 1, N
         DO 105 I = 1, M
            COMPVAL = DBTRAN( ISEED )
*
*           Now determine whether we actually need this value.  The
*           strategy is to chop out the proper triangle based on what
*           particular kind of trapezoidal matrix we're dealing with.
*
            USEIT = .TRUE.
            IF( UPLO .EQ. 'U' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF( UPLO .EQ. 'L' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J. GE. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J .GE. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            END IF
*
*           Compare the generated value to the one that's in the
*           received matrix.  If they don't match, tack another
*           error record onto what's already there.
*
            IF( USEIT ) THEN
               IF( A(I,J) .NE. COMPVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = 5
                     ERRDBUF(1, NERR) = A(I, J)
                     ERRDBUF(2, NERR) = COMPVAL
                  END IF
               END IF
            END IF
  105    CONTINUE
  100 CONTINUE
      RETURN
*
*     End of DCHKMAT.
*
      END
*
      SUBROUTINE DPRINTERRS( OUTNUM, MAXERR, NERR,
     $                       ERRIBUF, ERRDBUF, COUNTING, TFAILED )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL COUNTING
      INTEGER OUTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), TFAILED(*)
      DOUBLE PRECISION ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  DPRINTERRS: Print errors that have been recorded
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRDBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (input/ourput) INTEGER array, dimension NTESTS
*          Workspace used to keep track of which tests failed.
*          This array not accessed unless COUNTING is true.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      EXTERNAL IBTMYPROC, IBTNPROCS
*     ..
*     .. Local Scalars ..
      CHARACTER*1 MAT
      LOGICAL MATISINT
      INTEGER OLDTEST, NPROCS, PROW, PCOL, I, ERRTYPE
*     ..
*     .. Executable Statements ..
*
      IF( (IBTMYPROC().NE.0) .OR. (NERR.LE.0) ) RETURN
      OLDTEST = -1
      NPROCS = IBTNPROCS()
      PROW = ERRIBUF(3,1) / NPROCS
      PCOL = MOD( ERRIBUF(3,1), NPROCS )
      IF( NERR .GT. MAXERR ) WRITE(OUTNUM,13000)
*
      DO 20 I = 1, MIN( NERR, MAXERR )
         IF( ERRIBUF(1,I) .NE. OLDTEST ) THEN
            IF( OLDTEST .NE. -1 )
     $         WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM,1000) PROW, PCOL, ERRIBUF(1,I)
            IF( COUNTING ) TFAILED( ERRIBUF(1,I) ) = 1
            OLDTEST = ERRIBUF(1, I)
         END IF
*
*        Print out error message depending on type of error
*
         ERRTYPE = ERRIBUF(6, I)
         IF( ERRTYPE .LT. -10 ) THEN
            ERRTYPE = -ERRTYPE - 10
            MAT = 'C'
            MATISINT = .TRUE.
         ELSE IF( ERRTYPE .LT. 0 ) THEN
            ERRTYPE = -ERRTYPE
            MAT = 'R'
            MATISINT = .TRUE.
         ELSE
            MATISINT = .FALSE.
         END IF
*
*        RA/CA arrays from MAX/MIN have different printing protocol
*
         IF( MATISINT ) THEN
            IF( ERRIBUF(2, I) .EQ. -1 ) THEN
               WRITE(OUTNUM,11000) ERRIBUF(4,I), ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,7000) ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,8000) ERRIBUF(4,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,9000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,10000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $                             INT( ERRDBUF(2,I) ),
     $                             INT( ERRDBUF(1,I) )
            END IF
*
*        Have memory overwrites in matrix A
*
         ELSE
            IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,2000) ERRIBUF(5,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,3000) ERRIBUF(4,I), ERRDBUF(2,I),
     $                            ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,4000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE IF( ERRTYPE .EQ. ERR_TRI ) THEN
               WRITE(OUTNUM,5000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            ELSE
               WRITE(OUTNUM,6000) ERRIBUF(4,I), ERRIBUF(5,I),
     $                            ERRDBUF(2,I), ERRDBUF(1,I)
            END IF
         END IF
   20 CONTINUE
      WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
*
 1000 FORMAT('PROCESS {',I4,',',I4,'} REPORTS ERRORS IN TEST#',I6,':')
 2000 FORMAT('   Buffer overwrite ',I4,
     $       ' elements before the start of A:',/,
     $       '   Expected=',G22.15,
     $       '; Received=',G22.15)
 3000 FORMAT('   Buffer overwrite ',I4,' elements after the end of A:',
     $       /,'   Expected=',G22.15,
     $       '; Received=',G22.15)
 4000 FORMAT('   LDA-M gap overwrite at postion (',I4,',',I4,'):',/,
     $       '   Expected=',G22.15,
     $       '; Received=',G22.15)
 5000 FORMAT('   Complementory triangle overwrite at A(',I4,',',I4,
     $       '):',/,'   Expected=',G22.15,
     $       '; Received=',G22.15)
 6000 FORMAT('   Invalid element at A(',I4,',',I4,'):',/,
     $       '   Expected=',G22.15,
     $       '; Received=',G22.15)
 7000 FORMAT('   Buffer overwrite ',I4,' elements before the start of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
 8000 FORMAT('   Buffer overwrite ',I4,' elements after the end of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
*
 9000 FORMAT('   LD',A1,'A-M gap overwrite at postion (',I4,',',I4,'):'
     $       ,/,'   Expected=',I12,'; Received=',I12)
*
10000 FORMAT('   Invalid element at ',A1,'A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,'; Received=',I12)
11000 FORMAT('   Overwrite at position (',I4,',',I4,') of non-existent '
     $       ,A1,'A array.',/,'   Expected=',I12,'; Received=',I12)
12000 FORMAT('PROCESS {',I4,',',I4,'} DONE ERROR REPORT FOR TEST#',
     $       I6,'.')
13000 FORMAT('WARNING: There were more errors than could be recorded.',
     $       /,'Increase MEMELTS to get complete listing.')
      RETURN
*
*     End DPRINTERRS
*
      END
*
*
      SUBROUTINE CBTCHECKIN( NFTESTS, OUTNUM, MAXERR, NERR, IERR,
     $                       CVAL, TFAILED )
      INTEGER NFTESTS, OUTNUM, MAXERR, NERR
      INTEGER IERR(*), TFAILED(*)
      COMPLEX CVAL(*)
*
*  Purpose
*  =======
*  CBTCHECKIN: Process 0 receives error report from all processes.
*
*  Arguments
*  =========
*  NFTESTS  (input/output) INTEGER
*           if NFTESTS is <= 0 upon entry, NFTESTS is not written to.
*           Otherwise, on entry it specifies the total number of tests
*           run, and on exit it is the number of tests which failed.
*
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRCBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (workspace) INTEGER array, dimension NFTESTS
*          Workspace used to keep track of which tests failed.
*          If input of NFTESTS < 1, this array not accessed.
*
*  ===================================================================
*
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTMSGID
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID
*     ..
*     .. Local Scalars ..
      LOGICAL COUNTING
      INTEGER K, NERR2, IAM, NPROCS, NTESTS
*
*     Proc 0 collects error info from everyone
*
      IAM = IBTMYPROC()
      NPROCS = IBTNPROCS()
*
      IF( IAM .EQ. 0 ) THEN
*
*        If we are finding out how many failed tests there are, initialize
*        the total number of tests (NTESTS), and zero the test failed array
*
         COUNTING = NFTESTS .GT. 0
         IF( COUNTING ) THEN
            NTESTS = NFTESTS
            DO 10 K = 1, NTESTS
               TFAILED(K) = 0
   10       CONTINUE
         END IF
*
         CALL CPRINTERRS(OUTNUM, MAXERR, NERR, IERR, CVAL, COUNTING,
     $                   TFAILED)
*
         DO 20 K = 1, NPROCS-1
            CALL BTSEND(3, 0, K, K, IBTMSGID()+50)
            CALL BTRECV(3, 1, NERR2, K, IBTMSGID()+50)
            IF( NERR2 .GT. 0 ) THEN
               NERR = NERR + NERR2
               CALL BTRECV(3, NERR2*6, IERR, K, IBTMSGID()+51)
               CALL BTRECV(5, NERR2*2, CVAL, K, IBTMSGID()+51)
               CALL CPRINTERRS(OUTNUM, MAXERR, NERR2, IERR, CVAL,
     $                         COUNTING, TFAILED)
            END IF
   20    CONTINUE
*
*        Count up number of tests that failed
*
         IF( COUNTING ) THEN
            NFTESTS = 0
            DO 30 K = 1, NTESTS
               NFTESTS = NFTESTS + TFAILED(K)
   30       CONTINUE
         END IF
*
*     Send my error info to proc 0
*
      ELSE
         CALL BTRECV(3, 0, K, 0, IBTMSGID()+50)
         CALL BTSEND(3, 1, NERR, 0, IBTMSGID()+50)
         IF( NERR .GT. 0 ) THEN
            CALL BTSEND(3, NERR*6, IERR, 0, IBTMSGID()+51)
            CALL BTSEND(5, NERR*2, CVAL, 0, IBTMSGID()+51)
         END IF
      ENDIF
*
      RETURN
*
*     End of CBTCHECKIN
*
      END
*
      SUBROUTINE CINITMAT(UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL, TESTNUM, MYROW, MYCOL)
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST, TESTNUM, MYROW, MYCOL
      COMPLEX CHECKVAL
      COMPLEX MEM(*)
*
*     .. External Subroutines ..
      EXTERNAL CGENMAT, CPADMAT
*     ..
*     .. Executable Statements ..
*
      CALL CGENMAT( M, N, MEM(IPRE+1), LDA, TESTNUM, MYROW, MYCOL )
      CALL CPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST, CHECKVAL )
*
      RETURN
      END
*
      SUBROUTINE CGENMAT( M, N, A, LDA, TESTNUM, MYROW, MYCOL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA, TESTNUM, MYROW, MYCOL
*     ..
*     .. Array Arguments ..
      COMPLEX A(LDA,N)
*     ..
*
*  Purpose
*  =======
*  CGENMAT: Generates an M-by-N matrix filled with random elements.
*
*  Arguments
*  =========
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (output) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  TESTNUM  (input) INTEGER
*           Unique number for this test case, used as a basis for
*           the random seeds.
*
*  ====================================================================
*
*     .. External Functions ..
      INTEGER IBTNPROCS
      COMPLEX CBTRAN
      EXTERNAL CBTRAN, IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. Executable Statements ..
*
*     ISEED's four values must be positive integers less than 4096,
*     fourth one has to be odd. (see _LARND).  Use some goofy
*     functions to come up with seed values which together should
*     be unique.
*
      NPROCS = IBTNPROCS()
      SRC = MYROW * NPROCS + MYCOL
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
      DO 10 J = 1, N
         DO 10 I = 1, M
            A(I, J) = CBTRAN( ISEED )
   10 CONTINUE
*
      RETURN
*
*     End of CGENMAT.
*
      END
*
      COMPLEX FUNCTION CBTRAN(ISEED)
      INTEGER ISEED(*)
*
*     .. External Functions ..
      DOUBLE COMPLEX ZLARND
      EXTERNAL ZLARND
      CBTRAN = CMPLX( ZLARND(2, ISEED) )
*
      RETURN
*
*     End of Cbtran
*
      END
*
      SUBROUTINE CPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST
      COMPLEX CHECKVAL
*     ..
*     .. Array Arguments ..
      COMPLEX MEM( * )
*     ..
*
*  Purpose
*  =======
*
*  CPADMAT: Pad Matrix.
*  This routines surrounds a matrix with a guardzone initialized to the
*  value CHECKVAL.  There are three distinct guardzones:
*  - A contiguous zone of size IPRE immediately before the start
*    of the matrix.
*  - A contiguous zone of size IPOST immedately after the end of the
*    matrix.
*  - Interstitial zones within each column of the matrix, in the
*    elements A( M+1:LDA, J ).
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM      (output) complex array, dimension (IPRE+IPOST+LDA*N)
*           The address IPRE elements ahead of the matrix A you want to
*           pad, which is then of dimension (LDA,N).
*
*  IPRE     (input) INTEGER
*           The size of the guard zone ahead of the matrix A.
*
*  IPOST    (input) INTEGER
*           The size of the guard zone behind the matrix A.
*
*  CHECKVAL (input) complex
*           The value to insert into the guard zones.
*
*  ====================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            MEM( I ) = CHECKVAL
   10    CONTINUE
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            MEM( I ) = CHECKVAL
   20    CONTINUE
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K+LDA-M-1
               MEM( I ) = CHECKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
*     If the matrix is upper or lower trapezoidal, calculate the
*     additional triangular area which needs to be padded,  Each
*     element referred to is in the Ith row and the Jth column.
*
      IF( UPLO .EQ. 'U' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 41 I = 1, M
                  DO 42 J = 1, I
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   42             CONTINUE
   41          CONTINUE
            ELSE
               DO 43 I = 2, M
                  DO 44 J = 1, I-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   44             CONTINUE
   43          CONTINUE
            END IF
         ELSE
            IF( DIAG .EQ. 'U' ) THEN
               DO 45 I = M-N+1, M
                  DO 46 J = 1, I-(M-N)
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   46             CONTINUE
   45          CONTINUE
            ELSE
               DO 47 I = M-N+2, M
                  DO 48 J = 1, I-(M-N)-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   48             CONTINUE
   47          CONTINUE
            END IF
         END IF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 49 I = 1, M
                  DO 50 J = N-M+I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   50             CONTINUE
   49          CONTINUE
            ELSE
               DO 51 I = 1, M-1
                  DO 52 J = N-M+I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   52             CONTINUE
   51          CONTINUE
            END IF
         ELSE
            IF( UPLO .EQ. 'U' ) THEN
               DO 53 I = 1, N
                  DO 54 J = I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   54             CONTINUE
   53          CONTINUE
            ELSE
               DO 55 I = 1, N-1
                  DO 56 J = I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   56             CONTINUE
   55          CONTINUE
            END IF
         END IF
      END IF
*
*     End of CPADMAT.
*
      RETURN
      END
*
      SUBROUTINE CCHKPAD( UPLO, DIAG, M, N, MEM, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, IPRE, IPOST, CHECKVAL,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST
      INTEGER TESTNUM, MAXERR, NERR
      COMPLEX CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      COMPLEX MEM(*), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  CCHKPAD: Check padding put in by PADMAT.
*  Checks that padding around target matrix has not been overwritten
*  by the previous point-to-point or broadcast send.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM       (input) complex array, dimension(IPRE+IPOST+LDA*N).
*            Memory location IPRE elements in front of the matrix A.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*  IPRE     (input) INTEGER
*           The size of the guard zone before the start of A.
*
*  IPOST    (input) INTEGER
*           The size of guard zone after A.
*
*  CHECKVAL (input) complex
*           The value to pad matrix with.
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRCBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL ISTRAP
      INTEGER I, J, K, IRST, IRND, ICST, ICND, SRC, DEST
      INTEGER NPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE - I + 1
                  ERRIBUF(6, NERR) = ERR_PRE
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   10    CONTINUE
      END IF
*
*     Check buffer behind A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I - J + 1
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = ERR_POST
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         DO 40 J = 1, N
            DO 30 I = M+1, LDA
               K = IPRE + (J-1)*LDA + I
               IF( MEM(K) .NE. CHECKVAL) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_GAP
                     ERRDBUF(1, NERR) = MEM(K)
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Determine limits of trapezoidal matrix
*
      ISTRAP = .FALSE.
      IF( UPLO .EQ. 'U' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 2
            IRND = M
            ICST = 1
            ICND = M - 1
         ELSEIF( M .GT. N ) THEN
            IRST = ( M-N ) + 2
            IRND = M
            ICST = 1
            ICND = N - 1
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            IRST = IRST - 1
            ICND = ICND + 1
         ENDIF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 1
            IRND = 1
            ICST = ( N-M ) + 2
            ICND = N
         ELSEIF( M .GT. N ) THEN
            IRST = 1
            IRND = 1
            ICST = 2
            ICND = N
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            ICST = ICST - 1
         ENDIF
      ENDIF
*
*     Check elements and report any errors
*
      IF( ISTRAP ) THEN
         DO 100 J = ICST, ICND
            DO 105 I = IRST, IRND
               IF( MEM( IPRE + (J-1)*LDA + I ) .NE. CHECKVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_TRI
                     ERRDBUF(1, NERR) = MEM( IPRE + (J-1)*LDA + I )
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
  105       CONTINUE
*
*           Update the limits to allow filling in padding
*
            IF( UPLO .EQ. 'U' ) THEN
               IRST = IRST + 1
            ELSE
               IRND = IRND + 1
            ENDIF
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of CCHKPAD.
*
      END
*
      SUBROUTINE CCHKMAT( UPLO, DIAG, M, N, A, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, TESTNUM, MAXERR, NERR,
     $                    ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      COMPLEX A(LDA,N), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  cCHKMAT:  Check matrix to see whether there were any transmission
*            errors.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (input) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRCBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC, DEST
      LOGICAL USEIT
      COMPLEX COMPVAL
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      COMPLEX CBTRAN
      EXTERNAL CBTRAN, IBTNPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Initialize ISEED with the same values as used in CGENMAT.
*
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
*     Generate the elements randomly with the same method used in GENMAT.
*     Note that for trapezoidal matrices, we generate all elements in the
*     enclosing rectangle and then ignore the complementary triangle.
*
      DO 100 J = 1, N
         DO 105 I = 1, M
            COMPVAL = CBTRAN( ISEED )
*
*           Now determine whether we actually need this value.  The
*           strategy is to chop out the proper triangle based on what
*           particular kind of trapezoidal matrix we're dealing with.
*
            USEIT = .TRUE.
            IF( UPLO .EQ. 'U' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF( UPLO .EQ. 'L' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J. GE. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J .GE. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            END IF
*
*           Compare the generated value to the one that's in the
*           received matrix.  If they don't match, tack another
*           error record onto what's already there.
*
            IF( USEIT ) THEN
               IF( A(I,J) .NE. COMPVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = 5
                     ERRDBUF(1, NERR) = A(I, J)
                     ERRDBUF(2, NERR) = COMPVAL
                  END IF
               END IF
            END IF
  105    CONTINUE
  100 CONTINUE
      RETURN
*
*     End of CCHKMAT.
*
      END
*
      SUBROUTINE CPRINTERRS( OUTNUM, MAXERR, NERR,
     $                       ERRIBUF, ERRDBUF, COUNTING, TFAILED )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL COUNTING
      INTEGER OUTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), TFAILED(*)
      COMPLEX ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  CPRINTERRS: Print errors that have been recorded
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRCBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (input/ourput) INTEGER array, dimension NTESTS
*          Workspace used to keep track of which tests failed.
*          This array not accessed unless COUNTING is true.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      EXTERNAL IBTMYPROC, IBTNPROCS
*     ..
*     .. Local Scalars ..
      CHARACTER*1 MAT
      LOGICAL MATISINT
      INTEGER OLDTEST, NPROCS, PROW, PCOL, I, ERRTYPE
*     ..
*     .. Executable Statements ..
*
      IF( (IBTMYPROC().NE.0) .OR. (NERR.LE.0) ) RETURN
      OLDTEST = -1
      NPROCS = IBTNPROCS()
      PROW = ERRIBUF(3,1) / NPROCS
      PCOL = MOD( ERRIBUF(3,1), NPROCS )
      IF( NERR .GT. MAXERR ) WRITE(OUTNUM,13000)
*
      DO 20 I = 1, MIN( NERR, MAXERR )
         IF( ERRIBUF(1,I) .NE. OLDTEST ) THEN
            IF( OLDTEST .NE. -1 )
     $         WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM,1000) PROW, PCOL, ERRIBUF(1,I)
            IF( COUNTING ) TFAILED( ERRIBUF(1,I) ) = 1
            OLDTEST = ERRIBUF(1, I)
         END IF
*
*        Print out error message depending on type of error
*
         ERRTYPE = ERRIBUF(6, I)
         IF( ERRTYPE .LT. -10 ) THEN
            ERRTYPE = -ERRTYPE - 10
            MAT = 'C'
            MATISINT = .TRUE.
         ELSE IF( ERRTYPE .LT. 0 ) THEN
            ERRTYPE = -ERRTYPE
            MAT = 'R'
            MATISINT = .TRUE.
         ELSE
            MATISINT = .FALSE.
         END IF
*
*        RA/CA arrays from MAX/MIN have different printing protocol
*
         IF( MATISINT ) THEN
            IF( ERRIBUF(2, I) .EQ. -1 ) THEN
               WRITE(OUTNUM,11000) ERRIBUF(4,I), ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,7000) ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,8000) ERRIBUF(4,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,9000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,10000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $                             INT( ERRDBUF(2,I) ),
     $                             INT( ERRDBUF(1,I) )
            END IF
*
*        Have memory overwrites in matrix A
*
         ELSE
            IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,2000) ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), AIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), AIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,3000) ERRIBUF(4,I),
     $         REAL( ERRDBUF(2,I) ), AIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), AIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,4000)
     $         ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), AIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), AIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_TRI ) THEN
               WRITE(OUTNUM,5000) ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), AIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), AIMAG( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,6000) ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), AIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), AIMAG( ERRDBUF(1,I) )
            END IF
         END IF
   20 CONTINUE
      WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
*
 1000 FORMAT('PROCESS {',I4,',',I4,'} REPORTS ERRORS IN TEST#',I6,':')
 2000 FORMAT('   Buffer overwrite ',I4,
     $       ' elements before the start of A:',/,
     $       '   Expected=','[',G15.8,',',G15.8,']',
     $       '; Received=','[',G15.8,',',G15.8,']')
 3000 FORMAT('   Buffer overwrite ',I4,' elements after the end of A:',
     $       /,'   Expected=','[',G15.8,',',G15.8,']',
     $       '; Received=','[',G15.8,',',G15.8,']')
 4000 FORMAT('   LDA-M gap overwrite at postion (',I4,',',I4,'):',/,
     $       '   Expected=','[',G15.8,',',G15.8,']',
     $       '; Received=','[',G15.8,',',G15.8,']')
 5000 FORMAT('   Complementory triangle overwrite at A(',I4,',',I4,
     $       '):',/,'   Expected=','[',G15.8,',',G15.8,']',
     $       '; Received=','[',G15.8,',',G15.8,']')
 6000 FORMAT('   Invalid element at A(',I4,',',I4,'):',/,
     $       '   Expected=','[',G15.8,',',G15.8,']',
     $       '; Received=','[',G15.8,',',G15.8,']')
 7000 FORMAT('   Buffer overwrite ',I4,' elements before the start of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
 8000 FORMAT('   Buffer overwrite ',I4,' elements after the end of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
*
 9000 FORMAT('   LD',A1,'A-M gap overwrite at postion (',I4,',',I4,'):'
     $       ,/,'   Expected=',I12,'; Received=',I12)
*
10000 FORMAT('   Invalid element at ',A1,'A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,'; Received=',I12)
11000 FORMAT('   Overwrite at position (',I4,',',I4,') of non-existent '
     $       ,A1,'A array.',/,'   Expected=',I12,'; Received=',I12)
12000 FORMAT('PROCESS {',I4,',',I4,'} DONE ERROR REPORT FOR TEST#',
     $       I6,'.')
13000 FORMAT('WARNING: There were more errors than could be recorded.',
     $       /,'Increase MEMELTS to get complete listing.')
      RETURN
*
*     End CPRINTERRS
*
      END
*
*
      SUBROUTINE ZBTCHECKIN( NFTESTS, OUTNUM, MAXERR, NERR, IERR,
     $                       ZVAL, TFAILED )
      INTEGER NFTESTS, OUTNUM, MAXERR, NERR
      INTEGER IERR(*), TFAILED(*)
      DOUBLE COMPLEX ZVAL(*)
*
*  Purpose
*  =======
*  ZBTCHECKIN: Process 0 receives error report from all processes.
*
*  Arguments
*  =========
*  NFTESTS  (input/output) INTEGER
*           if NFTESTS is <= 0 upon entry, NFTESTS is not written to.
*           Otherwise, on entry it specifies the total number of tests
*           run, and on exit it is the number of tests which failed.
*
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRZBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (workspace) INTEGER array, dimension NFTESTS
*          Workspace used to keep track of which tests failed.
*          If input of NFTESTS < 1, this array not accessed.
*
*  ===================================================================
*
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTMSGID
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID
*     ..
*     .. Local Scalars ..
      LOGICAL COUNTING
      INTEGER K, NERR2, IAM, NPROCS, NTESTS
*
*     Proc 0 collects error info from everyone
*
      IAM = IBTMYPROC()
      NPROCS = IBTNPROCS()
*
      IF( IAM .EQ. 0 ) THEN
*
*        If we are finding out how many failed tests there are, initialize
*        the total number of tests (NTESTS), and zero the test failed array
*
         COUNTING = NFTESTS .GT. 0
         IF( COUNTING ) THEN
            NTESTS = NFTESTS
            DO 10 K = 1, NTESTS
               TFAILED(K) = 0
   10       CONTINUE
         END IF
*
         CALL ZPRINTERRS(OUTNUM, MAXERR, NERR, IERR, ZVAL, COUNTING,
     $                   TFAILED)
*
         DO 20 K = 1, NPROCS-1
            CALL BTSEND(3, 0, K, K, IBTMSGID()+50)
            CALL BTRECV(3, 1, NERR2, K, IBTMSGID()+50)
            IF( NERR2 .GT. 0 ) THEN
               NERR = NERR + NERR2
               CALL BTRECV(3, NERR2*6, IERR, K, IBTMSGID()+51)
               CALL BTRECV(7, NERR2*2, ZVAL, K, IBTMSGID()+51)
               CALL ZPRINTERRS(OUTNUM, MAXERR, NERR2, IERR, ZVAL,
     $                         COUNTING, TFAILED)
            END IF
   20    CONTINUE
*
*        Count up number of tests that failed
*
         IF( COUNTING ) THEN
            NFTESTS = 0
            DO 30 K = 1, NTESTS
               NFTESTS = NFTESTS + TFAILED(K)
   30       CONTINUE
         END IF
*
*     Send my error info to proc 0
*
      ELSE
         CALL BTRECV(3, 0, K, 0, IBTMSGID()+50)
         CALL BTSEND(3, 1, NERR, 0, IBTMSGID()+50)
         IF( NERR .GT. 0 ) THEN
            CALL BTSEND(3, NERR*6, IERR, 0, IBTMSGID()+51)
            CALL BTSEND(7, NERR*2, ZVAL, 0, IBTMSGID()+51)
         END IF
      ENDIF
*
      RETURN
*
*     End of ZBTCHECKIN
*
      END
*
      SUBROUTINE ZINITMAT(UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL, TESTNUM, MYROW, MYCOL)
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST, TESTNUM, MYROW, MYCOL
      DOUBLE COMPLEX CHECKVAL
      DOUBLE COMPLEX MEM(*)
*
*     .. External Subroutines ..
      EXTERNAL ZGENMAT, ZPADMAT
*     ..
*     .. Executable Statements ..
*
      CALL ZGENMAT( M, N, MEM(IPRE+1), LDA, TESTNUM, MYROW, MYCOL )
      CALL ZPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST, CHECKVAL )
*
      RETURN
      END
*
      SUBROUTINE ZGENMAT( M, N, A, LDA, TESTNUM, MYROW, MYCOL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER M, N, LDA, TESTNUM, MYROW, MYCOL
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX A(LDA,N)
*     ..
*
*  Purpose
*  =======
*  ZGENMAT: Generates an M-by-N matrix filled with random elements.
*
*  Arguments
*  =========
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (output) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  TESTNUM  (input) INTEGER
*           Unique number for this test case, used as a basis for
*           the random seeds.
*
*  ====================================================================
*
*     .. External Functions ..
      INTEGER IBTNPROCS
      DOUBLE COMPLEX ZBTRAN
      EXTERNAL ZBTRAN, IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. Executable Statements ..
*
*     ISEED's four values must be positive integers less than 4096,
*     fourth one has to be odd. (see _LARND).  Use some goofy
*     functions to come up with seed values which together should
*     be unique.
*
      NPROCS = IBTNPROCS()
      SRC = MYROW * NPROCS + MYCOL
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
      DO 10 J = 1, N
         DO 10 I = 1, M
            A(I, J) = ZBTRAN( ISEED )
   10 CONTINUE
*
      RETURN
*
*     End of ZGENMAT.
*
      END
*
      DOUBLE COMPLEX FUNCTION ZBTRAN(ISEED)
      INTEGER ISEED(*)
*
*     .. External Functions ..
      DOUBLE COMPLEX ZLARND
      EXTERNAL ZLARND
      ZBTRAN = ZLARND(2, ISEED)
*
      RETURN
*
*     End of Zbtran
*
      END
*
      SUBROUTINE ZPADMAT( UPLO, DIAG, M, N, MEM, LDA, IPRE, IPOST,
     $                    CHECKVAL )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, IPRE, IPOST
      DOUBLE COMPLEX CHECKVAL
*     ..
*     .. Array Arguments ..
      DOUBLE COMPLEX MEM( * )
*     ..
*
*  Purpose
*  =======
*
*  ZPADMAT: Pad Matrix.
*  This routines surrounds a matrix with a guardzone initialized to the
*  value CHECKVAL.  There are three distinct guardzones:
*  - A contiguous zone of size IPRE immediately before the start
*    of the matrix.
*  - A contiguous zone of size IPOST immedately after the end of the
*    matrix.
*  - Interstitial zones within each column of the matrix, in the
*    elements A( M+1:LDA, J ).
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM      (output) double complex array, dimension (IPRE+IPOST+LDA*N)
*           The address IPRE elements ahead of the matrix A you want to
*           pad, which is then of dimension (LDA,N).
*
*  IPRE     (input) INTEGER
*           The size of the guard zone ahead of the matrix A.
*
*  IPOST    (input) INTEGER
*           The size of the guard zone behind the matrix A.
*
*  CHECKVAL (input) double complex
*           The value to insert into the guard zones.
*
*  ====================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, K
*     ..
*     .. Executable Statements ..
*
*     Put check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            MEM( I ) = CHECKVAL
   10    CONTINUE
      END IF
*
*     Put check buffer in back of A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            MEM( I ) = CHECKVAL
   20    CONTINUE
      END IF
*
*     Put check buffer in all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         K = IPRE + M + 1
         DO 40 J = 1, N
            DO 30 I = K, K+LDA-M-1
               MEM( I ) = CHECKVAL
   30       CONTINUE
            K = K + LDA
   40    CONTINUE
      END IF
*
*     If the matrix is upper or lower trapezoidal, calculate the
*     additional triangular area which needs to be padded,  Each
*     element referred to is in the Ith row and the Jth column.
*
      IF( UPLO .EQ. 'U' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 41 I = 1, M
                  DO 42 J = 1, I
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   42             CONTINUE
   41          CONTINUE
            ELSE
               DO 43 I = 2, M
                  DO 44 J = 1, I-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   44             CONTINUE
   43          CONTINUE
            END IF
         ELSE
            IF( DIAG .EQ. 'U' ) THEN
               DO 45 I = M-N+1, M
                  DO 46 J = 1, I-(M-N)
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   46             CONTINUE
   45          CONTINUE
            ELSE
               DO 47 I = M-N+2, M
                  DO 48 J = 1, I-(M-N)-1
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   48             CONTINUE
   47          CONTINUE
            END IF
         END IF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         IF( M .LE. N ) THEN
            IF( DIAG .EQ. 'U' ) THEN
               DO 49 I = 1, M
                  DO 50 J = N-M+I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   50             CONTINUE
   49          CONTINUE
            ELSE
               DO 51 I = 1, M-1
                  DO 52 J = N-M+I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   52             CONTINUE
   51          CONTINUE
            END IF
         ELSE
            IF( UPLO .EQ. 'U' ) THEN
               DO 53 I = 1, N
                  DO 54 J = I, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   54             CONTINUE
   53          CONTINUE
            ELSE
               DO 55 I = 1, N-1
                  DO 56 J = I+1, N
                     K = IPRE + I + (J-1)*LDA
                     MEM( K ) = CHECKVAL
   56             CONTINUE
   55          CONTINUE
            END IF
         END IF
      END IF
*
*     End of ZPADMAT.
*
      RETURN
      END
*
      SUBROUTINE ZCHKPAD( UPLO, DIAG, M, N, MEM, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, IPRE, IPOST, CHECKVAL,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, IPRE, IPOST
      INTEGER TESTNUM, MAXERR, NERR
      DOUBLE COMPLEX CHECKVAL
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      DOUBLE COMPLEX MEM(*), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  ZCHKPAD: Check padding put in by PADMAT.
*  Checks that padding around target matrix has not been overwritten
*  by the previous point-to-point or broadcast send.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*  MEM       (input) double complex array, dimension(IPRE+IPOST+LDA*N).
*            Memory location IPRE elements in front of the matrix A.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*  IPRE     (input) INTEGER
*           The size of the guard zone before the start of A.
*
*  IPOST    (input) INTEGER
*           The size of guard zone after A.
*
*  CHECKVAL (input) double complex
*           The value to pad matrix with.
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRZBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      LOGICAL ISTRAP
      INTEGER I, J, K, IRST, IRND, ICST, ICND, SRC, DEST
      INTEGER NPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Check buffer in front of A
*
      IF( IPRE .GT. 0 ) THEN
         DO 10 I = 1, IPRE
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE - I + 1
                  ERRIBUF(6, NERR) = ERR_PRE
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   10    CONTINUE
      END IF
*
*     Check buffer behind A
*
      IF( IPOST .GT. 0 ) THEN
         J = IPRE + LDA*N + 1
         DO 20 I = J, J+IPOST-1
            IF( MEM(I) .NE. CHECKVAL ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = SRC
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I - J + 1
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = ERR_POST
                  ERRDBUF(1, NERR) = MEM(I)
                  ERRDBUF(2, NERR) = CHECKVAL
               END IF
            END IF
   20    CONTINUE
      END IF
*
*     Check all (LDA-M) gaps
*
      IF( LDA .GT. M ) THEN
         DO 40 J = 1, N
            DO 30 I = M+1, LDA
               K = IPRE + (J-1)*LDA + I
               IF( MEM(K) .NE. CHECKVAL) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_GAP
                     ERRDBUF(1, NERR) = MEM(K)
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
   30       CONTINUE
   40    CONTINUE
      END IF
*
*     Determine limits of trapezoidal matrix
*
      ISTRAP = .FALSE.
      IF( UPLO .EQ. 'U' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 2
            IRND = M
            ICST = 1
            ICND = M - 1
         ELSEIF( M .GT. N ) THEN
            IRST = ( M-N ) + 2
            IRND = M
            ICST = 1
            ICND = N - 1
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            IRST = IRST - 1
            ICND = ICND + 1
         ENDIF
      ELSE IF( UPLO .EQ. 'L' ) THEN
         ISTRAP = .TRUE.
         IF( M .LE. N ) THEN
            IRST = 1
            IRND = 1
            ICST = ( N-M ) + 2
            ICND = N
         ELSEIF( M .GT. N ) THEN
            IRST = 1
            IRND = 1
            ICST = 2
            ICND = N
         ENDIF
         IF( DIAG .EQ. 'U' ) THEN
            ICST = ICST - 1
         ENDIF
      ENDIF
*
*     Check elements and report any errors
*
      IF( ISTRAP ) THEN
         DO 100 J = ICST, ICND
            DO 105 I = IRST, IRND
               IF( MEM( IPRE + (J-1)*LDA + I ) .NE. CHECKVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = ERR_TRI
                     ERRDBUF(1, NERR) = MEM( IPRE + (J-1)*LDA + I )
                     ERRDBUF(2, NERR) = CHECKVAL
                  END IF
               END IF
  105       CONTINUE
*
*           Update the limits to allow filling in padding
*
            IF( UPLO .EQ. 'U' ) THEN
               IRST = IRST + 1
            ELSE
               IRND = IRND + 1
            ENDIF
  100    CONTINUE
      END IF
*
      RETURN
*
*     End of ZCHKPAD.
*
      END
*
      SUBROUTINE ZCHKMAT( UPLO, DIAG, M, N, A, LDA, RSRC, CSRC,
     $                    MYROW, MYCOL, TESTNUM, MAXERR, NERR,
     $                    ERRIBUF, ERRDBUF )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      CHARACTER*1 UPLO, DIAG
      INTEGER M, N, LDA, RSRC, CSRC, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR)
      DOUBLE COMPLEX A(LDA,N), ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  zCHKMAT:  Check matrix to see whether there were any transmission
*            errors.
*
*  Arguments
*  =========
*  UPLO     (input) CHARACTER*1
*           Is the matrix A 'U'pper or 'L'ower trapezoidal, or 'G'eneral
*           rectangular?
*
*  DIAG     (input) CHARACTER*1
*           For trapezoidal matrices, is the main diagonal included
*           ('N') or not ('U')?
*
*   M       (input) INTEGER
*           The number of rows of the matrix A.  M >= 0.
*
*   N       (input) INTEGER
*           The number of columns of the matrix A.  N >= 0.
*
*   A       (input) @up@(doctype) array, dimension (LDA,N)
*           The m by n matrix A.  Fortran77 (column-major) storage
*           assumed.
*
*   LDA     (input) INTEGER
*           The leading dimension of the array A.  LDA >= max(1, M).
*
*  RSRC     (input) INTEGER
*           The process row of the source of the matrix.
*
*  CSRC     (input) INTEGER
*           The process column of the source of the matrix.
*
*  MYROW    (input) INTEGER
*           Row of this process in the process grid.
*
*  MYCOL    (input) INTEGER
*           Column of this process in the process grid.
*
*
*  TESTNUM  (input) INTEGER
*           The number of the test being checked.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRZBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  ===================================================================
*
*     .. Local Scalars ..
      INTEGER I, J, NPROCS, SRC, DEST
      LOGICAL USEIT
      DOUBLE COMPLEX COMPVAL
*     ..
*     .. Local Arrays ..
      INTEGER ISEED(4)
*     ..
*     .. External Functions ..
      INTEGER IBTNPROCS
      DOUBLE COMPLEX ZBTRAN
      EXTERNAL ZBTRAN, IBTNPROCS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      SRC = RSRC * NPROCS + CSRC
      DEST = MYROW * NPROCS + MYCOL
*
*     Initialize ISEED with the same values as used in ZGENMAT.
*
      ISEED(1) = MOD( 1002 + TESTNUM*5 + SRC*3, 4096 )
      ISEED(2) = MOD( 2027 + TESTNUM*7 + SRC, 4096 )
      ISEED(3) = MOD( 1234 + TESTNUM + SRC*3, 4096 )
      ISEED(4) = MOD( 4311 + TESTNUM*10 + SRC*2, 4096 )
*
*     Generate the elements randomly with the same method used in GENMAT.
*     Note that for trapezoidal matrices, we generate all elements in the
*     enclosing rectangle and then ignore the complementary triangle.
*
      DO 100 J = 1, N
         DO 105 I = 1, M
            COMPVAL = ZBTRAN( ISEED )
*
*           Now determine whether we actually need this value.  The
*           strategy is to chop out the proper triangle based on what
*           particular kind of trapezoidal matrix we're dealing with.
*
            USEIT = .TRUE.
            IF( UPLO .EQ. 'U' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( I .GE. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( I .GT. M-N+J ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            ELSE IF( UPLO .EQ. 'L' ) THEN
               IF( M .LE. N ) THEN
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J. GE. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I+(N-M) ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               ELSE
                  IF( DIAG .EQ. 'U' ) THEN
                     IF( J .GE. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  ELSE
                     IF( J .GT. I ) THEN
                        USEIT = .FALSE.
                     END IF
                  END IF
               END IF
            END IF
*
*           Compare the generated value to the one that's in the
*           received matrix.  If they don't match, tack another
*           error record onto what's already there.
*
            IF( USEIT ) THEN
               IF( A(I,J) .NE. COMPVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = SRC
                     ERRIBUF(3, NERR) = DEST
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = J
                     ERRIBUF(6, NERR) = 5
                     ERRDBUF(1, NERR) = A(I, J)
                     ERRDBUF(2, NERR) = COMPVAL
                  END IF
               END IF
            END IF
  105    CONTINUE
  100 CONTINUE
      RETURN
*
*     End of ZCHKMAT.
*
      END
*
      SUBROUTINE ZPRINTERRS( OUTNUM, MAXERR, NERR,
     $                       ERRIBUF, ERRDBUF, COUNTING, TFAILED )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      LOGICAL COUNTING
      INTEGER OUTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), TFAILED(*)
      DOUBLE COMPLEX ERRDBUF(2, MAXERR)
*     ..
*
*  Purpose
*  =======
*  ZPRINTERRS: Print errors that have been recorded
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           Device number for output.
*
*  MAXERR   (input) INTEGER
*           Max number of errors that can be stored in ERRIBUFF or
*           ERRZBUFF
*
*  NERR     (output) INTEGER
*           The number of errors that have been found.
*
*  ERRIBUF  (output) INTEGER array, dimension (6,MAXERRS)
*           Buffer in which to store integer error information.  It will
*           be built up in the following format for the call to TSEND.
*           All integer information is recorded in the following 6-tuple
*           {TESTNUM, SRC, DEST, I, J, WHAT}. These values are figured:
*             SRC = RSRC * NPROCS + CSRC
*             DEST = RDEST * NPROCS + CDEST
*             WHAT
*              = 1 : Error in pre-padding
*              = 2 : Error in post-padding
*              = 3 : Error in LDA-M gap
*              = 4 : Error in complementory triangle
*              ELSE: Error in matrix
*           If there are more errors than can fit in the error buffer,
*           the error number will indicate the actual number of errors
*           found, but the buffer will be truncated to the maximum
*           number of errors which can fit.
*
*  ERRDBUF  (output) @(doctype) array, dimension (2, MAXERRS)
*           Buffer in which to store error data information.
*           {Incorrect, Predicted}
*
*  TFAILED (input/ourput) INTEGER array, dimension NTESTS
*          Workspace used to keep track of which tests failed.
*          This array not accessed unless COUNTING is true.
*
*  ===================================================================
*
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      EXTERNAL IBTMYPROC, IBTNPROCS
*     ..
*     .. Local Scalars ..
      CHARACTER*1 MAT
      LOGICAL MATISINT
      INTEGER OLDTEST, NPROCS, PROW, PCOL, I, ERRTYPE
*     ..
*     .. Executable Statements ..
*
      IF( (IBTMYPROC().NE.0) .OR. (NERR.LE.0) ) RETURN
      OLDTEST = -1
      NPROCS = IBTNPROCS()
      PROW = ERRIBUF(3,1) / NPROCS
      PCOL = MOD( ERRIBUF(3,1), NPROCS )
      IF( NERR .GT. MAXERR ) WRITE(OUTNUM,13000)
*
      DO 20 I = 1, MIN( NERR, MAXERR )
         IF( ERRIBUF(1,I) .NE. OLDTEST ) THEN
            IF( OLDTEST .NE. -1 )
     $         WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM,1000) PROW, PCOL, ERRIBUF(1,I)
            IF( COUNTING ) TFAILED( ERRIBUF(1,I) ) = 1
            OLDTEST = ERRIBUF(1, I)
         END IF
*
*        Print out error message depending on type of error
*
         ERRTYPE = ERRIBUF(6, I)
         IF( ERRTYPE .LT. -10 ) THEN
            ERRTYPE = -ERRTYPE - 10
            MAT = 'C'
            MATISINT = .TRUE.
         ELSE IF( ERRTYPE .LT. 0 ) THEN
            ERRTYPE = -ERRTYPE
            MAT = 'R'
            MATISINT = .TRUE.
         ELSE
            MATISINT = .FALSE.
         END IF
*
*        RA/CA arrays from MAX/MIN have different printing protocol
*
         IF( MATISINT ) THEN
            IF( ERRIBUF(2, I) .EQ. -1 ) THEN
               WRITE(OUTNUM,11000) ERRIBUF(4,I), ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,7000) ERRIBUF(5,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,8000) ERRIBUF(4,I), MAT,
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,9000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $            INT( ERRDBUF(2,I) ), INT( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,10000) MAT, ERRIBUF(4,I), ERRIBUF(5,I),
     $                             INT( ERRDBUF(2,I) ),
     $                             INT( ERRDBUF(1,I) )
            END IF
*
*        Have memory overwrites in matrix A
*
         ELSE
            IF( ERRTYPE .EQ. ERR_PRE ) THEN
               WRITE(OUTNUM,2000) ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), DIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), DIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_POST ) THEN
               WRITE(OUTNUM,3000) ERRIBUF(4,I),
     $         REAL( ERRDBUF(2,I) ), DIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), DIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_GAP ) THEN
               WRITE(OUTNUM,4000)
     $         ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), DIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), DIMAG( ERRDBUF(1,I) )
            ELSE IF( ERRTYPE .EQ. ERR_TRI ) THEN
               WRITE(OUTNUM,5000) ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), DIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), DIMAG( ERRDBUF(1,I) )
            ELSE
               WRITE(OUTNUM,6000) ERRIBUF(4,I), ERRIBUF(5,I),
     $         REAL( ERRDBUF(2,I) ), DIMAG( ERRDBUF(2,I) ),
     $         REAL( ERRDBUF(1,I) ), DIMAG( ERRDBUF(1,I) )
            END IF
         END IF
   20 CONTINUE
      WRITE(OUTNUM,12000) PROW, PCOL, OLDTEST
*
 1000 FORMAT('PROCESS {',I4,',',I4,'} REPORTS ERRORS IN TEST#',I6,':')
 2000 FORMAT('   Buffer overwrite ',I4,
     $       ' elements before the start of A:',/,
     $       '   Expected=','[',G22.15,',',G22.15,']',
     $       '; Received=','[',G22.15,',',G22.15,']')
 3000 FORMAT('   Buffer overwrite ',I4,' elements after the end of A:',
     $       /,'   Expected=','[',G22.15,',',G22.15,']',
     $       '; Received=','[',G22.15,',',G22.15,']')
 4000 FORMAT('   LDA-M gap overwrite at postion (',I4,',',I4,'):',/,
     $       '   Expected=','[',G22.15,',',G22.15,']',
     $       '; Received=','[',G22.15,',',G22.15,']')
 5000 FORMAT('   Complementory triangle overwrite at A(',I4,',',I4,
     $       '):',/,'   Expected=','[',G22.15,',',G22.15,']',
     $       '; Received=','[',G22.15,',',G22.15,']')
 6000 FORMAT('   Invalid element at A(',I4,',',I4,'):',/,
     $       '   Expected=','[',G22.15,',',G22.15,']',
     $       '; Received=','[',G22.15,',',G22.15,']')
 7000 FORMAT('   Buffer overwrite ',I4,' elements before the start of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
 8000 FORMAT('   Buffer overwrite ',I4,' elements after the end of ',
     $       A1,'A:',/,'   Expected=',I12,'; Received=',I12)
*
 9000 FORMAT('   LD',A1,'A-M gap overwrite at postion (',I4,',',I4,'):'
     $       ,/,'   Expected=',I12,'; Received=',I12)
*
10000 FORMAT('   Invalid element at ',A1,'A(',I4,',',I4,'):',/,
     $       '   Expected=',I12,'; Received=',I12)
11000 FORMAT('   Overwrite at position (',I4,',',I4,') of non-existent '
     $       ,A1,'A array.',/,'   Expected=',I12,'; Received=',I12)
12000 FORMAT('PROCESS {',I4,',',I4,'} DONE ERROR REPORT FOR TEST#',
     $       I6,'.')
13000 FORMAT('WARNING: There were more errors than could be recorded.',
     $       /,'Increase MEMELTS to get complete listing.')
      RETURN
*
*     End ZPRINTERRS
*
      END
*
*
      SUBROUTINE ISUMTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*)
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ITESTSUM:  Test integer SUM COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  MEM      (workspace) INTEGER array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, IGSUM2D
      EXTERNAL IINITMAT, ICHKPAD, IBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I, IAM,
     $        IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC, ISIZE, ISTART,
     $        ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1, ITR2, J, K, LDA,
     $        LDADST, LDASRC, M, MAXERR, MYCOL, MYROW, N, NERR, NFAIL,
     $        NPCOL, NPROW, NSKIP, PREAPTR, RDEST, RDEST2, SETWHAT,
     $        TESTNUM
      INTEGER CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -911
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      MAXERR = ( ISIZE * (MEMLEN-I) ) / ( ISIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SUM tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL IINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              CALL IGSUM2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RDEST2,
     $                                     CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ICHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ICHKSUM(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED)
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL IBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL IBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('INTEGER SUM TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,6I6,2I5)
 7000 FORMAT('INTEGER SUM TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('INTEGER SUM TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ITESTSUM.
*
      END
*
      INTEGER FUNCTION IBTABS(VAL)
      INTEGER VAL
      IBTABS = ABS(VAL)
      RETURN
      END
*
      SUBROUTINE ICHKSUM( SCOPE, ICTXT, M, N, A, LDA, TESTNUM, MAXERR,
     $                    NERR, ERRIBUF, ERRDBUF, ISEED )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), ISEED(*)
      INTEGER A(LDA,*), ERRDBUF(2, MAXERR)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      INTEGER IBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTRAN
*     ..
*     .. Local Scalars ..
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, NODE, NNODES, DEST
      INTEGER I, J, K
      INTEGER ANS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            ANS = 0
            DO 40 K = 0, NNODES-1
               ANS = ANS + IBTRAN( ISEED(K*4+1) )
   40       CONTINUE
*
*           The error bound is figured by
*           2 * eps * (nnodes-1) * max(|max element|, |ans|).
*           The 2 allows for errors in the distributed _AND_ local result.
*           The eps is machine epsilon.  The number of floating point adds
*           is (nnodes - 1).  We use the fact that 0.5 is the maximum element
*           in order to save ourselves some computation.
*
            IF( ANS .NE. A(I,J) ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = ANS
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ICHKSUM
*
      END
*
*
      SUBROUTINE SSUMTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*)
      REAL MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  STESTSUM:  Test real SUM COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  MEM      (workspace) REAL array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, SGSUM2D
      EXTERNAL SINITMAT, SCHKPAD, SBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I, IAM,
     $        IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC, ISIZE, ISTART,
     $        ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1, ITR2, J, K, LDA,
     $        LDADST, LDASRC, M, MAXERR, MYCOL, MYROW, N, NERR, NFAIL,
     $        NPCOL, NPROW, NSKIP, PREAPTR, RDEST, RDEST2, SETWHAT,
     $        SSIZE, TESTNUM
      REAL CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.61E0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      MAXERR = ( SSIZE * (MEMLEN-I) ) / ( SSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SUM tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL SINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              CALL SGSUM2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RDEST2,
     $                                     CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL SCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL SCHKSUM(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED)
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL SBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL SBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('REAL SUM TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,6I6,2I5)
 7000 FORMAT('REAL SUM TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('REAL SUM TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of STESTSUM.
*
      END
*
      REAL FUNCTION SBTABS(VAL)
      REAL VAL
      SBTABS = ABS(VAL)
      RETURN
      END
*
      REAL FUNCTION SBTEPS()
*
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTMSGID
      REAL SLAMCH
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID, SLAMCH
*     ..
*     .. Local Scalars ..
      INTEGER I, IAM, NNODES
      REAL EPS, EPS2
      SAVE EPS
      DATA EPS /-22.0E0/
*     ..
*     .. Executable Statements ..
*
*     First time called, must get max epsilon possessed by any
*     participating process
*
      IF( EPS .EQ. -22.0E0 ) THEN
         IAM = IBTMYPROC()
         NNODES = IBTNPROCS()
         EPS = SLAMCH('epsilon')
         IF( IAM .EQ. 0 ) THEN
            IF( NNODES .GT. 1 ) THEN
               DO 10 I = 1, NNODES-1
                  CALL BTRECV( 4, 1, EPS2, I, IBTMSGID()+20 )
                  IF( EPS .LT. EPS2 ) EPS = EPS2
   10          CONTINUE
            END IF
            CALL BTSEND( 4, 1, EPS, -1, IBTMSGID()+20 )
         ELSE
            CALL BTSEND( 4, 1, EPS, 0, IBTMSGID()+20 )
            CALL BTRECV( 4, 1, EPS, 0, IBTMSGID()+20 )
         ENDIF
      END IF
      SBTEPS = EPS
      RETURN
*
*     End SBTEPS
*
      END
*
      SUBROUTINE SCHKSUM( SCOPE, ICTXT, M, N, A, LDA, TESTNUM, MAXERR,
     $                    NERR, ERRIBUF, ERRDBUF, ISEED )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), ISEED(*)
      REAL A(LDA,*), ERRDBUF(2, MAXERR)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      REAL SBTEPS
      REAL SBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, SBTEPS, SBTRAN
*     ..
*     .. Local Scalars ..
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, NODE, NNODES, DEST
      INTEGER I, J, K
      REAL ANS, EPS, ERRBND, POSNUM, NEGNUM, TMP
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            ANS = 0
            POSNUM = 0
            NEGNUM = 0
            DO 40 K = 0, NNODES-1
               TMP = SBTRAN( ISEED(K*4+1) )
               IF( TMP .LT. 0 ) THEN
                  NEGNUM = NEGNUM + TMP
               ELSE
                  POSNUM = POSNUM + TMP
               END IF
               ANS = ANS + TMP
   40       CONTINUE
*
*           The error bound is figured by
*           2 * eps * (nnodes-1) * max(|max element|, |ans|).
*           The 2 allows for errors in the distributed _AND_ local result.
*           The eps is machine epsilon.  The number of floating point adds
*           is (nnodes - 1).  We use the fact that 0.5 is the maximum element
*           in order to save ourselves some computation.
*
            ERRBND = 2 * EPS * NNODES * MAX( POSNUM, -NEGNUM )
            IF( ABS( ANS - A(I,J) ) .GT. ERRBND ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = ANS
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of SCHKSUM
*
      END
*
*
      SUBROUTINE DSUMTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*)
      DOUBLE PRECISION MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  DTESTSUM:  Test double precision SUM COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  MEM      (workspace) DOUBLE PRECISION array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, DGSUM2D
      EXTERNAL DINITMAT, DCHKPAD, DBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CDEST, CDEST2, CONTEXT, DSIZE, ERRDPTR, ERRIPTR, I,
     $        IAM, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC, ISIZE, ISTART,
     $        ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1, ITR2, J, K, LDA,
     $        LDADST, LDASRC, M, MAXERR, MYCOL, MYROW, N, NERR, NFAIL,
     $        NPCOL, NPROW, NSKIP, PREAPTR, RDEST, RDEST2, SETWHAT,
     $        TESTNUM
      DOUBLE PRECISION CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.81D0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      MAXERR = ( DSIZE * (MEMLEN-I) ) / ( DSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SUM tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL DINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              CALL DGSUM2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RDEST2,
     $                                     CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL DCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL DCHKSUM(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED)
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL DBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL DBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE PRECISION SUM TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,6I6,2I5)
 7000 FORMAT('DOUBLE PRECISION SUM TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE PRECISION SUM TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of DTESTSUM.
*
      END
*
      DOUBLE PRECISION FUNCTION DBTABS(VAL)
      DOUBLE PRECISION VAL
      DBTABS = ABS(VAL)
      RETURN
      END
*
      DOUBLE PRECISION FUNCTION DBTEPS()
*
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTMSGID
      DOUBLE PRECISION DLAMCH
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTMSGID, DLAMCH
*     ..
*     .. Local Scalars ..
      INTEGER I, IAM, NNODES
      DOUBLE PRECISION EPS, EPS2
      SAVE EPS
      DATA EPS /-22.0D0/
*     ..
*     .. Executable Statements ..
*
*     First time called, must get max epsilon possessed by any
*     participating process
*
      IF( EPS .EQ. -22.0D0 ) THEN
         IAM = IBTMYPROC()
         NNODES = IBTNPROCS()
         EPS = DLAMCH('epsilon')
         IF( IAM .EQ. 0 ) THEN
            IF( NNODES .GT. 1 ) THEN
               DO 10 I = 1, NNODES-1
                  CALL BTRECV( 6, 1, EPS2, I, IBTMSGID()+20 )
                  IF( EPS .LT. EPS2 ) EPS = EPS2
   10          CONTINUE
            END IF
            CALL BTSEND( 6, 1, EPS, -1, IBTMSGID()+20 )
         ELSE
            CALL BTSEND( 6, 1, EPS, 0, IBTMSGID()+20 )
            CALL BTRECV( 6, 1, EPS, 0, IBTMSGID()+20 )
         ENDIF
      END IF
      DBTEPS = EPS
      RETURN
*
*     End DBTEPS
*
      END
*
      SUBROUTINE DCHKSUM( SCOPE, ICTXT, M, N, A, LDA, TESTNUM, MAXERR,
     $                    NERR, ERRIBUF, ERRDBUF, ISEED )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE PRECISION A(LDA,*), ERRDBUF(2, MAXERR)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      DOUBLE PRECISION DBTEPS
      DOUBLE PRECISION DBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, DBTEPS, DBTRAN
*     ..
*     .. Local Scalars ..
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, NODE, NNODES, DEST
      INTEGER I, J, K
      DOUBLE PRECISION ANS, EPS, ERRBND, POSNUM, NEGNUM, TMP
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            ANS = 0
            POSNUM = 0
            NEGNUM = 0
            DO 40 K = 0, NNODES-1
               TMP = DBTRAN( ISEED(K*4+1) )
               IF( TMP .LT. 0 ) THEN
                  NEGNUM = NEGNUM + TMP
               ELSE
                  POSNUM = POSNUM + TMP
               END IF
               ANS = ANS + TMP
   40       CONTINUE
*
*           The error bound is figured by
*           2 * eps * (nnodes-1) * max(|max element|, |ans|).
*           The 2 allows for errors in the distributed _AND_ local result.
*           The eps is machine epsilon.  The number of floating point adds
*           is (nnodes - 1).  We use the fact that 0.5 is the maximum element
*           in order to save ourselves some computation.
*
            ERRBND = 2 * EPS * NNODES * MAX( POSNUM, -NEGNUM )
            IF( ABS( ANS - A(I,J) ) .GT. ERRBND ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = ANS
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of DCHKSUM
*
      END
*
*
      SUBROUTINE CSUMTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*)
      COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  CTESTSUM:  Test complex SUM COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  MEM      (workspace) COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, CGSUM2D
      EXTERNAL CINITMAT, CCHKPAD, CBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CDEST, CDEST2, CONTEXT, CSIZE, ERRDPTR, ERRIPTR, I,
     $        IAM, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC, ISIZE, ISTART,
     $        ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1, ITR2, J, K, LDA,
     $        LDADST, LDASRC, M, MAXERR, MYCOL, MYROW, N, NERR, NFAIL,
     $        NPCOL, NPROW, NSKIP, PREAPTR, RDEST, RDEST2, SETWHAT,
     $        TESTNUM
      COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = CMPLX( -0.91E0, -0.71E0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      CSIZE = IBTSIZEOF('C')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      MAXERR = ( CSIZE * (MEMLEN-I) ) / ( CSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SUM tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL CINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              CALL CGSUM2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RDEST2,
     $                                     CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL CCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL CCHKSUM(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED)
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL CBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL CBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('COMPLEX SUM TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,6I6,2I5)
 7000 FORMAT('COMPLEX SUM TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('COMPLEX SUM TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of CTESTSUM.
*
      END
*
      REAL FUNCTION CBTABS(VAL)
      COMPLEX VAL
      CBTABS = ABS( REAL(VAL) ) + ABS( AIMAG(VAL) )
      RETURN
      END
*
      SUBROUTINE CCHKSUM( SCOPE, ICTXT, M, N, A, LDA, TESTNUM, MAXERR,
     $                    NERR, ERRIBUF, ERRDBUF, ISEED )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), ISEED(*)
      COMPLEX A(LDA,*), ERRDBUF(2, MAXERR)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      REAL SBTEPS
      COMPLEX CBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, SBTEPS, CBTRAN
*     ..
*     .. Local Scalars ..
      LOGICAL NUMOK
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, NODE, NNODES, DEST
      INTEGER I, J, K
      COMPLEX ANS, TMP
      REAL EPS, ERRBND, RPOSNUM, RNEGNUM, IPOSNUM, INEGNUM
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            ANS = 0
            RPOSNUM = 0
            RNEGNUM = 0
            IPOSNUM = 0
            INEGNUM = 0
            DO 40 K = 0, NNODES-1
               TMP = CBTRAN( ISEED(K*4+1) )
               IF( REAL( TMP ) .LT. 0 ) THEN
                  RNEGNUM = RNEGNUM + REAL( TMP )
               ELSE
                  RPOSNUM = RPOSNUM + REAL( TMP )
               END IF
               IF( AIMAG( TMP ) .LT. 0 ) THEN
                  INEGNUM = INEGNUM + AIMAG( TMP )
               ELSE
                  IPOSNUM = IPOSNUM + AIMAG( TMP )
               END IF
               ANS = ANS + TMP
   40       CONTINUE
*
*           The error bound is figured by
*           2 * eps * (nnodes-1) * max(|max element|, |ans|).
*           The 2 allows for errors in the distributed _AND_ local result.
*           The eps is machine epsilon.  The number of floating point adds
*           is (nnodes - 1).  We use the fact that 0.5 is the maximum element
*           in order to save ourselves some computation.
*
            TMP = ANS - A(I,J)
            ERRBND = 2 * EPS * NNODES * MAX( RPOSNUM, -RNEGNUM )
            NUMOK = ( REAL(TMP) .LE. ERRBND )
            ERRBND = 2 * EPS * NNODES * MAX( IPOSNUM, -INEGNUM )
            NUMOK = NUMOK .AND. ( AIMAG(TMP) .LE. ERRBND )
            IF( .NOT.NUMOK ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = ANS
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of CCHKSUM
*
      END
*
*
      SUBROUTINE ZSUMTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*)
      DOUBLE COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ZTESTSUM:  Test double complex SUM COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  MEM      (workspace) DOUBLE COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, ZGSUM2D
      EXTERNAL ZINITMAT, ZCHKPAD, ZBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I, IAM,
     $        IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC, ISIZE, ISTART,
     $        ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1, ITR2, J, K, LDA,
     $        LDADST, LDASRC, M, MAXERR, MYCOL, MYROW, N, NERR, NFAIL,
     $        NPCOL, NPROW, NSKIP, PREAPTR, RDEST, RDEST2, SETWHAT,
     $        TESTNUM, ZSIZE
      DOUBLE COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = DCMPLX( -9.11D0, -9.21D0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      ZSIZE = IBTSIZEOF('Z')
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      MAXERR = ( ZSIZE * (MEMLEN-I) ) / ( ZSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run SUM tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL ZINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
                              CALL ZGSUM2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RDEST2,
     $                                     CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ZCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ZCHKSUM(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED)
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL ZBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL ZBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE COMPLEX SUM TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,6I6,2I5)
 7000 FORMAT('DOUBLE COMPLEX SUM TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE COMPLEX SUM TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ZTESTSUM.
*
      END
*
      DOUBLE PRECISION FUNCTION ZBTABS(VAL)
      DOUBLE COMPLEX VAL
      ZBTABS = ABS( DBLE(VAL) ) + ABS( DIMAG(VAL) )
      RETURN
      END
*
      SUBROUTINE ZCHKSUM( SCOPE, ICTXT, M, N, A, LDA, TESTNUM, MAXERR,
     $                    NERR, ERRIBUF, ERRDBUF, ISEED )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE COMPLEX A(LDA,*), ERRDBUF(2, MAXERR)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS
      DOUBLE PRECISION DBTEPS
      DOUBLE COMPLEX ZBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, DBTEPS, ZBTRAN
*     ..
*     .. Local Scalars ..
      LOGICAL NUMOK
      INTEGER NPROCS, NPROW, NPCOL, MYROW, MYCOL, NODE, NNODES, DEST
      INTEGER I, J, K
      DOUBLE COMPLEX ANS, TMP
      DOUBLE PRECISION EPS, ERRBND, RPOSNUM, RNEGNUM, IPOSNUM, INEGNUM
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            ANS = 0
            RPOSNUM = 0
            RNEGNUM = 0
            IPOSNUM = 0
            INEGNUM = 0
            DO 40 K = 0, NNODES-1
               TMP = ZBTRAN( ISEED(K*4+1) )
               IF( DBLE( TMP ) .LT. 0 ) THEN
                  RNEGNUM = RNEGNUM + DBLE( TMP )
               ELSE
                  RPOSNUM = RPOSNUM + DBLE( TMP )
               END IF
               IF( DIMAG( TMP ) .LT. 0 ) THEN
                  INEGNUM = INEGNUM + DIMAG( TMP )
               ELSE
                  IPOSNUM = IPOSNUM + DIMAG( TMP )
               END IF
               ANS = ANS + TMP
   40       CONTINUE
*
*           The error bound is figured by
*           2 * eps * (nnodes-1) * max(|max element|, |ans|).
*           The 2 allows for errors in the distributed _AND_ local result.
*           The eps is machine epsilon.  The number of floating point adds
*           is (nnodes - 1).  We use the fact that 0.5 is the maximum element
*           in order to save ourselves some computation.
*
            TMP = ANS - A(I,J)
            ERRBND = 2 * EPS * NNODES * MAX( RPOSNUM, -RNEGNUM )
            NUMOK = ( DBLE(TMP) .LE. ERRBND )
            ERRBND = 2 * EPS * NNODES * MAX( IPOSNUM, -INEGNUM )
            NUMOK = NUMOK .AND. ( DIMAG(TMP) .LE. ERRBND )
            IF( .NOT.NUMOK ) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = ANS
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ZCHKSUM
*
      END
*
*
      SUBROUTINE IAMXTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ITESTAMX:  Test integer AMX COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) INTEGER array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, IGAMX2D
      EXTERNAL IINITMAT, ICHKPAD, IBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      INTEGER CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -911
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( ISIZE * (MEMLEN-I) ) / ( ISIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MAX tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL IINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL IGAMX2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ICHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ICHKAMX(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL IRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL IBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL IBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('INTEGER AMX TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('INTEGER AMX TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('INTEGER AMX TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ITESTAMX.
*
      END
*
      SUBROUTINE IBTSPCOORD( SCOPE, PNUM, MYROW, MYCOL, NPCOL,
     $                       PROW, PCOL )
      CHARACTER*1 SCOPE
      INTEGER PNUM, MYROW, MYCOL, NPCOL, PROW, PCOL
*
      IF( SCOPE .EQ. 'R' ) THEN
         PROW = MYROW
         PCOL = PNUM
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         PROW = PNUM
         PCOL = MYCOL
      ELSE
         PROW = PNUM / NPCOL
         PCOL = MOD( PNUM, NPCOL )
      END IF
      RETURN
*
*     End of ibtspcoord
*
      END
*
      INTEGER FUNCTION IBTSPNUM( SCOPE, PROW, PCOL, NPCOL )
      CHARACTER*1 SCOPE
      INTEGER PROW, PCOL, NPCOL
      IF( SCOPE .EQ. 'R' ) THEN
         IBTSPNUM = PCOL
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         IBTSPNUM = PROW
      ELSE
         IBTSPNUM = PROW*NPCOL + PCOL
      END IF
*
      RETURN
*
*     End of ibtscpnum
*
      END
*
      SUBROUTINE IRCCHK( IPRE, IPOST, PADVAL, M, N, RA, CA, LDI, MYROW,
     $                   MYCOL, TESTNUM, MAXERR, NERR,
     $                   ERRIBUF, ERRDBUF )
*
*     .. Scalar Arguments ..
      INTEGER IPRE, IPOST, PADVAL, M, N, LDI, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR)
      INTEGER ERRDBUF(2, MAXERR)
*     ..
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, K, IAM
*     ..
*     .. Executable Statements ..
*
      IAM = MYROW * IBTNPROCS() + MYCOL
*
*     Check pre padding
*
      IF( LDI .NE. -1 ) THEN
         IF( IPRE .GT. 0 ) THEN
            DO 10 I = 1, IPRE
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -ERR_PRE
                     ERRDBUF(1, NERR) = INT( RA(I) )
                     ERRDBUF(2, NERR) = INT( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -10 - ERR_PRE
                     ERRDBUF(1, NERR) = INT( CA(I) )
                     ERRDBUF(2, NERR) = INT( PADVAL )
                  END IF
               ENDIF
   10       CONTINUE
         END IF
*
*        Check post padding
*
         IF( IPOST .GT. 0 ) THEN
            K = IPRE + LDI*N
            DO 20 I = K+1, K+IPOST
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -ERR_POST
                     ERRDBUF(1, NERR) = INT( RA(I) )
                     ERRDBUF(2, NERR) = INT( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -10 - ERR_POST
                     ERRDBUF(1, NERR) = INT( CA(I) )
                     ERRDBUF(2, NERR) = INT( PADVAL )
                  END IF
               ENDIF
   20       CONTINUE
         END IF
*
*        Check all (LDI-M) gaps
*
         IF( LDI .GT. M ) THEN
            K = IPRE + M + 1
            DO 40 J = 1, N
               DO 30 I = M+1, LDI
                  K = IPRE + (J-1)*LDI + I
                  IF( RA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -ERR_GAP
                        ERRDBUF(1, NERR) = INT( RA(K) )
                        ERRDBUF(2, NERR) = INT( PADVAL )
                     END IF
                  END IF
                  IF( CA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -10 - ERR_GAP
                        ERRDBUF(1, NERR) = INT( CA(K) )
                        ERRDBUF(2, NERR) = INT( PADVAL )
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
*
*     if RA and CA don't exist, buffs better be untouched
*
      ELSE
         DO 50 I = 1, IPRE+IPOST
            IF( RA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -ERR_PRE
                  ERRDBUF(1, NERR) = INT( RA(I) )
                  ERRDBUF(2, NERR) = INT( PADVAL )
               END IF
            END IF
            IF( CA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -10 - ERR_PRE
                  ERRDBUF(1, NERR) = INT( CA(I) )
                  ERRDBUF(2, NERR) = INT( PADVAL )
               END IF
            END IF
   50    CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE ICHKAMX( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      INTEGER A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSPNUM, IBTRAN, IBTABS
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, IBTRAN
      EXTERNAL IBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMX, CAMX
      INTEGER IAMX, I, J, K, H, DEST, NODE
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = IBTRAN( ISEED )
            IAMX = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  IBTRAN( ISEED(K*4+1) )
                  IF( IBTABS( VALS(K+1) ) .GT. IBTABS( VALS(IAMX) ) )
     $               IAMX = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMX) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = IBTABS( VALS(K) ).NE.IBTABS( VALS(IAMX) )
                     IF( .NOT.ERROR ) IAMX = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( IBTABS( A(I,J) ) .NE. IBTABS( VALS(IAMX) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMX)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMX ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMX) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMX-1, MYROW, MYCOL,
     $                                NPCOL, RAMX, CAMX )
                     IF( RAMX .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMX
                     END IF
                     IF( CAMX .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMX
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ICHKAMX
*
      END
*
*
      SUBROUTINE SAMXTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      REAL MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  STESTAMX:  Test real AMX COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) REAL array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, SGAMX2D
      EXTERNAL SINITMAT, SCHKPAD, SBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, SSIZE, TESTNUM, VALPTR
      REAL CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.61E0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( SSIZE * (MEMLEN-I) ) / ( SSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MAX tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL SINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL SGAMX2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL SCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL SCHKAMX(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL SRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL SBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL SBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('REAL AMX TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('REAL AMX TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('REAL AMX TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of STESTAMX.
*
      END
*
      SUBROUTINE SRCCHK( IPRE, IPOST, PADVAL, M, N, RA, CA, LDI, MYROW,
     $                   MYCOL, TESTNUM, MAXERR, NERR,
     $                   ERRIBUF, ERRDBUF )
*
*     .. Scalar Arguments ..
      INTEGER IPRE, IPOST, PADVAL, M, N, LDI, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR)
      REAL ERRDBUF(2, MAXERR)
*     ..
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, K, IAM
*     ..
*     .. Executable Statements ..
*
      IAM = MYROW * IBTNPROCS() + MYCOL
*
*     Check pre padding
*
      IF( LDI .NE. -1 ) THEN
         IF( IPRE .GT. 0 ) THEN
            DO 10 I = 1, IPRE
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -ERR_PRE
                     ERRDBUF(1, NERR) = REAL( RA(I) )
                     ERRDBUF(2, NERR) = REAL( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -10 - ERR_PRE
                     ERRDBUF(1, NERR) = REAL( CA(I) )
                     ERRDBUF(2, NERR) = REAL( PADVAL )
                  END IF
               ENDIF
   10       CONTINUE
         END IF
*
*        Check post padding
*
         IF( IPOST .GT. 0 ) THEN
            K = IPRE + LDI*N
            DO 20 I = K+1, K+IPOST
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -ERR_POST
                     ERRDBUF(1, NERR) = REAL( RA(I) )
                     ERRDBUF(2, NERR) = REAL( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -10 - ERR_POST
                     ERRDBUF(1, NERR) = REAL( CA(I) )
                     ERRDBUF(2, NERR) = REAL( PADVAL )
                  END IF
               ENDIF
   20       CONTINUE
         END IF
*
*        Check all (LDI-M) gaps
*
         IF( LDI .GT. M ) THEN
            K = IPRE + M + 1
            DO 40 J = 1, N
               DO 30 I = M+1, LDI
                  K = IPRE + (J-1)*LDI + I
                  IF( RA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -ERR_GAP
                        ERRDBUF(1, NERR) = REAL( RA(K) )
                        ERRDBUF(2, NERR) = REAL( PADVAL )
                     END IF
                  END IF
                  IF( CA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -10 - ERR_GAP
                        ERRDBUF(1, NERR) = REAL( CA(K) )
                        ERRDBUF(2, NERR) = REAL( PADVAL )
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
*
*     if RA and CA don't exist, buffs better be untouched
*
      ELSE
         DO 50 I = 1, IPRE+IPOST
            IF( RA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -ERR_PRE
                  ERRDBUF(1, NERR) = REAL( RA(I) )
                  ERRDBUF(2, NERR) = REAL( PADVAL )
               END IF
            END IF
            IF( CA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -10 - ERR_PRE
                  ERRDBUF(1, NERR) = REAL( CA(I) )
                  ERRDBUF(2, NERR) = REAL( PADVAL )
               END IF
            END IF
   50    CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE SCHKAMX( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      REAL A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      REAL SBTEPS, SBTABS
      REAL SBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, SBTRAN, SBTEPS, SBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMX, CAMX
      INTEGER IAMX, I, J, K, H, DEST, NODE
      REAL EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = SBTRAN( ISEED )
            IAMX = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  SBTRAN( ISEED(K*4+1) )
                  IF( SBTABS( VALS(K+1) ) .GT. SBTABS( VALS(IAMX) ) )
     $               IAMX = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMX) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = SBTABS( VALS(K) ).NE.SBTABS( VALS(IAMX) )
                     IF( .NOT.ERROR ) IAMX = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( SBTABS( A(I,J) ) .NE. SBTABS( VALS(IAMX) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMX)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMX ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMX) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMX-1, MYROW, MYCOL,
     $                                NPCOL, RAMX, CAMX )
                     IF( RAMX .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMX
                     END IF
                     IF( CAMX .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMX
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of SCHKAMX
*
      END
*
*
      SUBROUTINE DAMXTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      DOUBLE PRECISION MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  DTESTAMX:  Test double precision AMX COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) DOUBLE PRECISION array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, DGAMX2D
      EXTERNAL DINITMAT, DCHKPAD, DBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, DSIZE, ERRDPTR,
     $        ERRIPTR, I, IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST,
     $        IPRE, ISC, ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO,
     $        ITR, ITR1, ITR2, J, K, LDA, LDADST, LDASRC, LDI, M,
     $        MAXERR, MYCOL, MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP,
     $        PREAPTR, RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      DOUBLE PRECISION CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.81D0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( DSIZE * (MEMLEN-I) ) / ( DSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MAX tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL DINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL DGAMX2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL DCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL DCHKAMX(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL DRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL DBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL DBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE PRECISION AMX TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('DOUBLE PRECISION AMX TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE PRECISION AMX TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of DTESTAMX.
*
      END
*
      SUBROUTINE DRCCHK( IPRE, IPOST, PADVAL, M, N, RA, CA, LDI, MYROW,
     $                   MYCOL, TESTNUM, MAXERR, NERR,
     $                   ERRIBUF, ERRDBUF )
*
*     .. Scalar Arguments ..
      INTEGER IPRE, IPOST, PADVAL, M, N, LDI, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR)
      DOUBLE PRECISION ERRDBUF(2, MAXERR)
*     ..
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, K, IAM
*     ..
*     .. Executable Statements ..
*
      IAM = MYROW * IBTNPROCS() + MYCOL
*
*     Check pre padding
*
      IF( LDI .NE. -1 ) THEN
         IF( IPRE .GT. 0 ) THEN
            DO 10 I = 1, IPRE
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -ERR_PRE
                     ERRDBUF(1, NERR) = DBLE( RA(I) )
                     ERRDBUF(2, NERR) = DBLE( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -10 - ERR_PRE
                     ERRDBUF(1, NERR) = DBLE( CA(I) )
                     ERRDBUF(2, NERR) = DBLE( PADVAL )
                  END IF
               ENDIF
   10       CONTINUE
         END IF
*
*        Check post padding
*
         IF( IPOST .GT. 0 ) THEN
            K = IPRE + LDI*N
            DO 20 I = K+1, K+IPOST
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -ERR_POST
                     ERRDBUF(1, NERR) = DBLE( RA(I) )
                     ERRDBUF(2, NERR) = DBLE( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -10 - ERR_POST
                     ERRDBUF(1, NERR) = DBLE( CA(I) )
                     ERRDBUF(2, NERR) = DBLE( PADVAL )
                  END IF
               ENDIF
   20       CONTINUE
         END IF
*
*        Check all (LDI-M) gaps
*
         IF( LDI .GT. M ) THEN
            K = IPRE + M + 1
            DO 40 J = 1, N
               DO 30 I = M+1, LDI
                  K = IPRE + (J-1)*LDI + I
                  IF( RA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -ERR_GAP
                        ERRDBUF(1, NERR) = DBLE( RA(K) )
                        ERRDBUF(2, NERR) = DBLE( PADVAL )
                     END IF
                  END IF
                  IF( CA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -10 - ERR_GAP
                        ERRDBUF(1, NERR) = DBLE( CA(K) )
                        ERRDBUF(2, NERR) = DBLE( PADVAL )
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
*
*     if RA and CA don't exist, buffs better be untouched
*
      ELSE
         DO 50 I = 1, IPRE+IPOST
            IF( RA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -ERR_PRE
                  ERRDBUF(1, NERR) = DBLE( RA(I) )
                  ERRDBUF(2, NERR) = DBLE( PADVAL )
               END IF
            END IF
            IF( CA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -10 - ERR_PRE
                  ERRDBUF(1, NERR) = DBLE( CA(I) )
                  ERRDBUF(2, NERR) = DBLE( PADVAL )
               END IF
            END IF
   50    CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE DCHKAMX( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE PRECISION A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      DOUBLE PRECISION DBTEPS, DBTABS
      DOUBLE PRECISION DBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, DBTRAN, DBTEPS, DBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMX, CAMX
      INTEGER IAMX, I, J, K, H, DEST, NODE
      DOUBLE PRECISION EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = DBTRAN( ISEED )
            IAMX = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  DBTRAN( ISEED(K*4+1) )
                  IF( DBTABS( VALS(K+1) ) .GT. DBTABS( VALS(IAMX) ) )
     $               IAMX = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMX) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = DBTABS( VALS(K) ).NE.DBTABS( VALS(IAMX) )
                     IF( .NOT.ERROR ) IAMX = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( DBTABS( A(I,J) ) .NE. DBTABS( VALS(IAMX) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMX)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMX ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMX) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMX-1, MYROW, MYCOL,
     $                                NPCOL, RAMX, CAMX )
                     IF( RAMX .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMX
                     END IF
                     IF( CAMX .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMX
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of DCHKAMX
*
      END
*
*
      SUBROUTINE CAMXTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  CTESTAMX:  Test complex AMX COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, CGAMX2D
      EXTERNAL CINITMAT, CCHKPAD, CBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, CSIZE, ERRDPTR,
     $        ERRIPTR, I, IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST,
     $        IPRE, ISC, ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO,
     $        ITR, ITR1, ITR2, J, K, LDA, LDADST, LDASRC, LDI, M,
     $        MAXERR, MYCOL, MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP,
     $        PREAPTR, RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = CMPLX( -0.91E0, -0.71E0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      CSIZE = IBTSIZEOF('C')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( CSIZE * (MEMLEN-I) ) / ( CSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MAX tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL CINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL CGAMX2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL CCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL CCHKAMX(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL CRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL CBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL CBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('COMPLEX AMX TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('COMPLEX AMX TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('COMPLEX AMX TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of CTESTAMX.
*
      END
*
      SUBROUTINE CRCCHK( IPRE, IPOST, PADVAL, M, N, RA, CA, LDI, MYROW,
     $                   MYCOL, TESTNUM, MAXERR, NERR,
     $                   ERRIBUF, ERRDBUF )
*
*     .. Scalar Arguments ..
      INTEGER IPRE, IPOST, PADVAL, M, N, LDI, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR)
      COMPLEX ERRDBUF(2, MAXERR)
*     ..
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, K, IAM
*     ..
*     .. Executable Statements ..
*
      IAM = MYROW * IBTNPROCS() + MYCOL
*
*     Check pre padding
*
      IF( LDI .NE. -1 ) THEN
         IF( IPRE .GT. 0 ) THEN
            DO 10 I = 1, IPRE
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -ERR_PRE
                     ERRDBUF(1, NERR) = CMPLX( RA(I) )
                     ERRDBUF(2, NERR) = CMPLX( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -10 - ERR_PRE
                     ERRDBUF(1, NERR) = CMPLX( CA(I) )
                     ERRDBUF(2, NERR) = CMPLX( PADVAL )
                  END IF
               ENDIF
   10       CONTINUE
         END IF
*
*        Check post padding
*
         IF( IPOST .GT. 0 ) THEN
            K = IPRE + LDI*N
            DO 20 I = K+1, K+IPOST
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -ERR_POST
                     ERRDBUF(1, NERR) = CMPLX( RA(I) )
                     ERRDBUF(2, NERR) = CMPLX( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -10 - ERR_POST
                     ERRDBUF(1, NERR) = CMPLX( CA(I) )
                     ERRDBUF(2, NERR) = CMPLX( PADVAL )
                  END IF
               ENDIF
   20       CONTINUE
         END IF
*
*        Check all (LDI-M) gaps
*
         IF( LDI .GT. M ) THEN
            K = IPRE + M + 1
            DO 40 J = 1, N
               DO 30 I = M+1, LDI
                  K = IPRE + (J-1)*LDI + I
                  IF( RA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -ERR_GAP
                        ERRDBUF(1, NERR) = CMPLX( RA(K) )
                        ERRDBUF(2, NERR) = CMPLX( PADVAL )
                     END IF
                  END IF
                  IF( CA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -10 - ERR_GAP
                        ERRDBUF(1, NERR) = CMPLX( CA(K) )
                        ERRDBUF(2, NERR) = CMPLX( PADVAL )
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
*
*     if RA and CA don't exist, buffs better be untouched
*
      ELSE
         DO 50 I = 1, IPRE+IPOST
            IF( RA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -ERR_PRE
                  ERRDBUF(1, NERR) = CMPLX( RA(I) )
                  ERRDBUF(2, NERR) = CMPLX( PADVAL )
               END IF
            END IF
            IF( CA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -10 - ERR_PRE
                  ERRDBUF(1, NERR) = CMPLX( CA(I) )
                  ERRDBUF(2, NERR) = CMPLX( PADVAL )
               END IF
            END IF
   50    CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE CCHKAMX( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      COMPLEX A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      REAL SBTEPS, CBTABS
      COMPLEX CBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, CBTRAN, SBTEPS, CBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMX, CAMX
      INTEGER IAMX, I, J, K, H, DEST, NODE
      REAL EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = CBTRAN( ISEED )
            IAMX = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  CBTRAN( ISEED(K*4+1) )
                  IF( CBTABS( VALS(K+1) ) .GT. CBTABS( VALS(IAMX) ) )
     $               IAMX = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMX) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = ABS( CBTABS(VALS(K)) - CBTABS(VALS(IAMX)) )
     $                       .GT. 3*EPS
                     IF( .NOT.ERROR ) IAMX = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ABS( CBTABS(A(I,J)) - CBTABS(VALS(IAMX)) )
     $                    .GT. 3*EPS
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMX)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMX ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMX) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMX-1, MYROW, MYCOL,
     $                                NPCOL, RAMX, CAMX )
                     IF( RAMX .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMX
                     END IF
                     IF( CAMX .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMX
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of CCHKAMX
*
      END
*
*
      SUBROUTINE ZAMXTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      DOUBLE COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ZTESTAMX:  Test double complex AMX COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) DOUBLE COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, ZGAMX2D
      EXTERNAL ZINITMAT, ZCHKPAD, ZBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR, ZSIZE
      DOUBLE COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = DCMPLX( -9.11D0, -9.21D0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      ZSIZE = IBTSIZEOF('Z')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( ZSIZE * (MEMLEN-I) ) / ( ZSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MAX tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL ZINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL ZGAMX2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ZCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ZCHKAMX(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL ZRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL ZBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL ZBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE COMPLEX AMX TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('DOUBLE COMPLEX AMX TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE COMPLEX AMX TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ZTESTAMX.
*
      END
*
      SUBROUTINE ZRCCHK( IPRE, IPOST, PADVAL, M, N, RA, CA, LDI, MYROW,
     $                   MYCOL, TESTNUM, MAXERR, NERR,
     $                   ERRIBUF, ERRDBUF )
*
*     .. Scalar Arguments ..
      INTEGER IPRE, IPOST, PADVAL, M, N, LDI, MYROW, MYCOL, TESTNUM
      INTEGER MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR)
      DOUBLE COMPLEX ERRDBUF(2, MAXERR)
*     ..
*     .. Parameters ..
      INTEGER ERR_PRE, ERR_POST, ERR_GAP, ERR_TRI, ERR_MAT
      PARAMETER( ERR_PRE = 1, ERR_POST = 2, ERR_GAP = 3, ERR_TRI = 4 )
      PARAMETER( ERR_MAT = 5 )
*     ..
*     .. External Functions ..
      INTEGER  IBTNPROCS
      EXTERNAL IBTNPROCS
*     ..
*     .. Local Scalars ..
      INTEGER I, J, K, IAM
*     ..
*     .. Executable Statements ..
*
      IAM = MYROW * IBTNPROCS() + MYCOL
*
*     Check pre padding
*
      IF( LDI .NE. -1 ) THEN
         IF( IPRE .GT. 0 ) THEN
            DO 10 I = 1, IPRE
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -ERR_PRE
                     ERRDBUF(1, NERR) = DCMPLX( RA(I) )
                     ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I
                     ERRIBUF(5, NERR) = IPRE - I + 1
                     ERRIBUF(6, NERR) = -10 - ERR_PRE
                     ERRDBUF(1, NERR) = DCMPLX( CA(I) )
                     ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                  END IF
               ENDIF
   10       CONTINUE
         END IF
*
*        Check post padding
*
         IF( IPOST .GT. 0 ) THEN
            K = IPRE + LDI*N
            DO 20 I = K+1, K+IPOST
               IF( RA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -ERR_POST
                     ERRDBUF(1, NERR) = DCMPLX( RA(I) )
                     ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                  END IF
               ENDIF
               IF( CA(I) .NE. PADVAL ) THEN
                  NERR = NERR + 1
                  IF( NERR .LE. MAXERR ) THEN
                     ERRIBUF(1, NERR) = TESTNUM
                     ERRIBUF(2, NERR) = LDI
                     ERRIBUF(3, NERR) = IAM
                     ERRIBUF(4, NERR) = I - K
                     ERRIBUF(5, NERR) = I
                     ERRIBUF(6, NERR) = -10 - ERR_POST
                     ERRDBUF(1, NERR) = DCMPLX( CA(I) )
                     ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                  END IF
               ENDIF
   20       CONTINUE
         END IF
*
*        Check all (LDI-M) gaps
*
         IF( LDI .GT. M ) THEN
            K = IPRE + M + 1
            DO 40 J = 1, N
               DO 30 I = M+1, LDI
                  K = IPRE + (J-1)*LDI + I
                  IF( RA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -ERR_GAP
                        ERRDBUF(1, NERR) = DCMPLX( RA(K) )
                        ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                     END IF
                  END IF
                  IF( CA(K) .NE. PADVAL) THEN
                     NERR = NERR + 1
                     IF( NERR .LE. MAXERR ) THEN
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = LDI
                        ERRIBUF(3, NERR) = IAM
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -10 - ERR_GAP
                        ERRDBUF(1, NERR) = DCMPLX( CA(K) )
                        ERRDBUF(2, NERR) = DCMPLX( PADVAL )
                     END IF
                  END IF
   30          CONTINUE
   40       CONTINUE
         END IF
*
*     if RA and CA don't exist, buffs better be untouched
*
      ELSE
         DO 50 I = 1, IPRE+IPOST
            IF( RA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -ERR_PRE
                  ERRDBUF(1, NERR) = DCMPLX( RA(I) )
                  ERRDBUF(2, NERR) = DCMPLX( PADVAL )
               END IF
            END IF
            IF( CA(I) .NE. PADVAL) THEN
               NERR = NERR + 1
               IF( NERR .LE. MAXERR ) THEN
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = LDI
                  ERRIBUF(3, NERR) = IAM
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = IPRE+IPOST
                  ERRIBUF(6, NERR) = -10 - ERR_PRE
                  ERRDBUF(1, NERR) = DCMPLX( CA(I) )
                  ERRDBUF(2, NERR) = DCMPLX( PADVAL )
               END IF
            END IF
   50    CONTINUE
      ENDIF
*
      RETURN
      END
*
      SUBROUTINE ZCHKAMX( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE COMPLEX A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      DOUBLE PRECISION DBTEPS, ZBTABS
      DOUBLE COMPLEX ZBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, ZBTRAN, DBTEPS, ZBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMX, CAMX
      INTEGER IAMX, I, J, K, H, DEST, NODE
      DOUBLE PRECISION EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = ZBTRAN( ISEED )
            IAMX = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  ZBTRAN( ISEED(K*4+1) )
                  IF( ZBTABS( VALS(K+1) ) .GT. ZBTABS( VALS(IAMX) ) )
     $               IAMX = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMX) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = ABS( ZBTABS(VALS(K)) - ZBTABS(VALS(IAMX)) )
     $                       .GT. 3*EPS
                     IF( .NOT.ERROR ) IAMX = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ABS( ZBTABS(A(I,J)) - ZBTABS(VALS(IAMX)) )
     $                    .GT. 3*EPS
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMX)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMX ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMX) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMX-1, MYROW, MYCOL,
     $                                NPCOL, RAMX, CAMX )
                     IF( RAMX .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMX
                     END IF
                     IF( CAMX .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMX
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ZCHKAMX
*
      END
*
*
      SUBROUTINE IAMNTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      INTEGER MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ITESTAMN:  Test integer AMN COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) INTEGER array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, IGAMN2D
      EXTERNAL IINITMAT, ICHKPAD, IBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      INTEGER CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -911
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( ISIZE * (MEMLEN-I) ) / ( ISIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MIN tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL IINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL IGAMN2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ICHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ICHKAMN(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL IRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL IBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL IBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('INTEGER AMN TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('INTEGER AMN TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('INTEGER AMN TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ITESTAMN.
*
      END
*
      SUBROUTINE ICHKAMN( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      INTEGER A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSPNUM, IBTRAN, IBTABS
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, IBTRAN
      EXTERNAL IBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMN, CAMN
      INTEGER IAMN, I, J, K, H, DEST, NODE
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = IBTRAN( ISEED )
            IAMN = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  IBTRAN( ISEED(K*4+1) )
                  IF( IBTABS( VALS(K+1) ) .LT. IBTABS( VALS(IAMN) ) )
     $               IAMN = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMN) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = IBTABS( VALS(K) ).NE.IBTABS( VALS(IAMN) )
                     IF( .NOT.ERROR ) IAMN = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( IBTABS( A(I,J) ) .NE. IBTABS( VALS(IAMN) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMN)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMN ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMN) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMN-1, MYROW, MYCOL,
     $                                NPCOL, RAMN, CAMN )
                     IF( RAMN .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMN
                     END IF
                     IF( CAMN .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMN
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ICHKAMN
*
      END
*
*
      SUBROUTINE SAMNTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      REAL MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  STESTAMN:  Test real AMN COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) REAL array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, SGAMN2D
      EXTERNAL SINITMAT, SCHKPAD, SBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, SSIZE, TESTNUM, VALPTR
      REAL CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.61E0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      SSIZE = IBTSIZEOF('S')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( SSIZE * (MEMLEN-I) ) / ( SSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MIN tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL SINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL SGAMN2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL SCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL SCHKAMN(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL SRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL SBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL SBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('REAL AMN TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('REAL AMN TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('REAL AMN TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of STESTAMN.
*
      END
*
      SUBROUTINE SCHKAMN( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      REAL A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      REAL SBTEPS, SBTABS
      REAL SBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, SBTRAN, SBTEPS, SBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMN, CAMN
      INTEGER IAMN, I, J, K, H, DEST, NODE
      REAL EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = SBTRAN( ISEED )
            IAMN = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  SBTRAN( ISEED(K*4+1) )
                  IF( SBTABS( VALS(K+1) ) .LT. SBTABS( VALS(IAMN) ) )
     $               IAMN = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMN) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = SBTABS( VALS(K) ).NE.SBTABS( VALS(IAMN) )
                     IF( .NOT.ERROR ) IAMN = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( SBTABS( A(I,J) ) .NE. SBTABS( VALS(IAMN) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMN)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMN ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMN) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMN-1, MYROW, MYCOL,
     $                                NPCOL, RAMN, CAMN )
                     IF( RAMN .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMN
                     END IF
                     IF( CAMN .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMN
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of SCHKAMN
*
      END
*
*
      SUBROUTINE DAMNTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      DOUBLE PRECISION MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  DTESTAMN:  Test double precision AMN COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) DOUBLE PRECISION array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, DGAMN2D
      EXTERNAL DINITMAT, DCHKPAD, DBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, DSIZE, ERRDPTR,
     $        ERRIPTR, I, IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST,
     $        IPRE, ISC, ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO,
     $        ITR, ITR1, ITR2, J, K, LDA, LDADST, LDASRC, LDI, M,
     $        MAXERR, MYCOL, MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP,
     $        PREAPTR, RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      DOUBLE PRECISION CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = -0.81D0
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      DSIZE = IBTSIZEOF('D')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( DSIZE * (MEMLEN-I) ) / ( DSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MIN tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL DINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL DGAMN2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL DCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL DCHKAMN(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL DRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL DBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL DBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE PRECISION AMN TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('DOUBLE PRECISION AMN TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE PRECISION AMN TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of DTESTAMN.
*
      END
*
      SUBROUTINE DCHKAMN( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE PRECISION A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      DOUBLE PRECISION DBTEPS, DBTABS
      DOUBLE PRECISION DBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, DBTRAN, DBTEPS, DBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMN, CAMN
      INTEGER IAMN, I, J, K, H, DEST, NODE
      DOUBLE PRECISION EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = DBTRAN( ISEED )
            IAMN = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  DBTRAN( ISEED(K*4+1) )
                  IF( DBTABS( VALS(K+1) ) .LT. DBTABS( VALS(IAMN) ) )
     $               IAMN = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMN) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = DBTABS( VALS(K) ).NE.DBTABS( VALS(IAMN) )
                     IF( .NOT.ERROR ) IAMN = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ( DBTABS( A(I,J) ) .NE. DBTABS( VALS(IAMN) ) )
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMN)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMN ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMN) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMN-1, MYROW, MYCOL,
     $                                NPCOL, RAMN, CAMN )
                     IF( RAMN .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMN
                     END IF
                     IF( CAMN .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMN
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of DCHKAMN
*
      END
*
*
      SUBROUTINE CAMNTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  CTESTAMN:  Test complex AMN COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, CGAMN2D
      EXTERNAL CINITMAT, CCHKPAD, CBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, CSIZE, ERRDPTR,
     $        ERRIPTR, I, IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST,
     $        IPRE, ISC, ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO,
     $        ITR, ITR1, ITR2, J, K, LDA, LDADST, LDASRC, LDI, M,
     $        MAXERR, MYCOL, MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP,
     $        PREAPTR, RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR
      COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = CMPLX( -0.91E0, -0.71E0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      CSIZE = IBTSIZEOF('C')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( CSIZE * (MEMLEN-I) ) / ( CSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MIN tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL CINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL CGAMN2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL CCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL CCHKAMN(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL CRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL CBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL CBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('COMPLEX AMN TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('COMPLEX AMN TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('COMPLEX AMN TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of CTESTAMN.
*
      END
*
      SUBROUTINE CCHKAMN( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      COMPLEX A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      REAL SBTEPS, CBTABS
      COMPLEX CBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, CBTRAN, SBTEPS, CBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMN, CAMN
      INTEGER IAMN, I, J, K, H, DEST, NODE
      REAL EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = SBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = CBTRAN( ISEED )
            IAMN = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  CBTRAN( ISEED(K*4+1) )
                  IF( CBTABS( VALS(K+1) ) .LT. CBTABS( VALS(IAMN) ) )
     $               IAMN = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMN) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = ABS( CBTABS(VALS(K)) - CBTABS(VALS(IAMN)) )
     $                       .GT. 3*EPS
                     IF( .NOT.ERROR ) IAMN = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ABS( CBTABS(A(I,J)) - CBTABS(VALS(IAMN)) )
     $                    .GT. 3*EPS
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMN)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMN ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMN) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMN-1, MYROW, MYCOL,
     $                                NPCOL, RAMN, CAMN )
                     IF( RAMN .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMN
                     END IF
                     IF( CAMN .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMN
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of CCHKAMN
*
      END
*
*
      SUBROUTINE ZAMNTEST( OUTNUM, VERB, TOPSREPEAT, TOPSCOHRNT, NSCOPE,
     $                     SCOPE0, NTOP, TOP0, NMAT, M0, N0, LDAS0,
     $                     LDAD0, LDI0, NDEST, RDEST0, CDEST0, NGRID,
     $                     CONTEXT0, P0, Q0, ISEED, RMEM, CMEM, RCLEN,
     $                     MEM, MEMLEN )
*
*  -- BLACS tester (version 1.0) --
*  University of Tennessee
*  December 15, 1994
*
*
*     .. Scalar Arguments ..
      INTEGER MEMLEN, NDEST, NGRID, NMAT, NSCOPE, NTOP, OUTNUM, RCLEN,
     $        TOPSCOHRNT, TOPSREPEAT, VERB
*     ..
*     .. Array Arguments ..
      CHARACTER*1 SCOPE0(NSCOPE), TOP0(NTOP)
      INTEGER M0(NMAT), N0(NMAT), LDAS0(NMAT), LDAD0(NMAT), LDI0(NMAT)
      INTEGER RDEST0(NDEST), CDEST0(NDEST), CONTEXT0(NGRID)
      INTEGER P0(NGRID), Q0(NGRID), ISEED(*), RMEM(RCLEN), CMEM(RCLEN)
      DOUBLE COMPLEX MEM(MEMLEN)
*     ..
*
*  Purpose
*  =======
*  ZTESTAMN:  Test double complex AMN COMBINE
*
*  Arguments
*  =========
*  OUTNUM   (input) INTEGER
*           The device number to write output to.
*
*  VERB     (input) INTEGER
*           The level of verbosity (how much printing to do).
*
*  NSCOPE   (input) INTEGER
*           The number of scopes to be tested.
*
*  SCOPE0   (input) CHARACTER*1 array of dimension (NSCOPE)
*           Values of the scopes to be tested.
*
*  NTOP     (input) INTEGER
*           The number of topologies to be tested.
*
*  TOP0     (input) CHARACTER*1 array of dimension (NTOP)
*           Values of the topologies to be tested.
*
*  NMAT     (input) INTEGER
*           The number of matrices to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  M0       (input) INTEGER array of dimension (NMAT)
*           Values of M to be tested.
*
*  N0       (input) INTEGER array of dimension (NMAT)
*           Values of N to be tested.
*
*  LDAS0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAS (leading dimension of A on source process)
*           to be tested.
*
*  LDAD0    (input) INTEGER array of dimension (NMAT)
*           Values of LDAD (leading dimension of A on destination
*           process) to be tested.
*  LDI0     (input) INTEGER array of dimension (NMAT)
*           Values of LDI (leading dimension of RA/CA) to be tested.
*           If LDI == -1, these RA/CA should not be accessed.
*
*  NDEST    (input) INTEGER
*           The number of destinations to be tested.
*
*  RDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of RDEST (row coordinate of destination) to be
*           tested.
*
*  CDEST0   (input) INTEGER array of dimension (NNDEST)
*           Values of CDEST (column coordinate of destination) to be
*           tested.
*
*  NGRID    (input) INTEGER
*           The number of process grids to be tested.
*
*  CONTEXT0 (input) INTEGER array of dimension (NGRID)
*           The BLACS context handles corresponding to the grids.
*
*  P0       (input) INTEGER array of dimension (NGRID)
*           Values of P (number of process rows, NPROW).
*
*  Q0       (input) INTEGER array of dimension (NGRID)
*           Values of Q (number of process columns, NPCOL).
*
*  ISEED    (workspace) INTEGER array of dimension ( MAX(NPROCS, NTESTS) )
*           Workspace used to hold each process's random number SEED.
*           This requires NPROCS (number of processor) elements.
*           If VERB < 2, this workspace also serves to indicate which
*           tests fail.  This requires workspace of NTESTS
*           (number of tests performed).
*
*  RMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all RA arrays, and their pre and post padding.
*
*  CMEM     (workspace) INTEGER array of dimension (RCLEN)
*           Used for all CA arrays, and their pre and post padding.
*
*  RCLEN    (input) INTEGER
*           The length, in elements, of RMEM and CMEM.
*
*  MEM      (workspace) DOUBLE COMPLEX array of dimension (MEMLEN)
*           Used for all other workspaces, including the matrix A,
*           and its pre and post padding.
*
*  MEMLEN   (input) INTEGER
*           The length, in elements, of MEM.
*
* =====================================================================
*
*     .. External Functions ..
      LOGICAL  ALLPASS, LSAME
      INTEGER  IBTMYPROC, IBTNPROCS, IBTSIZEOF
      EXTERNAL ALLPASS, LSAME, IBTMYPROC, IBTNPROCS, IBTSIZEOF
*     ..
*     .. External Subroutines ..
      EXTERNAL BLACS_GRIDINFO, ZGAMN2D
      EXTERNAL ZINITMAT, ZCHKPAD, ZBTCHECKIN
*     ..
*     .. Local Scalars ..
      CHARACTER*1 SCOPE, TOP
      LOGICAL INGRID, TESTOK, ALLRCV
      INTEGER APTR, CAPTR, CDEST, CDEST2, CONTEXT, ERRDPTR, ERRIPTR, I,
     $        IAM, ICHECKVAL, IDE, IGR, IMA, IPAD, IPOST, IPRE, ISC,
     $        ISIZE, ISTART, ISTOP, ITC, ITC1, ITC2, ITO, ITR, ITR1,
     $        ITR2, J, K, LDA, LDADST, LDASRC, LDI, M, MAXERR, MYCOL,
     $        MYROW, N, NERR, NFAIL, NPCOL, NPROW, NSKIP, PREAPTR,
     $        RAPTR, RDEST, RDEST2, SETWHAT, TESTNUM, VALPTR, ZSIZE
      DOUBLE COMPLEX CHECKVAL
*     ..
*     .. Executable Statements ..
*
*     Choose padding value, and make it unique
*
      CHECKVAL = DCMPLX( -9.11D0, -9.21D0 )
      IAM = IBTMYPROC()
      CHECKVAL = IAM * CHECKVAL
      ISIZE = IBTSIZEOF('I')
      ZSIZE = IBTSIZEOF('Z')
      ICHECKVAL = -IAM
*
*     Verify file parameters
*
      IF( IAM .EQ. 0 ) THEN
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, *) '  '
         WRITE(OUTNUM, 1000 )
         IF( VERB .GT. 0 ) THEN
            WRITE(OUTNUM,*) '  '
            WRITE(OUTNUM, 2000) 'NSCOPE:', NSCOPE
            WRITE(OUTNUM, 3000) ' SCOPE:', ( SCOPE0(I), I = 1, NSCOPE )
            WRITE(OUTNUM, 2000) 'TReps :', TOPSREPEAT
            WRITE(OUTNUM, 2000) 'TCohr :', TOPSCOHRNT
            WRITE(OUTNUM, 2000) 'NTOP  :', NTOP
            WRITE(OUTNUM, 3000) ' TOP  :', ( TOP0(I), I = 1, NTOP )
            WRITE(OUTNUM, 2000) 'NMAT  :', NMAT
            WRITE(OUTNUM, 2000) ' M    :', ( M0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' N    :', ( N0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAS :', ( LDAS0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDAD :', ( LDAD0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) ' LDI  :', ( LDI0(I), I = 1, NMAT )
            WRITE(OUTNUM, 2000) 'NDEST :', NDEST
            WRITE(OUTNUM, 2000) ' RDEST:',( RDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) ' CDEST:',( CDEST0(I), I = 1, NDEST )
            WRITE(OUTNUM, 2000) 'NGRIDS:', NGRID
            WRITE(OUTNUM, 2000) ' P    :', ( P0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) ' Q    :', ( Q0(I), I = 1, NGRID )
            WRITE(OUTNUM, 2000) 'VERB  :', VERB
            WRITE(OUTNUM,*) '  '
         END IF
         IF( VERB .GT. 1 ) THEN
            WRITE(OUTNUM,4000)
            WRITE(OUTNUM,5000)
         END IF
      END IF
      IF (TOPSREPEAT.EQ.0) THEN
         ITR1 = 0
         ITR2 = 0
      ELSE IF (TOPSREPEAT.EQ.1) THEN
         ITR1 = 1
         ITR2 = 1
      ELSE
         ITR1 = 0
         ITR2 = 1
      END IF
*
*     Find biggest matrix, so we know where to stick error info
*
      I = 0
      DO 10 IMA = 1, NMAT
         IPAD = 4 * M0(IMA)
         K = N0(IMA) * MAX0( LDAS0(IMA), LDAD0(IMA) ) + IPAD
         IF( K .GT. I ) I = K
   10  CONTINUE
      I = I + IBTNPROCS()
      MAXERR = ( ZSIZE * (MEMLEN-I) ) / ( ZSIZE*2 + ISIZE*6 )
      IF( MAXERR .LT. 1 ) THEN
         WRITE(OUTNUM,*) 'ERROR: Not enough memory to run MIN tests.'
         CALL BLACS_ABORT(-1, 1)
      END IF
      ERRDPTR = I + 1
      ERRIPTR = ERRDPTR + MAXERR
      NERR = 0
      TESTNUM = 0
      NFAIL = 0
      NSKIP = 0
*
*     Loop over grids of matrix
*
      DO 90 IGR = 1, NGRID
*
*        allocate process grid for the next batch of tests
*
         CONTEXT = CONTEXT0(IGR)
         CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
         INGRID = ( (MYROW.LT.NPROW) .AND. (MYCOL.LT.NPCOL) )
*
         DO 80 ISC = 1, NSCOPE
            SCOPE = SCOPE0(ISC)
            DO 70 ITO = 1, NTOP
               TOP = TOP0(ITO)
*
*              If testing multiring ('M') or general tree ('T'), need to
*              loop over calls to BLACS_SET to do full test
*
               IF( LSAME(TOP, 'M') ) THEN
                  SETWHAT = 13
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTART = -(NPCOL - 1)
                     ISTOP = -ISTART
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTART = -(NPROW - 1)
                     ISTOP = -ISTART
                  ELSE
                     ISTART = -(NPROW*NPCOL - 1)
                     ISTOP = -ISTART
                  ENDIF
               ELSE IF( LSAME(TOP, 'T') ) THEN
                  SETWHAT = 14
                  ISTART = 1
                  IF( SCOPE .EQ. 'R' ) THEN
                     ISTOP = NPCOL - 1
                  ELSE IF (SCOPE .EQ. 'C') THEN
                     ISTOP = NPROW - 1
                  ELSE
                     ISTOP = NPROW*NPCOL - 1
                  ENDIF
               ELSE
                  SETWHAT = 0
                  ISTART = 1
                  ISTOP = 1
               ENDIF
               DO 60 IMA = 1, NMAT
                  M = M0(IMA)
                  N = N0(IMA)
                  LDASRC = LDAS0(IMA)
                  LDADST = LDAD0(IMA)
                  LDI = LDI0(IMA)
                  IPRE  = 2 * M
                  IPOST = IPRE
                  PREAPTR = 1
                  APTR = PREAPTR + IPRE
*
                  DO 50 IDE = 1, NDEST
                     TESTNUM = TESTNUM + 1
                     RDEST2 = RDEST0(IDE)
                     CDEST2 = CDEST0(IDE)
*
*                    If everyone gets the answer, create some bogus rdest/cdest
*                    so IF's are easier
*
                     ALLRCV = ( (RDEST2.EQ.-1) .OR. (CDEST2.EQ.-1) )
                     IF( ALLRCV ) THEN
                        RDEST = NPROW - 1
                        CDEST = NPCOL - 1
                        IF (TOPSCOHRNT.EQ.0) THEN
                           ITR1 = 0
                           ITR2 = 0
                        ELSE IF (TOPSCOHRNT.EQ.1) THEN
                           ITR1 = 1
                           ITR2 = 1
                        ELSE
                           ITR1 = 0
                           ITR2 = 1
                        END IF
                     ELSE
                        RDEST = RDEST2
                        CDEST = CDEST2
                        ITC1 = 0
                        ITC2 = 0
                     END IF
                     IF( RDEST.GE.P0(IGR) .OR. CDEST.GE.Q0(IGR) ) THEN
                        NSKIP = NSKIP + 1
                        GOTO 50
                     END IF
*
                     IF( MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST ) THEN
                        LDA = LDADST
                     ELSE
                        LDA = LDASRC
                     END IF
                     VALPTR = APTR + IPOST + N * LDA
                     IF( VERB .GT. 1 ) THEN
                        IF( IAM .EQ. 0 ) THEN
                           WRITE(OUTNUM, 6000)
     $                     TESTNUM, 'RUNNING', SCOPE, TOP, M, N,
     $                     LDASRC, LDADST, LDI, RDEST2, CDEST2,
     $                     NPROW, NPCOL
                        END IF
                     END IF
*
*                    If I am in scope
*
                     TESTOK = .TRUE.
                     IF( INGRID ) THEN
                        IF( (MYROW.EQ.RDEST .AND. SCOPE.EQ.'R') .OR.
     $                      (MYCOL.EQ.CDEST .AND. SCOPE.EQ.'C') .OR.
     $                      (SCOPE .EQ. 'A') ) THEN
*
                           K = NERR
                           DO 40 ITR = ITR1, ITR2
                              CALL BLACS_SET(CONTEXT, 15, ITR)
                           DO 35 ITC = ITC1, ITC2
                              CALL BLACS_SET(CONTEXT, 16, ITC)
                           DO 30 J = ISTART, ISTOP
                              IF( J.EQ.0) GOTO 30
                              IF( SETWHAT.NE.0 )
     $                           CALL BLACS_SET(CONTEXT, SETWHAT, J)
*
*
*                             generate and pad matrix A
*
                              CALL ZINITMAT('G','-', M, N, MEM(PREAPTR),
     $                                      LDA, IPRE, IPOST,
     $                                      CHECKVAL, TESTNUM,
     $                                      MYROW, MYCOL )
*
*                             If they exist, pad RA and CA arrays
*
                              IF( LDI .NE. -1 ) THEN
                                 DO 15 I = 1, N*LDI + IPRE + IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   15                            CONTINUE
                                 RAPTR = 1 + IPRE
                                 CAPTR = 1 + IPRE
                              ELSE
                                 DO 20 I = 1, IPRE+IPOST
                                    RMEM(I) = ICHECKVAL
                                    CMEM(I) = ICHECKVAL
   20                            CONTINUE
                                 RAPTR = 1
                                 CAPTR = 1
                              END IF
*
                              CALL ZGAMN2D(CONTEXT, SCOPE, TOP, M, N,
     $                                     MEM(APTR), LDA, RMEM(RAPTR),
     $                                     CMEM(CAPTR), LDI,
     $                                     RDEST2, CDEST2)
*
*                             If I've got the answer, check for errors in
*                             matrix or padding
*
                              IF( (MYROW.EQ.RDEST .AND. MYCOL.EQ.CDEST)
     $                            .OR. ALLRCV ) THEN
                                 CALL ZCHKPAD('G','-', M, N,
     $                                        MEM(PREAPTR), LDA, RDEST,
     $                                        CDEST, MYROW, MYCOL,
     $                                        IPRE, IPOST, CHECKVAL,
     $                                        TESTNUM, MAXERR, NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR))
                                 CALL ZCHKAMN(SCOPE, CONTEXT, M, N,
     $                                        MEM(APTR), LDA,
     $                                        RMEM(RAPTR), CMEM(CAPTR),
     $                                        LDI, TESTNUM, MAXERR,NERR,
     $                                        MEM(ERRIPTR),MEM(ERRDPTR),
     $                                        ISEED, MEM(VALPTR))
                                 CALL ZRCCHK(IPRE, IPOST, ICHECKVAL,
     $                                       M, N, RMEM, CMEM, LDI,
     $                                       MYROW, MYCOL, TESTNUM,
     $                                       MAXERR, NERR,
     $                                       MEM(ERRIPTR), MEM(ERRDPTR))
                              END IF
   30                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 16, 0)
   35                      CONTINUE
                           CALL BLACS_SET(CONTEXT, 15, 0)
   40                      CONTINUE
                        TESTOK = ( K .EQ. NERR )
                        END IF
                     END IF
*
                     IF( VERB .GT. 1 ) THEN
                        I = NERR
                        CALL ZBTCHECKIN(0, OUTNUM, MAXERR, NERR,
     $                               MEM(ERRIPTR), MEM(ERRDPTR), ISEED)
                        IF( IAM .EQ. 0 ) THEN
                           IF( TESTOK .AND. NERR.EQ.I ) THEN
                              WRITE(OUTNUM,6000)TESTNUM,'PASSED ',
     $                              SCOPE, TOP, M, N, LDASRC,
     $                              LDADST, LDI, RDEST2, CDEST2,
     $                              NPROW, NPCOL
                           ELSE
                              NFAIL = NFAIL + 1
                              WRITE(OUTNUM,6000)TESTNUM,'FAILED ',
     $                             SCOPE, TOP, M, N, LDASRC,
     $                             LDADST, LDI, RDEST2, CDEST2,
     $                             NPROW, NPCOL
                           END IF
                        END IF
*
*                       Once we've printed out errors, can re-use buf space
*
                        NERR = 0
                     END IF
   50             CONTINUE
   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE
*
      IF( VERB .LT. 2 ) THEN
         NFAIL = TESTNUM
         CALL ZBTCHECKIN( NFAIL, OUTNUM, MAXERR, NERR, MEM(ERRIPTR),
     $                    MEM(ERRDPTR), ISEED )
      END IF
      IF( IAM .EQ. 0 ) THEN
         IF( VERB .GT. 1 ) WRITE(OUTNUM,*) '   '
         IF( NFAIL+NSKIP .EQ. 0 ) THEN
            WRITE(OUTNUM, 7000 ) TESTNUM
         ELSE
            WRITE(OUTNUM, 8000 ) TESTNUM, TESTNUM-NSKIP-NFAIL,
     $                           NSKIP, NFAIL
         END IF
      END IF
*
*     Log whether their were any failures
*
      TESTOK = ALLPASS( (NFAIL.EQ.0) )
*
 1000 FORMAT('DOUBLE COMPLEX AMN TESTS: BEGIN.' )
 2000 FORMAT(1X,A7,3X,10I6)
 3000 FORMAT(1X,A7,3X,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,5X,A1,
     $       5X,A1,5X,A1)
 4000 FORMAT(' TEST#  STATUS SCOPE TOP     M     N  LDAS  LDAD   LDI ',
     $       'RDEST CDEST    P    Q')
 5000 FORMAT(' ----- ------- ----- --- ----- ----- ----- ----- ----- ',
     $       '----- ----- ---- ----')
 6000 FORMAT(I6,1X,A7,5X,A1,3X,A1,7I6,2I5)
 7000 FORMAT('DOUBLE COMPLEX AMN TESTS: PASSED ALL',
     $       I5, ' TESTS.')
 8000 FORMAT('DOUBLE COMPLEX AMN TESTS:',I5,' TESTS;',I5,' PASSED,',
     $       I5,' SKIPPED,',I5,' FAILED.')
*
      RETURN
*
*     End of ZTESTAMN.
*
      END
*
      SUBROUTINE ZCHKAMN( SCOPE, ICTXT, M, N, A, LDA, RA, CA, LDI,
     $                    TESTNUM, MAXERR, NERR, ERRIBUF, ERRDBUF,
     $                    ISEED, VALS )
*
*     .. Scalar Arguments ..
      CHARACTER*1 SCOPE
      INTEGER ICTXT, M, N, LDA, LDI, TESTNUM, MAXERR, NERR
*     ..
*     .. Array Arguments ..
      INTEGER RA(*), CA(*), ERRIBUF(6, MAXERR), ISEED(*)
      DOUBLE COMPLEX A(LDA,*), ERRDBUF(2, MAXERR), VALS(*)
*     ..
*     .. External Functions ..
      INTEGER IBTMYPROC, IBTNPROCS, IBTSPNUM
      DOUBLE PRECISION DBTEPS, ZBTABS
      DOUBLE COMPLEX ZBTRAN
      EXTERNAL IBTMYPROC, IBTNPROCS, IBTSPNUM, ZBTRAN, DBTEPS, ZBTABS
*     ..
*     .. External Subroutines ..
      EXTERNAL IBTSPCOORD
*     ..
*     .. Local Scalars ..
      LOGICAL ERROR
      INTEGER NPROCS, NNODES, NPROW, NPCOL, MYROW, MYCOL, RAMN, CAMN
      INTEGER IAMN, I, J, K, H, DEST, NODE
      DOUBLE PRECISION EPS
*     ..
*     .. Executable Statements ..
*
      NPROCS = IBTNPROCS()
      EPS = DBTEPS()
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
      DEST = MYROW*NPROCS + MYCOL
*
*     Set up seeds to match those used by each proc's genmat call
*
      IF( SCOPE .EQ. 'R' ) THEN
         NNODES = NPCOL
         DO 10 I = 0, NNODES-1
            NODE = MYROW * NPROCS + I
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   10    CONTINUE
      ELSE IF( SCOPE .EQ. 'C' ) THEN
         NNODES = NPROW
         DO 20 I = 0, NNODES-1
            NODE = I * NPROCS + MYCOL
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   20    CONTINUE
      ELSE
         NNODES = NPROW * NPCOL
         DO 30 I = 0, NNODES-1
            NODE = (I / NPCOL) * NPROCS + MOD(I, NPCOL)
            ISEED(I*4+1) = MOD( 1002 + TESTNUM*5 + NODE*3, 4096 )
            ISEED(I*4+2) = MOD( 2027 + TESTNUM*7 + NODE, 4096 )
            ISEED(I*4+3) = MOD( 1234 + TESTNUM + NODE*3, 4096 )
            ISEED(I*4+4) = MOD( 4311 + TESTNUM*10 + NODE*2, 4096 )
   30    CONTINUE
      END IF
*
      DO 100 J = 1, N
         DO 90 I = 1, M
            H = (J-1)*LDI + I
            VALS(1) = ZBTRAN( ISEED )
            IAMN = 1
            IF( NNODES .GT. 1 ) THEN
               DO 40 K = 1, NNODES-1
                  VALS(K+1) =  ZBTRAN( ISEED(K*4+1) )
                  IF( ZBTABS( VALS(K+1) ) .LT. ZBTABS( VALS(IAMN) ) )
     $               IAMN = K + 1
   40          CONTINUE
            END IF
*
*           If BLACS have not returned same value we've chosen
*
            IF( A(I,J) .NE. VALS(IAMN) ) THEN
*
*              If we have RA and CA arrays
*
               IF( LDI .NE. -1 ) THEN
*
*                 Any number having the same absolute value is a valid max
*
                  K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
                  IF( K.GT.0 .AND. K.LE.NNODES ) THEN
                     ERROR = ABS( ZBTABS(VALS(K)) - ZBTABS(VALS(IAMN)) )
     $                       .GT. 3*EPS
                     IF( .NOT.ERROR ) IAMN = K
                  ELSE
                     ERROR = .TRUE.
                  END IF
               ELSE
*
*                 Error if BLACS answer not same absolute value, or if it
*                 was not really in the numbers being compared
*
                  ERROR = ABS( ZBTABS(A(I,J)) - ZBTABS(VALS(IAMN)) )
     $                    .GT. 3*EPS
                  IF( .NOT.ERROR ) THEN
                     DO 50 K = 1, NNODES
                        IF( VALS(K) .EQ. A(I,J) ) GOTO 60
   50                CONTINUE
                     ERROR = .TRUE.
   60                CONTINUE
                  ENDIF
               END IF
*
*              If the value is in error
*
               IF( ERROR ) THEN
                  NERR = NERR + 1
                  ERRIBUF(1, NERR) = TESTNUM
                  ERRIBUF(2, NERR) = NNODES
                  ERRIBUF(3, NERR) = DEST
                  ERRIBUF(4, NERR) = I
                  ERRIBUF(5, NERR) = J
                  ERRIBUF(6, NERR) = 5
                  ERRDBUF(1, NERR) = A(I,J)
                  ERRDBUF(2, NERR) = VALS(IAMN)
               END IF
            END IF
*
*           If they are defined, make sure coordinate entries are OK
*
            IF( LDI .NE. -1 ) THEN
               K = IBTSPNUM( SCOPE, RA(H), CA(H), NPCOL ) + 1
               IF( K.NE.IAMN ) THEN
*
*                 Make sure more than one proc doesn't have exact same value
*                 (and therefore there may be more than one valid coordinate
*                 for a single value)
*
                  IF( K.GT.NNODES .OR. K.LT.1 ) THEN
                     ERROR = .TRUE.
                  ELSE
                     ERROR = ( VALS(K) .NE. VALS(IAMN) )
                  END IF
                  IF( ERROR ) THEN
                     CALL IBTSPCOORD( SCOPE, IAMN-1, MYROW, MYCOL,
     $                                NPCOL, RAMN, CAMN )
                     IF( RAMN .NE. RA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -5
                        ERRDBUF(1, NERR) = RA(H)
                        ERRDBUF(2, NERR) = RAMN
                     END IF
                     IF( CAMN .NE. CA(H) ) THEN
                        NERR = NERR + 1
                        ERRIBUF(1, NERR) = TESTNUM
                        ERRIBUF(2, NERR) = NNODES
                        ERRIBUF(3, NERR) = DEST
                        ERRIBUF(4, NERR) = I
                        ERRIBUF(5, NERR) = J
                        ERRIBUF(6, NERR) = -15
                        ERRDBUF(1, NERR) = CA(H)
                        ERRDBUF(2, NERR) = CAMN
                     END IF
                  END IF
               END IF
            END IF
   90    CONTINUE
  100 CONTINUE
*
      RETURN
*
*     End of ZCHKAMN
*
      END
*
