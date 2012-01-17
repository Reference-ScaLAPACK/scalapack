      INTEGER FUNCTION PIPARMQ( ICTXT, ISPEC, NAME, OPTS, N, ILO, IHI,
     $                          LWORKNB )
*
*     Contribution from the Department of Computing Science and HPC2N,
*     Umea University, Sweden
*
*  -- ScaLAPACK auxiliary routine (version 2.0.1) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     Univ. of Colorado Denver and University of California, Berkeley.
*     January, 2012
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      INTEGER            ICTXT, IHI, ILO, ISPEC, LWORKNB, N
      CHARACTER          NAME*( * ), OPTS*( * )
*
*  Purpose
*  =======
*
*       This program sets problem and machine dependent parameters
*       useful for PxHSEQR and its subroutines. It is called whenever
*       PILAENVX is called with 12 <= ISPEC <= 16
*
*  Arguments
*  =========
*
*       ICTXT  (local input) INTEGER
*              On entry,  ICTXT  specifies the BLACS context handle,
*              indicating the global  context of the operation. The
*              context itself is global, but the value of ICTXT is
*              local.
*
*       ISPEC  (global input) INTEGER
*              ISPEC specifies which tunable parameter PIPARMQ should
*              return.
*
*              ISPEC=12: (INMIN)  Matrices of order nmin or less
*                        are sent directly to PxLAHQR, the implicit
*                        double shift QR algorithm.  NMIN must be
*                        at least 11.
*
*              ISPEC=13: (INWIN)  Size of the deflation window.
*                        This is best set greater than or equal to
*                        the number of simultaneous shifts NS.
*                        Larger matrices benefit from larger deflation
*                        windows.
*
*              ISPEC=14: (INIBL) Determines when to stop nibbling and
*                        invest in an (expensive) multi-shift QR sweep.
*                        If the aggressive early deflation subroutine
*                        finds LD converged eigenvalues from an order
*                        NW deflation window and LD.GT.(NW*NIBBLE)/100,
*                        then the next QR sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  Setting
*                        PIPARMQ(ISPEC=14) = 0 causes PxLAQR0 to skip a
*                        multi-shift QR sweep whenever early deflation
*                        finds a converged eigenvalue.  Setting
*                        PIPARMQ(ISPEC=14) greater than or equal to 100
*                        prevents PxLAQR0 from skipping a multi-shift
*                        QR sweep.
*
*              ISPEC=15: (NSHFTS) The number of simultaneous shifts in
*                        a multi-shift QR iteration.
*
*              ISPEC=16: (IACC22) PIPARMQ is set to 1 or 2 with the
*                        following meanings.
*                        1:  During the multi-shift QR sweep,
*                            PxLAQR5 and/or xLAQR6 accumulates reflections
*                            and uses matrix-matrix multiply to update
*                            the far-from-diagonal matrix entries.
*                        2:  During the multi-shift QR sweep.
*                            PxLAQR5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*
*                        ( IACC22=0 is valid in LAPACK but not here.
*                        Householder reflections are always accumulated
*                        for the performance consideration.
*                        If xTRMM is slower than xGEMM or NB is small,
*                        PIPARMQ(ISPEC=16)=1 may be more efficient than
*                        PIPARMQ(ISPEC=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice. )
*
*       NAME    (global input) character string
*               Name of the calling subroutine
*
*       OPTS    (global input) character string
*               This is a concatenation of the string arguments to
*               TTQRE.
*
*       N       (global input) integer scalar
*               N is the order of the Hessenberg matrix H.
*
*       ILO     (global input) INTEGER
*       IHI     (global input) INTEGER
*               It is assumed that H is already upper triangular
*               in rows and columns 1:ILO-1 and IHI+1:N.
*
*       LWORKNB   (global input) INTEGER
*               The amount of workspace available or the blockfactor.
*
*  Further Details
*  ===============
*
*       Little is known about how best to choose these parameters.
*       It is possible to use different values of the parameters
*       for each of PCHSEQR, PDHSEQR, PSHSEQR and PZHSEQR.
*
*       It is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been fully
*       implemented --- yet.
*
*
*       The best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of PxLAQR3 and PxLAQR5 and on the nature of each
*       particular eigenvalue problem.  Experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       Following is a list of default values supplied by PIPARMQ.
*       These defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       PIPARMQ(ISPEC=12) The PxLAQR1 vs PxLAQR0 crossover point.
*                         Default: 220. (Must be at least 11.)
*
*       PIPARMQ(ISPEC=13) Recommended deflation window size.
*                         This depends on ILO, IHI and NS, the
*                         number of simultaneous shifts returned
*                         by PIPARMQ(ISPEC=15).  The default for
*                         (IHI-ILO+1).LE.500 is NS.  The default
*                         for (IHI-ILO+1).GT.500 is 3*NS/2.
*
*       PIPARMQ(ISPEC=14) Nibble crossover point.
*                         The default for the serial case is 14.
*                         The default for the parallel case is
*                         335 * N**(-0.44) * NPROCS.
*
*       PIPARMQ(ISPEC=15) Number of simultaneous shifts, NS.
*                         a multi-shift QR iteration.
*
*                         If IHI-ILO+1 is ...
*
*                         greater than      ...but less    ... the
*                         or equal to ...      than        default is
*
*                                 0               30       NS =    2+
*                                30               60       NS =    4+
*                                60              150       NS =   10
*                               150              590       NS =   **
*                               590             3000       NS =   64
*                              3000             6000       NS =  128
*                              6000            12000       NS =  256
*                             12000            24000       NS =  512
*                             24000            48000       NS = 1024
*                             48000            96000       NS = 2048
*                             96000         INFINITY       NS = 4096
*
*                     (+)  By default matrices of this order are
*                          passed to the implicit double shift routine
*                          PxLAQR1.  See PIPARMQ(ISPEC=12) above. These
*                          values of NS are used only in case of a rare
*                          PxLAQR1 failure.
*
*                     (**) The asterisks (**) indicate an ad-hoc
*                          function increasing from 10 to 64.
*
*       PIPARMQ(ISPEC=16) Select structured matrix multiply.
*                         (See ISPEC=16 above for details.)
*                         Default: 3.
*
*     ================================================================
*     .. Parameters ..
      INTEGER            INMIN, INWIN, INIBL, ISHFTS, IACC22
      PARAMETER          ( INMIN = 12, INWIN = 13, INIBL = 14,
     $                   ISHFTS = 15, IACC22 = 16 )
      INTEGER            NMIN, NMIN2, K22MIN, KACMIN, NIBBLE, KNWSWP
      PARAMETER          ( NMIN = 220, K22MIN = 14, KACMIN = 14,
     $                   NIBBLE = 14, KNWSWP = 500, NMIN2 = 770 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0 )
*     ..
*     .. Local Scalars ..
      INTEGER            NH, NS, MYROW, MYCOL, NPROW, NPCOL, NP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LOG, MAX, MOD, NINT, REAL
*     ..
*     .. External functions ..
      INTEGER            ICEIL
      EXTERNAL           ICEIL
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO
*     ..
*     .. Executable Statements ..
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( ( ISPEC.EQ.ISHFTS ) .OR. ( ISPEC.EQ.INWIN ) .OR.
     $    ( ISPEC.EQ.IACC22 ) ) THEN
*
*        ==== Set the number simultaneous shifts ====
*
         NH = IHI - ILO + 1
         NS = 2
         IF( NH.GE.30 )
     $        NS = 4
         IF( NH.GE.60 )
     $        NS = 10
         IF( NH.GE.150 )
     $        NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ))
         IF( NH.GE.590 )
     $        NS = 64
         IF( NH.GE.3000 )
     $        NS = 128
         IF( NH.GE.6000 )
     $        NS = 256
         IF( NH.GE.12000 )
     $        NS = 512
         IF( NH.GE.24000 )
     $        NS = 1024
         IF( NH.GE.48000 )
     $        NS = 2048
         IF( NH.GE.96000 )
     $        NS = 4096
         IF( NH.GE.192000 )
     $        NS = 8192
         IF( NH.GE.384000 )
     $        NS = 16384
         IF( NH.GE.768000 )
     $        NS = 32768
         IF( NH.GE.1000000 )
     $        NS = ICEIL( NH, 25 )
         NS = MAX( NS, 2*MIN(NPROW,NPCOL) )
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      END IF
*
      IF( ISPEC.EQ.INMIN ) THEN
*
*
*        ===== Submatrices of order smaller than NMIN*min(P_r,P_c)
*        .     get sent to PxLAHQR, the classic ScaLAPACK algorithm.
*        .     This must be at least 11. ====
*
         PIPARMQ = NMIN * MIN( NPROW, NPCOL )
*
      ELSE IF( ISPEC.EQ.INIBL ) THEN
*
*        ==== INIBL: skip a multi-shift QR iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (NIBBLE*(window size)/100) deflations. ====
*
         NP = MIN( NPROW, NPCOL )
         IF( NP.EQ.1 ) THEN
            PIPARMQ = NIBBLE
         ELSE
            NH = IHI - ILO + 1
            PIPARMQ = MIN( 100,
     $           CEILING( 335.0D+0 * NH**(-0.44D+0) * NP ) )
         END IF
*
      ELSE IF( ISPEC.EQ.ISHFTS ) THEN
*
*        ==== NSHFTS: The number of simultaneous shifts =====
*
         PIPARMQ = NS
*
      ELSE IF( ISPEC.EQ.INWIN ) THEN
*
*        ==== NW: deflation window size.  ====
*
         IF( NH.LE.KNWSWP ) THEN
            PIPARMQ = NS
         ELSE
            PIPARMQ = 3*NS / 2
         END IF
*
      ELSE IF( ISPEC.EQ.IACC22 ) THEN
*
*        ==== IACC22: Whether to use 2-by-2 block structure while
*        .     doing it.  A small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     NH=IHI-ILO+1.
*
         PIPARMQ = 1
c         PIPARMQ = 0
c         IF( NS.GE.KACMIN )
c     $      PIPARMQ = 1
         IF( NS.GE.K22MIN )
     $      PIPARMQ = 2
*
      ELSE
*        ===== invalid value of ispec =====
         PIPARMQ = -1
*
      END IF
*
*     ==== End of PIPARMQ ====
*
      END
