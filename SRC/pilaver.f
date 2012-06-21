      SUBROUTINE PILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
C
C  -- ScaLAPACK computational routine (version 2.0.1 ) --
C  -- ScaLAPACK is a software package provided by Univ. of Tennessee,    --
C  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
C     January 2012
C
C  Purpose
C  =======
C
C  This subroutine return the ScaLAPACK version.
C
C  Arguments
C  =========
C  VERS_MAJOR   (output) INTEGER
C      return the scalapack major version
C  VERS_MINOR   (output) INTEGER
C      return the scalapack minor version from the major version
C  VERS_PATCH   (output) INTEGER
C      return the scalapack patch version from the minor version
C  =====================================================================
C
      INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
C  =====================================================================
      VERS_MAJOR = 2
      VERS_MINOR = 0
      VERS_PATCH = 2
C  =====================================================================
C
      RETURN
      END

