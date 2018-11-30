MODULE nrtype
  !Symbolic names for kind types of 4-, 2-, and 1-byte integers

  IMPLICIT NONE

  PUBLIC

  !Parameter-----------------------------------------------------------------
  INTEGER, PARAMETER :: I8B = SELECTED_INT_KIND(18) ! double precision integer
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)  ! single precision integer
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  !Symbolic names for kind types of single- and double-precision reals
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  !Symbolic names for kind types of 4-, 2-, and 1-byte complex
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  !Symbolic names for kind type of default logical
!  INTEGER, PARAMETER :: LGT = KIND(.true.)
  !Frequently used mathematical constants (with precision to spare)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp


  INTEGER, PARAMETER  :: StringLength         = 128
  INTEGER, PARAMETER  :: PathLength           = 256
  INTEGER,       PARAMETER :: null_int   = -9999999
  REAL,           PARAMETER :: null_real  = -9.9E15
  CHARACTER(LEN = 5), PARAMETER :: null_str   = '*****'

  INTEGER, PARAMETER :: verbose = 0
    !Errors
  INTEGER, parameter :: SUCCESS_               = 0

  ! Global variables
  CHARACTER(len=PathLength) :: dsnsim          !I parent folder for ensemble simulation
  CHARACTER(len=5)          :: pid             !I 5-digit pid to the ensemble member


!Derived data types for sparse matrices, single and double precision (see use in Chapter B2)
TYPE sprs2_sp
 INTEGER(I4B) :: n,len
 REAL(SP), DIMENSION(:), POINTER :: val
 INTEGER(I4B), DIMENSION(:), POINTER :: irow
 INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_sp
TYPE sprs2_dp
 INTEGER(I4B) :: n,len
 REAL(DP), DIMENSION(:), POINTER :: val
 INTEGER(I4B), DIMENSION(:), POINTER :: irow
 INTEGER(I4B), DIMENSION(:), POINTER :: jcol
END TYPE sprs2_dp
END MODULE nrtype
