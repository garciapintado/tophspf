PROGRAM test_srch

  ! Purpose:
  ! test the SUBROUTINE srch for mapping time series of inflows into a general grid
  ! by approximating the area of the inflow boundary condition and then
  ! interpolating the inflow time series in time

  ! compilation
  ! gfortran -g -O -c  nrtype.F90             -o nrtype.o
  ! gfortran -g -O -c  iso_varying_string.F90 -o iso_varying_string.o
  ! gfortran -g -O -c  ModuleFloodTools.F90   -o ModuleFloodTools.o
  ! gfortran -g -O -c  test_srch.F90          -o test_srch.o
  ! gfortran -o test_srch test_srch.o nrtype.o ModuleJGPFunctions.o ModuleGlobalData.o iso_varying_string.o

  USE nrtype
  !USE iso_varying_string
  !USE ModuleGlobalData
  USE ModuleFloodTools

  INTEGER                :: meqn = 3                     !1 number of equations in the matrix of conserved variables (I believe this is 5 in Geoclaw)
  INTEGER                :: mbc  = 2                     !n-pixels extension of the grid in domain boundaries
  INTEGER                :: mx   = 7                     !W-E number of pixels
  INTEGER                :: my   = 20                    !S-N number of pixels
  !REAL(DP)               :: xlower = 400.0               !ll corner of the domain W
  !REAL(DP)               :: ylower = 100.0               !ll corner of the domain S
  REAL(DP)               :: xlower = 373550.0            ! for the hydrodydamic domain of the severn
  REAL(DP)               :: ylower = 227925.0
  
  REAL(DP)               :: dx     = 15.0
  REAL(DP)               :: dy     = 15.0
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: q     ! = (/ 2.0, 5.0, 4.0 /) !initial state vector
  !INTEGER                :: maux   = 3                   !number of auxiliary variables (1st is topography) - not needed for srch
  REAL(DP)               :: t      = 0.0                !instantaneous time
  REAL(DP)               :: dt     =  5.0             !timestep
  TYPE(T_bci), DIMENSION(:), POINTER :: bci => null()
  TYPE(T_irts), POINTER  :: bdy => null()
  INTEGER                :: unit   = 11
  INTEGER                :: status, i, it

  CHARACTER(len=PathLength) :: bcifile = 'bci.dat'
  CHARACTER(len=PathLength) :: bdyfile = 'bdy.dat'
  CHARACTER(len=StringLength), DIMENSION(:), POINTER :: bcinames

  INTEGER               :: nt = 500                       ! advance 500 timesteps


  ! Read bci file
  CALL readBCI(bcifile, bci, unit, status)

  ALLOCATE(bcinames(SIZE(bci)))
  DO i=1, SIZE(bci)
     bcinames(i) = bci(i)%name
     WRITE(*,*)'bci: ',i
     WRITE(*,*) bci(i)
  END DO
  WRITE(*,*)'bci names: ',bcinames

  ! Read bdy file foo REAL array for now
  CALL readIRTS(bdyfile, bcinames, bdy, unit, status)

  !write timeseries input to the screen
  WRITE(*,*)'bdy:'
  CALL writeIRTS(bdy)

  ! init q
  ALLOCATE(q(meqn, 1-mbc:mx+mbc,1-mbc:my+mbc))
  q = 0.0
  q(1,1:mx,1:my) = 1.0_dp

  DO i=1-mbc,my+mbc
     WRITE(*,'(11F5.1)') q(1,:,i)
  END DO

  ! actual srch test
  DO it=1,2
   WRITE(*,*)'t :',t
   CALL srch(meqn, mbc, mx, my, xlower, ylower, dx, dy, t, dt, bci, bdy, q, verbose=1)
    t = t + dt
    DO i=1-mbc,my+mbc
       WRITE(*,'(11F9.5)') q(1,:,i)
    END DO
  END DO

  DEALLOCATE(bci)
  DEALLOCATE(bcinames)
  DEALLOCATE(bdy)

END PROGRAM test_srch
