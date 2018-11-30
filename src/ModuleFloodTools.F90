!--------------------------------------------
!        UoR, Hydrologic Modelling Group
!--------------------------------------------
!
!
! PROJECT       : DEMON, SINATRA
! MODULE        : FloodTools
! AUTHOR        : Javier Garcia-Pintado
! Last version  : v1.0
! Date          : 01/01/2015

MODULE ModuleFloodTools

  USE nrtype
  USE iso_varying_string

  !USE ModuleGlobalData

  IMPLICIT NONE

  PRIVATE

  !Subroutines-----------------------------------------------------------------

  !Constructor
  PUBLIC  :: ConstructInflows
  PRIVATE ::   ReadBCI                                     ! BCI input without the MOHID ModuleEnterData facilities
  PUBLIC  :: ReadTimes2IRTS
  PUBLIC  :: ReadIRTS                                      ! BDY input as multivariate irts object
  PUBLIC  :: ConstructKin1D
  PRIVATE :: kin1Dorder
  PUBLIC  :: readGmeta
  PRIVATE :: getHcPc

  !Selector
  PUBLIC  :: approxIRTS
  PUBLIC  :: getNNcell                                     ! get the nearest neigbour value from topography datasets for a specific location
  PUBLIC  :: writeGmeta
  PRIVATE :: getDdis
  PRIVATE :: getDslo0
  PRIVATE :: summaryAD ! 2D array DOUBLE
  PRIVATE :: summaryAI ! 2D array INTEGER
  PUBLIC  :: which     ! SUBROUTINE interface to which, including automatic allocation

  !Modifier
  PUBLIC  :: writeIRTS
  PUBLIC  :: srch
  PUBLIC  :: swap
  PUBLIC  :: sorti
  PUBLIC  :: kin1D

  !Destructor
  PUBLIC  :: DestructInflows
  PUBLIC  :: DestructKin1D

  !FUNCTIONS
  PUBLIC  :: str
  PUBLIC  :: ij2k    ! ELEMENTAL : [i,j] -> k [column-major oder vectorised 2D array]
  PUBLIC  :: k2ij
  PUBLIC  :: whichf   ! PACK LOGICAL into INTEGER index array
  PRIVATE :: AdjustSlope !ELEMENTAL

  !Interfaces
  INTERFACE approxIRTS
     MODULE PROCEDURE scalarApproxIRTS
     MODULE PROCEDURE vectorApproxIRTS
  END INTERFACE

  !Types
  !PUBLIC :: T_POSIXct
  !TYPE T_POSIXct
  !  REAL(DP), DIMENSION, POINTER :: 
  !END TYPE T_POSIXct

  PUBLIC :: T_bci                                         ! metadata for local mass inflows
  TYPE T_bci
    REAL                        :: w,e,s,n
    CHARACTER(len=1)            :: to   = '-'
    CHARACTER(len=StringLength) :: name = null_str
  END TYPE T_bci

  PUBLIC :: T_irts
  TYPE T_irts
    REAL(DP), DIMENSION(:),   POINTER                  :: time   => null()  ! [nt]
    REAL(DP), DIMENSION(:,:), POINTER                  :: value  => null()  ! [nt,m]
    CHARACTER(len=StringLength), DIMENSION(:), POINTER :: vname  => null()  ! [m] field names
  END TYPE T_irts

  PUBLIC :: T_gmeta
  TYPE T_gmeta
    CHARACTER(len=PathLength)   :: location = 'foo'
    CHARACTER(len=StringLength) :: mapset   = 'foo'
    CHARACTER(len=PathLength)   :: proj4    = ''
    REAL(DP)                    :: n, s, w, e        ! corners
    REAL(DP)                    :: nsres, ewres      ! resolution
    INTEGER(I4B)                :: rows, cols
    INTEGER(I4B)                :: cells
    REAL(DP), DIMENSION(2)      :: xlims             ! corners: eastings
    REAL(DP), DIMENSION(2)      :: ylims             ! corners: northings
    REAL(DP), DIMENSION(:), POINTER :: xseq => null() ! pixel center sequences
    REAL(DP), DIMENSION(:), POINTER :: yseq => null()
    REAL(DP), DIMENSION(:), POINTER :: ryseq => null()
  END TYPE T_gmeta

  PRIVATE :: T_SGCpar                              ! data.frame for subgrid channnel
  TYPE T_SGCpar
     INTEGER(I4B)                :: idx           ! index matching 'cmsk' classification
     INTEGER(I4B)                :: gtype         ! channel cross-section geometry
     REAL(DP)                    :: r             ! geometry parameters
     REAL(DP)                    :: p             !   "
     REAL(DP)                    :: s             !   "
     REAL(DP)                    :: n             ! manning coefficient
     REAL(DP)                    :: mc            ! meandering coefficient
  END TYPE T_SGCpar

  PRIVATE :: T_ldnetwork                          ! a local drainage network contains for a cell it corresponding upstream cells
  TYPE T_ldnetwork                                ! according to a local drainage direction map
    INTEGER(I4B), DIMENSION(:,:), POINTER :: up  => null()  ! [:,2] matrix with upstream pixel indexes
  END TYPE T_ldnetwork

  !Global Module Variables
  TYPE(T_bci), DIMENSION(:), POINTER, PUBLIC                 :: bci => null()
  TYPE(T_irts)             , POINTER, PUBLIC                 :: bdy => null()

  CHARACTER(len=StringLength), DIMENSION(:), POINTER, PUBLIC :: bcinames
  CHARACTER(len=PathLength)                          :: bcifile   = 'bci.data'
  CHARACTER(len=PathLength)                          :: bdyfile   = 'bdy.data'
  CHARACTER(len=StringLength)                        :: gmetafile = 'G.data'
  INTEGER                                            :: unit

  ! k01 [kinematic wave]
  ! k01 input
  TYPE(T_gmeta)           , POINTER, PUBLIC :: G   => null()
  TYPE(T_SGCpar), DIMENSION(:), POINTER     :: SGCpar => null()
  TYPE(T_ldnetwork), DIMENSION(:), POINTER  :: ldnet  => null()
  REAL(DP),     DIMENSION(:,:), ALLOCATABLE :: dem    ! [m] digital elevation model 
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: ldd
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: aflx   ! [-] accumulated flux in pixel number
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE :: cmsk   ! [-] channel classes [0 in overland]
  REAL(DP),     DIMENSION(:,:), ALLOCATABLE :: cwid   ! [m] channel width   [0 in overland]
  REAL(DP),     DIMENSION(:,:), ALLOCATABLE :: cdep   ! [m] channel depth   [0 in overland]
  REAL(DP)                                  :: dsl    ! outlet downstream slope

  ! k01 operation
  REAL(DP),     DIMENSION(:,:),   ALLOCATABLE  :: slo0   ! [-] bottom slope (kinematic) / surface slope (difussive apprx.)
  REAL(DP),     DIMENSION(:,:),   ALLOCATABLE  :: nman   ! [cmplx] Manning coefficient
  INTEGER(I4B), DIMENSION(:,:),   ALLOCATABLE  :: kinor  ! [-] packed (nk,2)) kinematic order cascade
  REAL(DP),     DIMENSION(:,:),   ALLOCATABLE  :: ddis   ! [m] downstream distance
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE  :: didx   ! [-] downstream pixel indexes
  INTEGER(I4B), DIMENSION(2)                   :: ouij   ! [i,j] outlet coordinates

  CONTAINS

  SUBROUTINE ConstructInflows(STAT)
    INTEGER, OPTIONAL, INTENT(OUT)           :: STAT

    !local
    INTEGER :: STAT_CALL, unit, i

    unit = 11                                               ! Warning : unit should be better automatically selected
    CALL ReadBCI(bcifile, bci, unit, STAT_CALL)
    IF (STAT_CALL /= SUCCESS_) STOP 'ModuleFloodTools - ConstructInflows - ERR01'

    ALLOCATE(bcinames(SIZE(bci)))
    DO i=1, SIZE(bci)
      bcinames(i) = bci(i)%name
    END DO
    IF (verbose > 0) WRITE(*,*)'bci names: ',bcinames

    CALL ReadIRTS(bdyfile, bcinames, bdy, unit, STAT_CALL)
    IF (STAT_CALL /= SUCCESS_) STOP 'ModuleFloodTools - ConstructInflows - ERR02'
    IF (verbose > 0) CALL writeIRTS(bdy)

    STAT = STAT_CALL

   END SUBROUTINE ConstructInflows

   SUBROUTINE DestructInflows(STAT)
     INTEGER, OPTIONAL, INTENT(OUT)           :: STAT

     INTEGER, DIMENSION(3) :: status

     status = 0

     DEALLOCATE(bci,      STAT=status(1))
     DEALLOCATE(bcinames, STAT=status(2))
     DEALLOCATE(bdy,      STAT=status(3))

     STAT = sum(status)

   END SUBROUTINE DestructInflows

   SUBROUTINE ReadBCI(fname, bci, unit, STAT)
     ! allocate the pointer BCI containing metadata of inflow boundary conditions
     ! and read input bci data file
     IMPLICIT NONE

     CHARACTER(len=PathLength), INTENT(IN) :: fname
     TYPE(T_bci), DIMENSION(:), POINTER    :: bci     ! allocated and filled within this subroutine
     INTEGER, INTENT(IN)                   :: unit
     INTEGER, INTENT(OUT), OPTIONAL        :: STAT

     ! local
     INTEGER             :: nbci = 0
     INTEGER             :: STAT_CALL
     CHARACTER(len=1000) :: auxstring
     INTEGER             :: line_length = 1000
     INTEGER             :: i

     NULLIFY(bci)

     OPEN(unit, file=fname, status='OLD', action='READ')
     !Count number of lines
     STAT_CALL = 0
     DO
        READ(unit,"(A)", IOSTAT=STAT_CALL) auxstring
        IF (len_trim(trim(auxstring)) > line_length) THEN
                    write(*,*) 'Maximum of ', line_length,' characters is supported.'
                    write(*,*) 'String: '//trim(auxstring)
                    write(*,*) 'File  : '//trim(adjustl(fname))
                    write(*,*) 'Line  : ', nbci + 1
                    STOP 'ModuleJGPFunctions - ReadBCI - ERR01'
        END IF
        IF (STAT_CALL==-1) EXIT
        nbci = nbci + 1
      END DO

    ! allocates bci
    IF (nbci .GT. 0) THEN
       ALLOCATE(bci(nbci))
    END IF

    ! copy data to bci pointer
    REWIND(unit)
    DO i = 1, nbci
       READ (unit, *, IOSTAT = STAT_CALL) bci(i)%w, bci(i)%e, bci(i)%s, bci(i)%n, bci(i)%to, bci(i)%name
    END DO

    STAT = STAT_CALL
    ! close unit
    CLOSE(unit)
   END SUBROUTINE ReadBCI

   SUBROUTINE ReadTimes2IRTS(fname, irts, unit, STAT, vname)
     ! initialises a T_irts to 0.0d0
     CHARACTER(len=PathLength),   INTENT(IN)  :: fname
     TYPE(T_irts), POINTER                    :: irts  
     INTEGER, INTENT(IN)                      :: unit
     INTEGER, INTENT(OUT), OPTIONAL           :: STAT
     CHARACTER(len=StringLength), DIMENSION(:), INTENT(IN), OPTIONAL :: vname

     ! local
     INTEGER(I4B)                :: nt = 0                           ! time series length
     INTEGER(I4B)                :: it, n, stat_call

     IF (ASSOCIATED(irts)) DEALLOCATE(irts)
     ALLOCATE(irts, STAT = stat_call)
     IF (stat_call /= SUCCESS_) STOP  'ModuleFloodTools - ReadTimes2IRTS - ERR01'
     NULLIFY(irts%time)
     NULLIFY(irts%value)
     NULLIFY(irts%vname)

     IF (PRESENT(vname)) THEN
      n = SIZE(vname)
     ELSE
      n = 1
     END IF

     ! count number of records in input file
     OPEN(unit, file=fname, status='OLD', action='READ')
     DO
       READ(unit, *, IOSTAT=stat_call)
       IF (stat_call==-1) EXIT
       nt = nt + 1
     END DO
     REWIND(unit)

     ! allocate irts and read data
     ALLOCATE(irts%time(nt), STAT=stat_call)
     ALLOCATE(irts%value(nt,n))
     ALLOCATE(irts%vname(n))

     DO it = 1, nt
       READ(unit, *, IOSTAT=stat_call) irts%time(it)
       IF (stat_call==-1) EXIT
     END DO
     CLOSE(unit)
     irts%value = 0.d0
     IF (PRESENT(vname)) THEN
       irts%vname = vname
     ELSE
       irts%vname = 'x'
     END IF

   END SUBROUTINE ReadTimes2IRTS

   SUBROUTINE ReadIRTS(fname, vname, irts, unit, STAT)
     CHARACTER(len=PathLength),   INTENT(IN)  :: fname
     CHARACTER(len=StringLength), DIMENSION(:), POINTER  :: vname ! variable names in the multivariate time series
     TYPE(T_irts), POINTER                  :: irts
     INTEGER, INTENT(IN)                    :: unit
     INTEGER, INTENT(OUT), OPTIONAL         :: STAT

     ! local
     INTEGER                     :: i, n                             ! number of variables in the multivariate time series
     INTEGER                     :: STAT_CALL
     CHARACTER(len=StringLength) :: tlabel                           ! time label in input header (later disregarded)
     CHARACTER(len=StringLength), DIMENSION(SIZE(vname))  :: lname   ! local names as read from input irts file
     INTEGER                     :: nt = 0                           ! number of records in irts

     n = SIZE(vname)

     OPEN(unit, file=fname, status='OLD', action='READ')
     ! Obtain number of variables and timeseries length
     READ (unit,*, IOSTAT = STAT_CALL) tlabel, lname

     !check names match input names
     DO i=1,n
        IF (lname(i) /= vname(i)) THEN
           WRITE(*,*)'mismatch in metadata and named fields in input irts, field: ',i
           STOP 'ModuleJGPFunctions - ReadIRTS - ERR01'
        END IF
     END DO

     IF (ASSOCIATED(irts)) DEALLOCATE(irts)

     ALLOCATE(irts, STAT = STAT_CALL)
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR02'
     NULLIFY(irts%time)
     NULLIFY(irts%value)
     NULLIFY(irts%vname)

     ! count number of records in input file
     DO
       READ(unit, *, IOSTAT=STAT_CALL)
       IF (STAT_CALL==-1) EXIT
       nt = nt + 1
     END DO
     REWIND(unit)

     ! allocate irts and read data
     ALLOCATE(irts%time(nt))
     ALLOCATE(irts%value(nt,n))
     ALLOCATE(irts%vname(n))

     READ(unit, *, IOSTAT=STAT_CALL)
     DO i = 1,nt
        READ(unit, *, IOSTAT=STAT_CALL) irts%time(i), irts%value(i,:)
                IF (STAT_CALL==-1) EXIT
     END DO
     CLOSE(unit)

     irts%vname = vname

   END SUBROUTINE ReadIRTS

   SUBROUTINE scalarApproxIRTS(itsi, to, itso)
     TYPE(T_irts), POINTER     :: itsi    ! irts input object
     REAL(DP), INTENT(IN)      :: to      ! output times
     TYPE(T_irts), POINTER     :: itso    ! irts output object

     !local automatic
     REAL(DP), DIMENSION(1)    :: tov
     tov = to
     CALL vectorApproxIRTS(itsi, tov, itso)
   END SUBROUTINE scalarApproxIRTS

   SUBROUTINE vectorApproxIRTS(itsi, to, itso)
     TYPE(T_irts), POINTER               :: itsi     ! irts input object
     REAL(DP), DIMENSION(:), INTENT(IN)  :: to      ! output times
     TYPE(T_irts), POINTER               :: itso     ! irts output object

     ! local
     INTEGER :: i, n, nti, nto, it
     INTEGER :: istat
     REAL(DP), DIMENSION(SIZE(itsi%time)-1) :: dtime
     !INTEGER, DIMENSION(SIZE(to))       :: i1,i2       ! bound index arrays
     !REAL,    DIMENSION(SIZE(to))       :: t1,t2       ! bound time and weight arrays
     INTEGER                            :: i1, i2
     REAL(DP)                           :: t1, t2, w1, w2

     n   = SIZE(itsi%value,2)
     nti = SIZE(itsi%time)
     nto = SIZE(to)

     IF (.NOT. ASSOCIATED(itsi)) THEN
       STOP 'approxIRTS - ERR01'
     ELSE IF (ASSOCIATED(itso)) THEN
       DEALLOCATE(itso)
       ! STOP 'approxIRTS - ERR02'
     END IF

     ALLOCATE(itso, STAT=istat)
     ALLOCATE(itso%time(nto))
     ALLOCATE(itso%value(nto,n))
     ALLOCATE(itso%vname(n))

     dtime = itsi%time(2:nti) - itsi%time(1:nti-1)
     IF (ANY(dtime < 0.0)) THEN
        STOP 'approxIRTS - ERR03: non-monotonic increments in input timeseries'
     END IF

     itso%vname = itsi%vname
     itso%time  = to

     !compute interpolation
     DO i=1,nto
        i1 = 1
        i2 = 1
        DO it = 1,nti
           IF (itsi%time(it) >= to(i)) EXIT
           i1 = it
           IF (it < nti) i2 = it + 1
        END DO
        t1 = itsi%time(i1)
        t2 = itsi%time(i2)
       ! WRITE(*,*)'i1, i2: ',i1,',',i2
       ! WRITE(*,*)'t1, t2: ',to(i),',',t1,',',t2

        IF (i1 == i2) THEN
           itso%value(i,:) = itsi%value(i1,:)
        ELSE
           w1 = abs((to(i) - t2) / (t1 - t2))
           w2 = abs((to(i) - t1) / (t1 - t2))
           itso%value(i,:) = w1 * itsi%value(i1,:) + w2 * itsi%value(i2,:)
        END IF
     END DO
   END SUBROUTINE vectorApproxIRTS

   SUBROUTINE writeIRTS(its)
      TYPE(T_irts), POINTER           :: its     ! irts input object

      INTEGER :: i,nt,n
      CHARACTER(len=7) :: fmt                    ! max 99 fields in output timeseries allowed

      nt = SIZE(its%time)
      n  = SIZE(its%value,2)

      fmt = '(' // trim(str(n+1)) // 'A20)'
      WRITE(*,fmt) 'time : ',its%vname
      DO i = 1,nt
         WRITE(*,*), its%time(i),':',its%value(i,:)
      END DO
   END SUBROUTINE writeIRTS

SUBROUTINE srch(meqn, mbc, mx, my, xlower, ylower, dx, dy, t, dt, bci, bdy, q)
  ! subroutine to include a source of external mass into clawpack
  ! the current assumption is that the source of mass has a velocity
  ! according to manning's equation and a parameterized slope and manning coefficient

  IMPLICIT NONE

  ! bdy are the timeseries of mass inflows

  ! input
  INTEGER,  INTENT(IN)        :: meqn, mbc, mx, my
  REAL(DP), INTENT(IN)        :: xlower, ylower, dx, dy, t, dt       ! xlower, y lower are the corners of the lower-left pixel in the input computational grid
  TYPE(T_bci), DIMENSION(:), POINTER :: bci
  TYPE(T_irts),              POINTER :: bdy

  ! inout
  REAL(DP), INTENT(INOUT)     :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)   ! automatic array

  ! local
  INTEGER                     :: nbci, ib, i, j
  REAL(DP)                    :: bciarea, xupper, yupper
  !REAL, DIMENSION(mx)         :: xseq                                      ! automatic arrays
  !REAL, DIMENSION(my)         :: yseq, ryseq
  INTEGER                     :: wgi, egi, sgi, ngi, nxbci, nybci
  REAL(DP), DIMENSION(2)      :: tdt
  TYPE(T_irts),               POINTER :: sbdy => null()        ! to store interpolated timeseries at t,t+dt
  REAL(DP)                    :: wins, dw, dh                  ! instantaneous flow, and increment in water volume and level for the timestep & bci
  REAL(DP), PARAMETER         :: inflow_tol = 1.0D-3           ! minimum inflow rate [m3/s] to be taken into account
  REAL(DP)                    :: dhu,dhv                       ! increment in momentum for the timestep & bci
  REAL(DP)                    :: hin = 1.D0                    ! warning this makes the strong assumtions that inflow comes through a 1-m depth water column

  REAL(DP), PARAMETER         :: inmann = 0.035                ! inflow assumption: manning coefficient
  REAL(DP), PARAMETER         :: inslp  = 0.05                 ! inflow assumption: slope
  REAL(DP)                    :: inhu                          ! new inflow momentum

  tdt = (/ t, t+dt /)

  !xseq   = xlower + (/(i, i=1,mx)/) * dx - 0.5d0 * dx
  !yseq   = ylower + (/(j, j=1,my)/) * dy - 0.5d0 * dy
  xupper = xlower + mx * dx
  yupper = ylower + my * dy

  IF (verbose > 1) THEN
    WRITE(*,*)'srch - verbosity level:',verbose
    WRITE(*,*) 'll   : ', xlower, ylower
    WRITE(*,*) 'ur   : ', xupper, yupper
    WRITE(*,*) 't t+dt  : ', tdt
  END IF

  nbci = SIZE(bci)

  ! warning: assumption: dt < available timesteps in input time
  CALL approxIRTS(bdy,tdt,sbdy)   ! get instant inflow by linear interpolation in time at the bounds of the time interval

  bcido: DO ib = 1, nbci
    ! map area into grid
    IF (bci(ib)%w < xlower .OR. bci(ib)%e > xupper .OR. &             ! warning: user should check the inflow location are either point-located or within a computational subgrid
        bci(ib)%s < ylower .OR. bci(ib)%n > yupper) CYCLE

    wgi = NINT((bci(ib)%w - (xlower + dx / 2)) / dx) + 1
    egi = NINT((bci(ib)%e - (xlower + dx / 2)) / dx) + 1
    sgi = NINT((bci(ib)%s - (ylower + dy / 2)) / dy) + 1
    ngi = NINT((bci(ib)%n - (ylower + dy / 2)) / dy) + 1

    ! get mapped area, dh, an update water level at the inflow location
    nxbci = abs(wgi - egi) + 1
    nybci = abs(ngi - sgi) + 1
    bciarea = nxbci * nybci * dx * dy
    wins = SUM(sbdy%value(:,ib)) / 2.d0

    dhif: IF (wins >= inflow_tol) THEN
      dw = wins * dt
      dh = dw / bciarea                                           ! mass input (assumed density=1) m3

      IF (verbose > 0) THEN
        CALL writeIRTS(sbdy)
        WRITE(*,*)'bciarea, t, dt, dw, dh: ',bciarea, t, dt, dw, dh
      END IF

      inhu = dh**(5.d0/3.d0) * sqrt(inslp) / inmann

      SELECT CASE (bci(ib)%to)
        CASE ('n')
          dhu = 0.d0
          dhv = inhu
          !dhv =   dh * wins / (dx * nxbci * hin)                            ! dh * v
        CASE ('s')
          dhu = 0.d0
          dhv = - inhu
          !dhv = - dh * wins / (dx * nxbci * hin)
        CASE ('e')
          !dhu =   dh * wins / (dy * nybci * hin)
          dhu = inhu
          dhv = 0.d0
        CASE ('w')
          !dhu = - dh * wins / (dy * nybci * hin)
          dhu = - inhu
          dhv = 0.d0
        CASE DEFAULT
          dhu = 0.d0
          dhv = 0.d0
      END SELECT
      !dhu = 0.d0
      !dhv = 0.d0
      DO j = sgi,ngi
        DO i = wgi,egi
          q(1,i,j) = q(1,i,j) + dh
          q(2,i,j) = q(2,i,j) + dhu
          q(3,i,j) = q(3,i,j) + dhv
        END DO
      END DO
    END IF dhif

  END DO bcido

END SUBROUTINE srch

SUBROUTINE getNNcell(topoint, xim,xcell,xip,yjm,ycell,yjp,            &
                     xlowtopo,ylowtopo,xhitopo,yhitopo,dxtopo,dytopo, &
                     mxtopo,mytopo,mtopo,i0topo,mtopoorder,topo)
  !====================================================================
  !
  !   cellgetnn get the nearest neighbour value
  !     defined from data from multiple regular Cartesian grids
  !     (using the finest data available in any region)
  !
  !     The rectangle has coords:
  !     xim <= x <= xip, yjm <= y <= yjp, with center (x,y) = (xcell, ycell)
  !
  !     The intersection (with one particular grid has coords:
  !     xintlo <= x <= xinthi, yintlo <= y <= yinthi
  !

      IMPLICIT NONE

      !input
      REAL(DP), INTENT(OUT) :: topoint                                         ! value to return: sample from the topo datasets

      REAL(DP), INTENT(IN)  :: xim,xcell,xip,yjm,ycell,yjp                     ! coordinates of the pixel center to be queried [xcell,ycell] plus borders of pixel
      REAL(DP), DIMENSION(:), INTENT(IN) :: xlowtopo, ylowtopo                 ! len=mtopofiles : centers of lower-left pixels for every topography dataset
      REAL(DP), DIMENSION(:), INTENT(IN) :: xhitopo, yhitopo                   ! len=mtopofiles : centers of upper-right pixels for every topography dataset
      REAL(DP), DIMENSION(:), INTENT(IN) :: dxtopo, dytopo                     ! len=mtopofiles : horizontal, vertical resolution for each topography dataset
      INTEGER,  DIMENSION(:), INTENT(IN) :: mxtopo, mytopo, mtopo              ! len=mtopofiles : number of pixels in topography dataset in x,y directions, and total
      INTEGER,  DIMENSION(:), INTENT(IN) :: i0topo                             ! len=mtopofiles :
      INTEGER,  DIMENSION(:), INTENT(IN) :: mtopoorder                         ! len=mtopofiles :
      REAL(DP), DIMENSION(:), INTENT(IN) :: topo                               ! len=mtoposize  :

      ! local
      INTEGER  :: m, mfid, i0, iid, jid, pix
      INTEGER  :: mtopofiles

      mtopofiles = SIZE(xlowtopo)

      topoint = 0.d0                                                           ! initialize the integral of the surface

      !IF ((xlowtopo(1) + (mxtopo(1) - 1) * dxtopo(1) .NE. xhitopo(1)) .OR. &
      !    (ylowtopo(1) + (mytopo(1) - 1) * dytopo(1) .NE. yhitopo(1))) THEN
      !   WRITE(*,*) 'getNNcell:: error in input topography'
      !   STOP
      !END IF

      DO m = 1,mtopofiles
         !look at topofiles, from fine to coarse
         mfid = mtopoorder(m)
         i0 = i0topo(mfid)
         !check grid cell is within the bbox of this topofile
         IF ( xcell >= xlowtopo(mfid) - dxtopo(mfid) .AND. &
              xcell <= xhitopo(mfid)  + dxtopo(mfid) .AND. &
              ycell >= ylowtopo(mfid) - dytopo(mfid) .AND. &
              ycell <= yhitopo(mfid)  + dytopo(mfid)) THEN

            iid = NINT((xcell - xlowtopo(mfid)) / dxtopo(mfid)) + 1
            jid = NINT((yhitopo(mfid) - ycell) / dytopo(mfid))+1
            pix   = (jid - 1) * mxtopo(mfid) + iid
            topoint = topo(i0-1+pix)
            RETURN
         END IF
      END DO

      write(6,601) xcell,ycell
 601  format('*** Error, grid cell center not within any topo grid',/, &
             '  xcell = ',e24.14,'  ycell = ',e24.14)
      stop

END SUBROUTINE getNNcell

CHARACTER(len=20) FUNCTION str(x)
   !   "Convert an integer to string."
   INTEGER, intent(in) :: x
   WRITE (str, *) x
   str = adjustl(str)
END FUNCTION str


SUBROUTINE swap(a,b)
  IMPLICIT NONE

  INTEGER(I4B), INTENT(INOUT) :: a, b
  INTEGER(I4B)                :: tmp

  tmp = a
  a   = b
  b   = tmp

END SUBROUTINE swap

SUBROUTINE sorti(x)
  ! sort a 1D array
  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:), INTENT(INOUT) :: x
  ! local
  INTEGER(I4B) :: i, n, loc

  n = SIZE(x)
  DO i = 1, n-1   ! except for the last
    loc = i - 1 + MINLOC(x(i:n), DIM=1)
    CALL swap(x(i), x(loc))
  END DO

END SUBROUTINE sorti

SUBROUTINE ConstructKin1D(fk01p, mask, STAT)
  ! Construct overland flow routing with a 1D kinematic wave approximation
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN) :: fk01p         !file with global data
    LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask   ![nx,ny] matrix
    INTEGER, OPTIONAL, INTENT(OUT)      :: STAT

    !local
    REAL(DP)     :: dx
    INTEGER(I4B) :: STAT_CALL, unit, ic, ikin
    INTEGER(I4B), DIMENSION(11) :: status
    INTEGER(I4B)                :: nk         ! number of pixels within the valid domain
    CHARACTER(len=StringLength) :: dem_f, ldd_f, aflx_f, cmsk_f, cwid_f, cdep_f, SGCpar_f
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: ouk
    INTEGER(I4B)                :: nx, ny !grid dimension nx: eastings, ny: northings

    unit = 11                                ! input unit: TODO automate selection
    nk = COUNT(mask)
    nx = SIZE(mask,1)
    ny = SIZE(mask,2)

    ! input
    ALLOCATE(dem(nx,ny),  STAT=status(1))    ! local drainage direction PCRaster codes: 7 8 9
    ALLOCATE(ldd(nx,ny),  STAT=status(2))    ! local drainage direction PCRaster codes: 7 8 9
    ALLOCATE(aflx(nx,ny), STAT=status(3))    ! accumulated flux pixels [-]              4 5 6
    ALLOCATE(cmsk(nx,ny), STAT=status(4))    !                                          1 2 3
    ALLOCATE(cwid(nx,ny), STAT=status(5))    !
    ALLOCATE(cdep(nx,ny), STAT=status(6))

    ! operation
    ALLOCATE(slo0(nx,ny),   STAT=status(7))    ! pixel bottom slope      [0/1]
    ALLOCATE(nman(nx,ny),   STAT=status(8))    ! manning coefficient     [-] [replaces nsur in ModuleHspfWater for distributed routing]
    ALLOCATE(kinor(nk,2),   STAT=status(9))    ! packed (i,j) indexes for kinematic cascade
    ALLOCATE(ddis(nx,ny),   STAT=status(10))   ! downstream distance
    ALLOCATE(didx(2,nx,ny), STAT=status(11))   ! downstream pixels

    !read main k01 input file
    OPEN(UNIT=unit, FILE=fk01p, STATUS='OLD', ACTION='READ', &
       IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit,*) dem_f
       READ(unit,*) ldd_f
       READ(unit,*) aflx_f
       READ(unit,*) cmsk_f
       READ(unit,*) cwid_f
       READ(unit,*) cdep_f
       READ(unit,*) SGCpar_f
       READ(unit,*) dsl                                 ! module scope
       CLOSE(unit)
    END IF
    !read geographical data
    WRITE(*,*) 'ModuleFloodToolds::ConstructKin1D::dem_f:',dem_f
    OPEN(UNIT=unit, FILE=dem_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) dem
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open dem ERR'
    END IF
    
    ! forced stop
    !WRITE(*,*)'dem(10,:):',dem(10,:) 
    !STOP 'ModuleFloodTools - forced STOP 001'

    OPEN(UNIT=unit, FILE=ldd_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) ldd
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open ldd ERR'
    END IF

    OPEN(UNIT=unit, FILE=aflx_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) aflx
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open aflx ERR'
    END IF
    OPEN(UNIT=unit, FILE=cmsk_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) cmsk
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open cmsk ERR'
    END IF
    OPEN(UNIT=unit, FILE=cwid_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) cwid
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open cwid ERR'
    END IF
    OPEN(UNIT=unit, FILE=cdep_f, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    IF ( STAT_CALL == SUCCESS_ ) THEN
       READ(unit) cdep
       CLOSE(unit)
    ELSE
       STOP 'ModuleFloodTools - open cdep ERR'
    END IF

    WHERE (.NOT. mask)
      ldd   = 0
      aflx  = 0
      cmsk  = 0       ! as floodplain
      cwid  = 0.D0
      cdep  = 0.D0
    END WHERE

    CALL readGmeta(gmetafile, G, unit, STAT_CALL)
    IF (STAT_CALL /= SUCCESS_) STOP 'ModuleFloodTools - ConstructKin1D - ERR01'
    CALL writeGmeta(G)

    CALL readSGCpar(SGCpar_f, SGCpar, unit, STAT_CALL)
    IF (STAT_CALL /= SUCCESS_) STOP 'ModuleFloodTools - ConstructKin1D - ERR01'
    CALL writeSGCpar(SGCpar)


    !outlet indexes
    CALL which(RESHAPE(ldd,[G%cells]) == 5, ouk)
    WRITE(*,*) 'ouk:',ouk
    IF (SIZE(ouk) > 1) STOP 'ConstructKin1D :: ldd > 1 oulet'
    ouij = RESHAPE(k2ij(ouk,G%cols),SHAPE(ouij))
    DEALLOCATE(ouk)
    IF (verbose > 0) THEN
      !WRITE(*,*) 'ouij:',ouij,' | east:',G%xseq(ouij(1,1)),',north:',G%ryseq(ouij(1,2))
      WRITE(*,*) 'ouij:',ouij,' | east:',G%xseq(ouij(1)),',north:',G%ryseq(ouij(2))
    
      WRITE(*,*) 'ldd(ouij):', ldd(ouij(1),ouij(2))


      !WRITE(*,*) 'summary(mask)'
      !CALL summaryAI((mask .EQV. .TRUE.) * 1)

      WRITE(*,*) NEW_LINE('A'),'summary(cmsk)'
      CALL summaryAI(cmsk)
    END IF
    ! expand manning coefficients as distributed
    nman = 0.d0
    DO  ic = 1,SIZE(SGCpar)
      WHERE(cmsk == SGCpar(ic)%idx)
        nman = SGCpar(ic)%n
      END WHERE
    END DO
    IF (COUNT(nman == 0.d0) > 0) THEN
      STOP 'ModuleFloodToold - ERR unassigned pixel in manning expansion'
    END IF
    IF (verbose > 0) THEN
      WRITE(*,*) NEW_LINE('A'),'summary(nman, mask)'
      CALL summaryAD(nman, mask)
    END IF
    !WRITE(*,*) 'dem(150,:)', dem(150,:)      !         N-S section at pixel 150
    !WRITE(*,*) 'ldd(150,:)', ldd(150,:)      !         N-S section at pixel 150
    !WRITE(*,*) 'aflx(150,:)', aflx(150,:)    !         N-S section at pixel 150
    !WRITE(*,*) 'cmsk(150,40:100)', cmsk(150,40:100)    !         N-S transect
    !WRITE(*,*) 'nman(150,40:100)', nman(150,40:100)    !         N-S transect
    !OPEN(UNIT=15, FILE='nman_out.for', &
    !   STATUS='REPLACE', ACTION='WRITE', &
    !   FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
    !   WRITE(15) nman
    !   CLOSE(15)
    ! into R:
    ! cc <- readBin(file.path(dsnsim,'input','00001','nman_out.for'), n=G$cells, what='double')
    ! tmpSP@data[,1] <- cc; plotGmeta(layer=tmpSP, col=matlab.like(100))
    ! WRITE(*,*) 'cwid(150,:)', cwid(150,:)    !         N-S section at pixel 100
    !WRITE(*,*) 'cdep(150,:)', cdep(150,:)    !         N-S section at pixel 100
    !WRITE(*,*) 'DEM(i=150,j=100):',dem(150,100)
    ! i+(j-1)*nx

    !calculate optimal k01 routing order
    CALL kin1Dorder(aflx, kinor)

    IF (verbose > 0) THEN
      WRITE(*,*) NEW_LINE('A'),'head(kinor) :'
      DO ikin=1,5
        WRITE(*,*) kinor(ikin,:), ij2k(kinor(ikin,1),kinor(ikin,2),G%cols)
      END DO
      WRITE(*,*) NEW_LINE('A'),'tail(kinor) :'
      DO ikin=nk-4,nk
        WRITE(*,*) kinor(ikin,:), ij2k(kinor(ikin,1),kinor(ikin,2),G%cols)
      END DO
    END IF
    !WRITE(*,*) ij2k(kinor(1:10,1),kinor(1:10,2),G%cols)
    !WRITE(*,*) 'test k2ij:'
    !WRITE(*,*) 'k2ij(286)', k2ij((/286/), G%cols)
    !WRITE(*,*) 'k2ij(kinor(1:10,:))', k2ij( ij2k(kinor(1:10,1),kinor(1:10,2),G%cols),G%cols)

    CALL getLdnet(kinor, ldd, ldnet)

    CALL getDdis(ldd, G%nsres, ddis)
    IF (verbose > 0) THEN
      !WRITE(*,*) 'ddis(150,20:100)', ddis(150,20:100)    !         N-S section at pixel 150
      WRITE(*,*) NEW_LINE('A'),'summary(ddis,mask)'
      CALL summaryAD(ddis,mask)
    END IF

    CALL getLDdidx(ldd,didx)
    !WRITE(*,*) 'summary(didx,mask)'
    !CALL summaryAI(didx,mask)

    CALL getDslo0(dem,ldd,ddis,dsl,slo0)
    IF (verbose > 0) THEN
      WRITE(*,*) NEW_LINE('A'),'summary(slo0, mask)'
      CALL summaryAD(slo0, mask)
      WRITE(*,*) NEW_LINE('A'),'ddis(ouij):', ddis(ouij(1),ouij(2))
      WRITE(*,*)               'slo0(ouij):', slo0(ouij(1),ouij(2))
    END IF
    IF (sum(status) == 0) THEN
      STAT_CALL = 0
    END IF
    IF (STAT_CALL /= SUCCESS_) STOP 'ModuleFloodTools - ConstructKin1D - ERR01'

    STAT = STAT_CALL

END SUBROUTINE ConstructKin1D

SUBROUTINE DestructKin1D(STAT)
     INTEGER, OPTIONAL, INTENT(OUT)           :: STAT

     INTEGER, DIMENSION(14) :: status

     status = 0

     DEALLOCATE(dem,   STAT=status(1))
     DEALLOCATE(ldd,   STAT=status(2))
     DEALLOCATE(aflx,  STAT=status(3))
     DEALLOCATE(cmsk,  STAT=status(4))
     DEALLOCATE(cwid,  STAT=status(5))
     DEALLOCATE(cdep,  STAT=status(6))

     DEALLOCATE(slo0,  STAT=status(7))
     DEALLOCATE(nman,  STAT=status(8))
     DEALLOCATE(kinor, STAT=status(9))
     DEALLOCATE(ddis,  STAT=status(10))
     DEALLOCATE(didx,  STAT=status(11))

     ! DEALLOCATE pointers
     DEALLOCATE(G,      STAT=status(12))
     DEALLOCATE(SGCpar, STAT=status(13))
     DEALLOCATE(ldnet,  STAT=status(14))

     STAT = sum(status)

END SUBROUTINE DestructKin1D

SUBROUTINE kin1Dorder(accuflux, kinor)
  !=====================
  !
  ! obtain a pixel simulation order for a upstream-downstream simulation cascade with 1D kinematic wave
  ! the assumption is that pixels in a matrix are numerated by rows starting from the upper-left corner and the advancing by rows
  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:,:), INTENT(IN)  :: accuflux
  INTEGER(I4B), DIMENSION(:,:), INTENT(OUT) :: kinor

  ! local variables
  INTEGER(I4B) :: nx,ny,n,i,j, k, acval, nlev, nl
  INTEGER(I4B), DIMENSION(SIZE(accuflux))   :: accsrt, accdif
  INTEGER(I4B), DIMENSION(SIZE(accuflux,1),SIZE(accuflux,2)) :: kinorm
  !INTEGER(I4B), DIMENSION(SIZE(accuflux)-1) :: accsr
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: acclev

  nx = SIZE(accuflux,1)
  ny = SIZE(accuflux,2)
  n  = nx*ny

  accsrt = RESHAPE(accuflux,[n])
  !WRITE(*,*) 'accuflux 1D:',accsrt
  CALL sorti(accsrt)
  !WRITE(*,*) 'accsrt:',accsrt
  accdif = (/ 1, accsrt(2:n) - accsrt(1:n-1) /)
  IF (ANY(accuflux==0)) THEN
    accdif(1) = 0
  END IF
  nlev = COUNT(accdif /= 0)
  !WRITE(*,*) 'nlev:',  nlev
  ALLOCATE(acclev(nlev))
  acclev = PACK(accsrt,accdif /= 0,acclev)  ! accumulated influx classes

  !WRITE(*,*) 'acclev:',acclev

  k = 0
  kinorm = 0
  DO nl = 1,nlev
    northDO: DO j = 1,ny
      eastDO: DO i = 1,nx
        IF (accuflux(i,j) == acclev(nl)) THEN
          k = k + 1
          kinorm(i,j) = k
        END IF
      END DO eastDO
    END DO northDO
  END DO
  DEALLOCATE(acclev)

  !WRITE(*,*) '1D kinematic cascade order:'
  !DO j=1,ny               ! geographical display
  !  WRITE(*,*)  kinorm(:,j)
  !END DO

  DO k=1,MAXVAL(kinorm)
   kinor(k,:) = MINLOC(abs(kinorm - k))
  END DO

END SUBROUTINE kin1Dorder


SUBROUTINE readGmeta(fname, G, unit, STAT)
  IMPLICIT NONE

  CHARACTER(len=StringLength),   INTENT(IN)  :: fname
  TYPE(T_gmeta), POINTER                   :: G
  INTEGER(I4B), INTENT(IN)                 :: unit
  INTEGER(I4B), INTENT(OUT), OPTIONAL      :: STAT

  ! local
  INTEGER(I4B)                           :: STAT_CALL, i, j


  NULLIFY(G)
  ALLOCATE(G, STAT = STAT_CALL)
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR02'
  NULLIFY(G%xseq)
  NULLIFY(G%yseq)
  NULLIFY(G%ryseq)

  OPEN(unit, file=fname, status='OLD', action='READ')
  ! Obtain number of variables and timeseries length
  READ(unit,*, IOSTAT = STAT_CALL) G%location
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR03'
  READ(unit,*, IOSTAT = STAT_CALL) G%mapset
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR04'
  READ(unit,*, IOSTAT = STAT_CALL) G%proj4
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR05'
  READ(unit,*, IOSTAT = STAT_CALL) G%n, G%s, G%w, G%e
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR06'
  READ(unit,*, IOSTAT = STAT_CALL) G%nsres, G%ewres
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR07'
  READ(unit,*, IOSTAT = STAT_CALL) G%rows, G%cols
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR08'
  READ(unit,*, IOSTAT = STAT_CALL) G%cells
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR09'
  READ(unit,*, IOSTAT = STAT_CALL) G%xlims
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR10'
  READ(unit,*, IOSTAT = STAT_CALL) G%ylims
     IF (STAT_CALL /= SUCCESS_) STOP  'ModuleFloodTools - ReadIRTS - ERR11'
  CLOSE(unit)

  ALLOCATE(G%xseq(G%cols))     ! W -> E
  ALLOCATE(G%yseq(G%rows))     ! S -> N
  ALLOCATE(G%ryseq(G%rows))    ! N -> S

  G%xseq = G%w + (/(i, i=1,G%rows)/) * G%ewres - 0.5d0 * G%ewres
  G%yseq = G%s + (/(j, j=1,G%rows)/) * G%nsres - 0.5d0 * G%nsres
  G%ryseq = G%yseq(G%rows:1:-1)

END SUBROUTINE readGmeta

SUBROUTINE writeGmeta(G)
  IMPLICIT NONE

  TYPE(T_gmeta), POINTER                   :: G

    WRITE(*,*) ''
    WRITE(*,*) 'Gmeta ::'
    WRITE(*,*) 'G%location:', trim(G%location)
    WRITE(*,*) 'G%mapset  :', trim(G%mapset)
    WRITE(*,*) 'G%proj4   :', trim(G%proj4)
    WRITE(*,*) 'G%n       :', G%n
    WRITE(*,*) 'G%s       :', G%s
    WRITE(*,*) 'G%w       :', G%w
    WRITE(*,*) 'G%e       :', G%e
    WRITE(*,*) 'G$nsres   :', G%nsres
    WRITE(*,*) 'G%ewres   :', G%ewres
    WRITE(*,*) 'G%rows    :', G%rows
    WRITE(*,*) 'G%cols    :', G%cols
    WRITE(*,*) 'G%cells   :', G%cells
    WRITE(*,*) 'G%xlims   :', G%xlims
    WRITE(*,*) 'G%ylims   :', G%ylims
    WRITE(*,*) 'G%xseq    :', G%xseq(1:3),'...'
    WRITE(*,*) 'G%yseq    :', G%yseq(1:3),'...'
    WRITE(*,*) 'G%ryseq   :', G%ryseq(1:3),'...'

END SUBROUTINE writeGmeta

SUBROUTINE readSGCpar(fname, p, unit, STAT)
  ! allocate and read SGCpar [subgrid channel parameters] structure
  IMPLICIT NONE

  CHARACTER(len=StringLength),   INTENT(IN)  :: fname
  TYPE(T_SGCpar), DIMENSION(:), POINTER      :: p
  INTEGER(I4B), INTENT(IN)                   :: unit
  INTEGER(I4B), INTENT(OUT), OPTIONAL        :: STAT

  ! local
  INTEGER(I4B)  :: STAT_CALL, i
  INTEGER(I4B)  :: nch = 0

  NULLIFY(p)

  !count input lines
  OPEN(unit, file=fname, status='OLD', action='READ')
  READ(unit,*, IOSTAT=STAT_CALL)               ! header
  DO
    READ(unit,*, IOSTAT=STAT_CALL)             ! 1 line / channel class
    IF  (STAT_CALL==-1) EXIT
    nch = nch + 1
  END DO
  WRITE(*,*) 'readSGCpar: nch=',nch

  IF (nch > 0) THEN
    ALLOCATE(p(nch))
  END IF
  ! read data into pointer
  REWIND(unit)
  READ(unit,*, IOSTAT=STAT_CALL)               ! header
  DO i = 1, nch
     READ (unit, *, IOSTAT = STAT_CALL) p(i)%idx, p(i)%gtype, p(i)%r, p(i)%p, &
           p(i)%s, p(i)%n, p(i)%mc
  END DO
  STAT = STAT_CALL
  ! close unit
  CLOSE(unit)
END SUBROUTINE readSGCpar

SUBROUTINE writeSGCpar(p)
  IMPLICIT NONE

  TYPE(T_SGCpar), DIMENSION(:), POINTER  :: p


  ! local
  INTEGER(I4B)    :: nch, i

  nch = SIZE(p)

  WRITE(*,*) ''
  WRITE(*,*) 'SGCpar ::'
  WRITE(*,*) 'idx gtype r p s n mc'
  DO i=1,nch
    WRITE(*,*) p(i)%idx, p(i)%gtype, p(i)%r, p(i)%p, &
               p(i)%s, p(i)%n, p(i)%mc
  END DO

END SUBROUTINE writeSGCpar

SUBROUTINE summaryAD(x2D, mask)
 IMPLICIT NONE

 REAL(DP), DIMENSION(:,:), INTENT(IN) :: x2D
 LOGICAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: mask

 !local
 INTEGER(I4B) :: n, nx, ny
 REAL(DP)     :: xmin, xbar, xmax

 nx = SIZE(x2D,1)
 ny = SIZE(x2D,2)

 IF (.NOT. PRESENT(mask)) THEN
   n = nx * ny
   xmin = MINVAL(x2D)
   xbar = SUM(x2D) / (nx * ny)
   xmax = MAXVAL(x2D)
 ELSE
   n = COUNT(mask)
   xmin = MINVAL(x2D, MASK=mask)
   xbar = SUM(x2D, MASK=mask) / COUNT(mask)
   xmax = MAXVAL(x2D, MASK=mask)
 END IF


 WRITE(*,*) 'head : ',x2D(1:MIN(5,nx),1)
 WRITE(*,*) 'tail : ',x2D((nx-MIN(4,nx-1)):nx,ny)
 WRITE(*,*) 'n min mean max'
 WRITE(*,*) n, xmin, xbar, xmax

END SUBROUTINE summaryAD

SUBROUTINE summaryAI(x2D, mask)
 IMPLICIT NONE

 INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: x2D
 LOGICAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: mask

 !local
 INTEGER(I4B) :: n, nx, ny, xmin, xbar, xmax

 nx = SIZE(x2D,1)
 ny = SIZE(x2D,2)

 IF (.NOT. PRESENT(mask)) THEN
   n = nx * ny
   xmin = MINVAL(x2D)
   xbar = SUM(x2D) / (nx * ny)
   xmax = MAXVAL(x2D)
 ELSE
   n = COUNT(mask)
   xmin = MINVAL(x2D, MASK=mask)
   xbar = SUM(x2D, MASK=mask) / COUNT(mask)
   xmax = MAXVAL(x2D, MASK=mask)
 END IF

 WRITE(*,*) 'head : ',x2D(1:MIN(5,nx),1)
 WRITE(*,*) 'tail : ',x2D((nx-4):nx,ny)
 WRITE(*,*) 'n min mean max'
 WRITE(*,*) n, xmin, xbar, xmax

END SUBROUTINE summaryAI


INTEGER(I4B) ELEMENTAL FUNCTION ij2k (i,j,ncol)
 ! map 2D array [i.e. matrix] indexes into column-major array 1D index
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN) :: i,j,ncol

  ij2k = i + (j-1)*ncol
END FUNCTION ij2k

FUNCTION k2ij (k,ncol)
 ! map 2D array [i.e. matrix] indexes into column-major array 1D index
  IMPLICIT NONE

  INTEGER(I4B), INTENT(IN), DIMENSION(:)   :: k
  INTEGER(I4B), INTENT(IN)                 :: ncol
  INTEGER(I4B), DIMENSION(SIZE(k),2) :: k2ij

  ! local
  INTEGER(I4B), DIMENSION(SIZE(k)) :: i,j
  j = CEILING(REAL(k)/ncol)
  i = k - (j-1)*ncol
  k2ij(:,1) = i
  k2ij(:,2) = j
END FUNCTION k2ij

SUBROUTINE getLdnet(kinor, ldd, ldnet)
  ! Allocate and obtain a local drainage network object
  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:,:), INTENT(IN)    :: kinor        ! [ncell,2]
  INTEGER(I4B), DIMENSION(:,:), INTENT(IN)    :: ldd          ! [nx,ny]       7 8 9
  TYPE(T_ldnetwork), DIMENSION(:), POINTER :: ldnet        !               4 5 6
                                                              !               1 2 3
  ! local
  INTEGER(I4B) :: nx, ny, nkin, ikin, nup, i
  LOGICAL, DIMENSION(8) :: intome
  INTEGER(I4B), DIMENSION(8,2) :: inds

  nx = SIZE(ldd,1)
  ny = SIZE(ldd,2)
  nkin = SIZE(kinor,1)

  NULLIFY(ldnet)
  ALLOCATE(ldnet(nkin))

  !WRITE(*,*) 'nkin:',nkin
  DO ikin=1,nkin
    inds   = 0
    intome = .FALSE.
    inds(1,:) = (/ kinor(ikin,1)-1, kinor(ikin,2)-1 /)    ! ul
    inds(2,:) = (/ kinor(ikin,1)  , kinor(ikin,2)-1 /)    ! uc
    inds(3,:) = (/ kinor(ikin,1)+1, kinor(ikin,2)-1 /)    ! ur
    inds(4,:) = (/ kinor(ikin,1)-1, kinor(ikin,2)   /)    ! cl
    inds(5,:) = (/ kinor(ikin,1)+1, kinor(ikin,2)   /)    ! cr
    inds(6,:) = (/ kinor(ikin,1)-1, kinor(ikin,2)+1 /)   ! ll
    inds(7,:) = (/ kinor(ikin,1)  , kinor(ikin,2)+1 /)   ! lc
    inds(8,:) = (/ kinor(ikin,1)+1, kinor(ikin,2)+1 /)   ! lr

    !  WRITE(*,*) 'k =',k,'| inds(:,:)'
    !  DO i = 1, 8
    !    WRITE(*,*) inds(i,:)
    !  END DO

    IF (ldd(inds(1,1),inds(1,2)) == 3) intome(1) = .TRUE.              ! ul
    IF (ldd(inds(2,1),inds(2,2)) == 2) intome(2) = .TRUE.              ! uc
    IF (ldd(inds(3,1),inds(3,2)) == 1) intome(3) = .TRUE.              ! ur
    IF (ldd(inds(4,1),inds(4,2)) == 6) intome(4) = .TRUE.              ! cl
    IF (ldd(inds(5,1),inds(5,2)) == 4) intome(5) = .TRUE.              ! cr
    IF (ldd(inds(6,1),inds(6,2)) == 9) intome(6) = .TRUE.              ! ll
    IF (ldd(inds(7,1),inds(7,2)) == 8) intome(7) = .TRUE.              ! lc
    IF (ldd(inds(8,1),inds(8,2)) == 7) intome(8) = .TRUE.              ! lr

    nup = COUNT(intome)
    NULLIFY(ldnet(ikin)%up)
    ALLOCATE(ldnet(ikin)%up(nup,2))
    !WRITE(*,*) 'intome:',intome
    !WRITE(*,*) 'nup | SIZE(whichf(intome)) | whichf(intome):',nup,'|',SIZE(whichf(intome)), '|',whichf(intome)
    ldnet(ikin)%up = inds(whichf(intome),:)
    !WRITE(*,*) 'SIZE(ldnet(ikin)%up)',SIZE(ldnet(ikin)%up)   ! SIZE(ldnet%up) == 0 for 1st order pixel [most upstream]
    !WRITE(*,*) 'ldnet(ikin)%up'
    !DO i = 1, SIZE(ldnet(ikin)%up,1)
    !  WRITE(*,*) ldnet(ikin)%up(i,:)
    !END DO
 END DO

END SUBROUTINE getLDnet

SUBROUTINE getLDdidx(ldd,didx)
 IMPLICIT NONE
 ! get array of downstream (i,j) pixels for a LDD
 INTEGER(I4B), DIMENSION(:,:), INTENT(IN)    :: ldd
 INTEGER(I4B), DIMENSION(:,:,:), INTENT(OUT) :: didx

 ! local
 INTEGER(I4B) :: i,j,nx,ny

 nx = SIZE(ldd,1)
 ny = SIZE(ldd,2)

 DO j=1,ny
  didx(1,:,j) = (/ (i, i=1,nx) /)
  didx(2,:,j) = j
 END DO

 WHERE (ldd == 1)
   didx(1,:,:) = didx(1,:,:) - 1
   didx(2,:,:) = didx(2,:,:) + 1
 END WHERE
 WHERE (ldd == 2)
   didx(2,:,:) = didx(2,:,:) + 1
 END WHERE
 WHERE (ldd == 3)
   didx(1,:,:) = didx(1,:,:) + 1
   didx(2,:,:) = didx(2,:,:) + 1
 END WHERE
 WHERE (ldd == 4)
   didx(1,:,:) = didx(1,:,:) - 1
 END WHERE
 WHERE (ldd == 6)
   didx(1,:,:) = didx(1,:,:) + 1
 END WHERE
 WHERE (ldd == 7)
   didx(1,:,:) = didx(1,:,:) - 1
   didx(2,:,:) = didx(2,:,:) - 1
 END WHERE
 WHERE (ldd == 8)
   didx(2,:,:) = didx(2,:,:) - 1
 END WHERE
 WHERE (ldd == 9)
   didx(1,:,:) = didx(1,:,:) + 1
   didx(2,:,:) = didx(2,:,:) - 1
 END WHERE


 !WRITE(*,*) NEW_LINE('A'),'ldd(100:105,100:105)'
 !DO j=100,105
 !  WRITE(*,*) ldd(100:105,j)
 !END DO

 !WRITE(*,*) 'didx(100:105,100:105)'
 !WRITE(*,*) 'i:'
 !DO j=100,105
 !  WRITE(*,*) didx(1,100:105,j)
 !END DO
 !WRITE(*,*) 'j:'
 !DO j=100,105
 !  WRITE(*,*) didx(2,100:105,j)
 !END DO

END SUBROUTINE getLDdidx


SUBROUTINE getDdis(ldd, dx, ddis)
  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:,:), INTENT(IN)  :: ldd          ! [nx,ny]
  REAL(DP),                     INTENT(IN)  :: dx           ! [m] cell resolution
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: ddis         ! [nx,ny]

  ddis = 0.d0
  WHERE(ldd == 7 .OR. ldd == 9 .OR. ldd == 1 .OR. ldd == 3)
    ddis = dx * sqrt(2.d0)
  END WHERE
  WHERE(ldd == 8 .OR. ldd == 4 .OR. ldd == 6 .OR. ldd == 2 .OR. ldd == 5)
    ddis = dx
  END WHERE

END SUBROUTINE getDdis

SUBROUTINE getDslo0(dem,ldd,ddis,dsl,slo0)
  IMPLICIT NONE

  REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: dem         ! [nx,ny]
  INTEGER(I4B), DIMENSION(:,:), INTENT(IN)  :: ldd         ! [nx,ny]
  REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: ddis        ! [m] downstream distance
  REAL(DP),                     INTENT(IN)  :: dsl         ! downstream boundary condition [slope 0/1]
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: slo0        ! [nx,ny]
  !LOGICAL,      OPTIONAL,       INTENT(IN)  :: positive              

  ! local
  INTEGER(I4B) :: i, j, nx, ny

  nx = SIZE(dem,1)
  ny = SIZE(dem,2)

  DO j = 1,ny
     DO i = 1,nx
       slo0(i,j) = dem(i,j) - dem(didx(1,i,j),didx(2,i,j)) ! > 0 for downstream slope
     END DO
  END DO

  slo0 = slo0 / ddis

  WHERE (ldd == 5)
    slo0 = dsl
  END WHERE

  CALL AdjustSlope(slo0)

  !IF (PRESENT(positive)) THEN
    WHERE (slo0 < 1.0E-04)       ! 0.1 m/km min slope
      slo0 = 1.0E-04
    END WHERE  
  !END IF

END SUBROUTINE getDslo0

SUBROUTINE which(msk,x)
  IMPLICIT NONE
  LOGICAL, DIMENSION(:), INTENT(IN) :: msk
  INTEGER(I4B), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: x 
 
  IF (ALLOCATED(x)) THEN
    DEALLOCATE(x)
  END IF
  ALLOCATE(x(SIZE(whichf(msk))))
  x = whichf(msk) 
END SUBROUTINE which

FUNCTION whichf(msk)
  IMPLICIT NONE
  LOGICAL, DIMENSION(:), INTENT(IN) :: msk
  INTEGER(I4B), DIMENSION(MAX(1,COUNT(msk))) :: whichf

  !local
  INTEGER(I4B) :: i,n

  whichf = 0
  n = SIZE(msk)
  IF (n > 0) THEN
    whichf = PACK((/ (i,i=1,n) /),msk)
  END IF
END FUNCTION whichf

SUBROUTINE kin1D(t0, dti, mask, h)
  IMPLICIT NONE

  REAL(DP),                     INTENT(IN)    :: t0      ! real time [from model POW]
  REAL(DP),                     INTENT(IN)    :: dti     ! integration time
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask   ![nx,ny] matrix
  REAL(DP), DIMENSION(:,:),     INTENT(INOUT) :: h       ! [m] water depth

  ! local
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: Qb     ! [m3/s] Q^t     for the FD scheme
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: Qa     ! [m3/s] Q^{t+1} for the FD scheme
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: h_c    ! SGC water depth
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: eta    ! water surface level (z+h)
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: a      ! [-] alpha coefficient (a^t) in A=a*Q^b
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: ck     ! kinematic celerity (ck^t)
  REAL(DP),     DIMENSION(SIZE(h,1),SIZE(h,2)) :: P      ! wetted perimeter in SGC
  REAL(DP)                                     :: dx
  REAL(DP)                                     :: b      ! beta
  INTEGER(I4B)                                 :: i,j,k,t,nk,ik
  REAL(DP)                                     :: Qaup, Aaup
  REAL(DP)                                     :: Qbar         ! for linear Q estimate
  REAL(DP)                                     :: fQ,gQ,Ct     ! for non-linear Q estimate
  REAL(DP)                                     :: eps
  INTEGER(I4B)                                 :: STAT_CALL
  REAL(DP)                                     :: inthb,intha  ! integral of MASS pre/pos routing

    eps = 1.0E-10
    b   = 0.6D0
    dx  = G%nsres
    Qa  = 0.d0
    Qb  = 0.d0                                    
    P   = 0.D0
    nk  = SIZE(kinor,1)
                                                           ! WARNING: should alpha depend on flow direction?


  ! from initial conditions obtain alpha [a_b] and kinematic celerity
    inthb = SUM(h, MASK=mask) * dx**2
   ! WRITE(*,*) 'kin1D :: h mass pre-route [m3]:',inthb
   !    OPEN(UNIT=15, FILE='hpre.for', &
   !    STATUS='REPLACE', ACTION='WRITE', &
   !    FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
   !    WRITE(15) h
   !    CLOSE(15)

  CALL getHcPc(cmsk, cwid, cdep, G%nsres, h, h_c, P) ! map initial condition into SGC h and P

  !eta = dem - cdep + h_c                            does not seem adequate for big pixels
  eta = dem + h                                      ! smoother
  CALL getDslo0(eta,ldd,ddis,dsl,slo0)               ! difussive approximation but just downstream flow allowed

  !WRITE(*,*) NEW_LINE('A'),'summary(h_c, mask)'
  !CALL summaryAD(h_c, mask)
  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),'summary(P, mask)'
    CALL summaryAD(P, mask)
  END IF
  !     OPEN(UNIT=15, FILE='P.for', &
  !     STATUS='REPLACE', ACTION='WRITE', &
  !     FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
  !     WRITE(15) P
  !     CLOSE(15)


  a  = (nman * P **(2.0D0/3.0D0) / SQRT(slo0))**0.6D0 ! alpha: estimate for overland
  !WRITE(*,*) NEW_LINE('A'),'summary(a, mask)'
  !CALL summaryAD(a, mask)

  ck = sqrt(slo0) / nman * 5.0D0 / 3.0D0 * h_c**(2.0D0/3.0D0)    ! c_k approximation (9.3.15)
  !WRITE(*,*) NEW_LINE('A'),'summary(ck, mask)'
  !CALL summaryAD(ck, mask)
  !     OPEN(UNIT=15, FILE='ck.for', &
  !     STATUS='REPLACE', ACTION='WRITE', &
  !     FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
  !     WRITE(15) ck
  !     CLOSE(15)

  ! get Qb
  Qb = ((h * dx)/a)**(1.0D0/b) ! assume overland is concentrated in the channel for channel-cells

  IF (verbose > 0) THEN
    WRITE(*,*) 'dt:',dti,' | dx/ck:',G%nsres/MAXVAL(ck),' | ck:',MAXVAL(ck)
    WRITE(*,*) NEW_LINE('A'),'summary(Qb, mask)'
    CALL summaryAD(Qb, mask)
  END IF
  !    OPEN(UNIT=15, FILE='Qb.for', &
  !    STATUS='REPLACE', ACTION='WRITE', &
  !    FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
  !    WRITE(15) Qb
  !    CLOSE(15)

  ! get Qa
  DO k=1,nk

    i = kinor(k,1)
    j = kinor(k,2)
    Qaup = 0.0D0
    Aaup = 0.0D0
    IF (SIZE(ldnet(k)%up) > 0) THEN
      DO ik = 1,SIZE(ldnet(k)%up,1)
        Qaup = Qaup + Qa(ldnet(k)%up(ik,1),ldnet(k)%up(ik,2))   ! updated flow upstream
        Aaup = Aaup + dti / ddis(ldnet(k)%up(ik,1),ldnet(k)%up(ik,2)) * &
                      Qa(ldnet(k)%up(ik,1),ldnet(k)%up(ik,2))
      END DO
    END IF
    Qbar = (Qb(i,j) + Qaup) / 2.0D0
    Qa(i,j) = (Aaup + a(i,j) * b * Qbar**(b-1.0D0) * Qb(i,j)) / & ! (9.6.7) linear estimate
              (dti/ddis(i,j) + a(i,j) * b * Qbar**(b-1.0D0))

    Ct = Aaup + a(i,j) * Qb(i,j)**b
    fQ = 1000.0D0
    DO                           ! Newton iteration
      fQ = dti / ddis(i,j) * Qa(i,j) + a(i,j) * Qa(i,j)**b - Ct
      gQ = dti / ddis(i,j) + a(i,j) * b * Qa(i,j)**(b-1.0D0)
      Qa(i,j) = Qa(i,j) - fQ / gQ                          ! (9.6.17)
      IF (isnan(Qa(i,j))) THEN
        Qa(i,j) = 0.0D0
        fQ = 0.0D0
      END IF
      IF (ABS(fQ) <= eps) EXIT
    END DO
    IF (Qa(i,j) < 0.D0) Qa(i,j) = 0.D0
  END DO
  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),'summary(Qa, mask)'
    CALL summaryAD(Qa, mask)
  END IF
  !       OPEN(UNIT=15, FILE='Qa.for', &
  !     STATUS='REPLACE', ACTION='WRITE', &
  !     FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
  !     WRITE(15) Qa
  !     CLOSE(15)

  ! get updated alpha [a_a]


  ! get updated h [h_a]
  h = a * Qa** b / dx        ! assume alpha has not changed significantly
  !WRITE(*,*) 'summary(h, mask)'
  !CALL summaryAD(h, mask)

  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),' Qb(outlet):',Qb(ouij(1),ouij(2))
    WRITE(*,*) 'Qa(outlet):',Qa(ouij(1),ouij(2))
    WRITE(*,*) 'h(outlet):',h(ouij(1),ouij(2))
    intha = SUM(h, MASK=mask) * dx**2
    WRITE(*,*) 'kin1D :: h mass pre-route [m3]:', inthb
    WRITE(*,*) 'kin1D :: h mass pos-route [m3]:', intha
    WRITE(*,*) 'kin1D :: h mass loss [m3]:',inthb - intha
  END IF
  !OPEN(UNIT=15, FILE='hpos.for', &
  !     STATUS='REPLACE', ACTION='WRITE', &
  !     FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=STAT_CALL)
  !     WRITE(15) h
  !     CLOSE(15)


  !STOP 'kin1D Forced STOP'
END SUBROUTINE kin1D

ELEMENTAL SUBROUTINE AdjustSlope (Slope)
  IMPLICIT NONE

  REAL(DP), INTENT(INOUT)      :: Slope

  ! local
  REAL(DP)                     :: sign

        !Slope correction given by City of Albuquerque, 1997, p.22-26
        !http://www.hkh-friend.net.np/rhdc/training/lectures/HEGGEN/Tc_3.pdf


        if (Slope .LT. 0.0D0) then
            sign = -1.0
        else
            sign = 1.0
        end if

        Slope = abs (Slope)

        if (Slope .GE. 0.04) then
            Slope = 0.05247 + 0.06363 * Slope - 0.182 * exp (-62.38 * Slope)
        end if
        if (Slope .LE. 1.0E-6) then
            Slope = 1.0E-06
        end if

        Slope = sign * Slope

END SUBROUTINE AdjustSlope

SUBROUTINE getHcPc(cmsk,cwid,cdep,dx, h,h_c,P)
  ! get SGC water depth and wetted Perimeter
  IMPLICIT NONE

  INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: cmsk
  REAL(DP),     DIMENSION(:,:), INTENT(IN) :: cwid
  REAL(DP),     DIMENSION(:,:), INTENT(IN) :: cdep
  REAL(DP),                     INTENT(IN) :: dx
  REAL(DP),     DIMENSION(:,:), INTENT(IN) :: h
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: h_c  ! SGC : water depth from channel bottom to water surface [within or out of bank]
  REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: P

  !local

  WHERE (cmsk == 0)
     h_c = h
     P   = dx
  ELSEWHERE
     WHERE (cwid * cdep > dx * h) ! within bank
        h_c = dx * h / cwid
        P   = cwid + 2*h_c
     ELSEWHERE                    ! overbank flow
        h_c = cdep + (dx * h - cwid * cdep) / dx
        P   = dx + 2 * cdep
     END WHERE
  END WHERE

END SUBROUTINE getHcPc


END MODULE ModuleFloodTools
