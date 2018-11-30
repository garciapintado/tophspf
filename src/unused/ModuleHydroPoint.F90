! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! ::::: Parameters and variables related to gauges
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
MODULE ModuleHydroPoint

    USE ModuleGlobalData

    IMPLICIT NONE

    PRIVATE

    !Subroutines--------------------------------------------------------------

    !Constructor
    PUBLIC  :: ConstructHydropts
    !PRIVATE ::   AllocateInstance
    !PRIVATE ::   ReadFileNames
    !PRIVATE ::   ReadDataFile
    PUBLIC  :: set_gauges

    !Management
    PRIVATE :: Ready

    !Type----------------
    PUBLIC :: T_HydroPoint
    TYPE T_HydroPoint
      REAL,                      DIMENSION(2)     :: coords     = null_real
      CHARACTER(len=PathLength)                   :: name       = null_str
      INTEGER                                     :: ihru       = null_int
      LOGICAL                                     :: isou       = .FALSE.       ! outlet for this ihru catchment
    END TYPE T_HydroPoint

    PUBLIC :: T_HydroPoints
    TYPE      T_HydroPoints
       INTEGER, DIMENSION(:), POINTER :: data
    END TYPE T_HydroPoints




    !Global Module Variables
    TYPE (T_HydroPoint), DIMENSION(:), ALLOCATABLE   :: hydropts
    !TYPE (T_HydroPoint), POINTER                    :: FirstObjHydroPoint => null()
    !TYPE (T_HydroPoint), POINTER                    :: Me                 => null()



    !integer, parameter :: OUTGAUGEUNIT=89
    !integer :: num_gauges
    !real(kind=8), allocatable :: xgauge(:), ygauge(:), t1gauge(:), t2gauge(:)
    !integer, allocatable ::  mbestsrc(:), mbestorder(:), igauge(:)

CONTAINS

  SUBROUTINE ConstructHydropts(ObjHydroptsID, ObjTime, ModelName, STAT)

        !Arguments---------------------------------------------------------------
        INTEGER, INTENT(INOUT)                          :: ObjHydroptsID
        integer                                         :: ObjTime
        character(len=*)                                :: ModelName
        integer, intent(OUT), OPTIONAL                  :: STAT

        STAT_ = UNKNOWN_


  END SUBROUTINE ConstructHydropts

    subroutine set_gauges(fname)

!        use amr_module

        IMPLICIT NONE

        ! Input
        character(len=*), intent(inout), optional :: fname

        ! Locals
        integer :: i
        integer, parameter :: iunit = 7

        ! Open file
        IF (.NOT. PRESENT(fname)) fname = 'hydroloc.dat'


!            call opendatafile(iunit,fname)
!        else
!            call opendatafile(iunit,'gauges.data')
!        endif

!        read(iunit,*) num_gauges

!        allocate(xgauge(num_gauges), ygauge(num_gauges))
!        allocate(t1gauge(num_gauges), t2gauge(num_gauges))
!        allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
!        allocate(igauge(num_gauges))

!        do i=1,num_gauges
!            read(iunit,*) igauge(i),xgauge(i),ygauge(i),t1gauge(i),t2gauge(i)
!        enddo

!        close(iunit)

        ! initialize for starters
!        mbestsrc = 0

        ! open file for output of gauge data all data is output in one binary
        ! file with format gauge number, level, time, depth by dumpgauge.
!        open(unit=OUTGAUGEUNIT, file='fort.gauge', status='unknown', &
!                                form='formatted')

    end subroutine set_gauges

end module ModuleHydroPoint
