PROGRAM hydroLand
! Author: Javier Garcia-Pintado. Last version: 08/02/2013
! Purpose:
!
  ! Compilation. Option a) static linking (options -g -O produce debugging information and default optimizization)
  ! note the -cpp 'preprocessor directive' option is automatically included for .F90 files (not for .f90 , which explicitly need '-ccp')

  ! gfortran -g -O -c  nrtype.F90 -o nrtype.o
  ! gfortran -g -O -c  ModuleGlobalData.F90 -o ModuleGlobalData.o
  ! gfortran -g -O -c  ModuleTime.F90 -o ModuleTime.o
  ! gfortran -g -O -c  ModuleDrawing.F90 -o ModuleDrawing.o
  ! gfortran -g -O -c  ModuleTimeSerie.F90 -o ModuleTimeSerie.o
  ! gfortran -g -O -c  ModuleEnterData.F90 -o ModuleEnterData.o
  ! gfortran -g -O -c  ModuleStopWatch.F90 -o ModuleStopWatch.o
  ! gfortran -g -O -c  ModuleFunctions.F90 -o ModuleFunctions.o
  ! gfortran -g -O -c  ModuleJGPFunctions.F90 -o ModuleJGPFunctions.o
  ! gfortran -g -O -c  ModuleFillMatrix.F90 -o ModuleFillMatrix.o
  ! gfortran -g -O -c  ModuleHorizontalGrid.F90 -o ModuleHorizontalGrid.o
  ! gfortran -g -O -c  ModuleHorizontalMap.F90 -o ModuleHorizontalMap.o
  ! gfortran -g -O -c  ModuleGridData.F90 -o ModuleGridData.o
  ! gfortran -g -O -c  ModuleBasinGeometry.F90 -o ModuleBasinGeometry.o
  ! gfortran -g -O -c  ModuleRunOff.F90 -o ModuleRunOff.o
  ! gfortran -g -O -c  ModuleRunOffProperties.F90 -o ModuleRunOffProperties.o
  ! gfortran -g -O -c  ModuleDrainageNetwork.F90 -o ModuleDrainageNetwork.o
  ! gfortran -g -O -c  ModulePorousMedia.F90 -o ModulePorousMedia.o
  ! gfortran -g -O -c  ModulePorousMediaProperties.F90 -o ModulePorousMediaProperties.o
  ! gfortran -g -O -c  ModuleGeometry.F90 -o ModuleGeometry.o
  ! gfortran -g -O -c  ModuleHydroPoint.F90 -o ModuleHydroPoint.o
  ! gfortran -g -O -c  ModuleHspfPwater.F90 -o ModuleHspfPwater.o

  ! gfortran -g -O -c  ModuleBasin.F90 -o ModuleBasin.o
  ! gfortran -g -O -c  hydroLand.F90 -o hydroLand.o
  ! gfortran -o hydroLand hydroLand.o nrtype.o ModuleGlobalData.o ModuleTime.o ModuleEnterData.o ModuleHydroPoint.o ModuleHspfPwater.o ModuleStopWatch.o ModuleFunctions.o

  !
! Compilation. Option b) build module as dynamic library and just link hydroLand against it
!  notes:         -Wall :: enables commonly used warning options
!
! gfortran  -fpic -g -O -c  -Wall nrtype.f90 -o nrtype.o
! gfortran  -fpic -g -O -c  -Wall hydrosim_hspf_mod.f90 -o hydrosim_hspf_mod.o
!
!  now, online manuals say to do:
! gfortran  -shared -Wl,-soname,libdisthspf.so.1 -o libdisthspf.so.1.0.0 nrtype.o hydrosim_hspf_mod.o
!
!  and to install the shared library, copy it to one standard directory (e.g., /home/pt902904/lib),
!  and create a soft link named as a stripped version "*.so" and run ldconfig on the library folder; e.g.,
! ldconfig -n /home/pt902904/lib
!
!  but all this is not working later when calling the programm. So, having /home/pt902904/lib in my $LD_LIBRARY_PATH,
!  I am just doing
! gfortran  -shared -o libdisthspf.so nrtype.o hydrosim_hspf_mod.o
! mv libdisthspf.so ~/lib
!
!  Finally, compile against the dynamic library
! gfortran -L/home/pt902904/lib -o runhydro runhydro.f90 -ldisthspf
! mv runhydro ~/bin

    USE ModuleGlobalData
    USE ModuleEnterData
    USE ModuleTime
    USE ModuleFunctions,      only : ReadTimeKeyWords
    !USE ModuleBasin
    USE ModuleHydroPoint
    USE ModuleStopWatch,      only : CreateWatchGroup, KillWatchGroup

    IMPLICIT NONE

    !Parameters

    !Instance IDs
    INTEGER                             :: ObjComputeTime       = 0
    INTEGER                             :: ObjBasin             = 0
    INTEGER                             :: ObjHydropts          = 0

    !Time Variables
    TYPE (T_Time)                       :: BeginTime, EndTime, CurrentTime
    REAL                                :: DT                   = null_real
    REAL                                :: MaxDT                = null_real
    LOGICAL                             :: VariableDT           = .FALSE.
    REAL                                :: GmtReference         = null_real
    REAL                                :: DTPredictionInterval = null_real

    CHARACTER(len=StringLength)  :: modelname   = 'HSPF_PWATER'   ! model selected to run. Possible values: 'HSPF_PWATER'
    CHARACTER(len=PathLength)    :: dsnfile     = null_str        ! main input data file name, unique argument to run hydroland (just contains 'FilesName' and 'dsn')
    LOGICAL                      :: SaveOutput  = .TRUE.          ! TODO: init to .FALSE. after debugging
    CHARACTER(len=PathLength)    :: cwd         = null_str        !

    !Other Stuff
    TYPE (T_Time)                       :: InitialSystemTime, FinalSystemTime
    INTEGER, DIMENSION(8)               :: F95Time              = null_int

    LOGICAL                             :: ModelConstructed     = .FALSE.

    INTEGER       :: status   = 0                                                    ! O
    INTEGER       :: verbolev = 1                                                    ! standard output information level

    CALL get_command_argument(1, dsnfile)                                            ! intrinsic
    IF (LEN_TRIM(dsnfile) == 0) THEN
      STOP 'runHydro - input file required - ERR00'
    ELSE
      dsnfile = trim(adjustl(dsnfile))
    END IF

    OPEN( UNIT=11, FILE=dsnfile, STATUS='OLD', ACTION='READ', &
          IOSTAT=status)
    IF ( status == 0 ) THEN
       READ(11,'(A)') FilesName                                                      ! overwrites FilesName default set in ModuleGlobalData
       READ(11,'(A)') dsn                                                            !    "       dsn='.'   default set in ModuleGlobalData
       CLOSE(11)
       FilesName = trim(adjustl(FilesName))
       dsn       = trim(adjustl(dsn))
       !datafile0 = dsn(1:LEN_TRIM(dsn)) // "/" // datafile0(1:LEN_TRIM(datafile0))
    ELSE
       WRITE (*,1002) status
       1002 FORMAT (1X,'open runfile failed--status = ', I6)
       STOP
    END IF
    IF (verbolev > 0) THEN
         WRITE (*,*) 'model name:   ', modelname
         WRITE (*,*) 'FilesName:    ', FilesName
         WRITE (*,*) 'dsn:          ', dsn
    END IF
    CALL getcwd(cwd, status)
    IF (status .NE. SUCCESS_) STOP 'Hydroland - ERR01'
    CALL chdir(dsn, status)
    IF (status .NE. SUCCESS_) STOP 'Hydroland - ERR02'

    CALL ConstructHydroLand ! TODO
    ! CALL ModifyMohidLand    TODO
    CALL KillHydroLand      ! TODO


SELECT CASE (modelname)
CASE('HSPF_PWATER')
 IF (verbolev > 0) THEN
   WRITE (*,'(1X,A15,A8,A1/)') 'Running model: ',modelname ,'!'
 END IF
 !CALL HSPF_Pwater(DataFile)
CASE DEFAULT
   WRITE (*,'(A,A)') 'Model not implemented: ', modelname
END SELECT

    CONTAINS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONSTRUCTOR CONS

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    SUBROUTINE ConstructHydroLand

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
      INTEGER                                       :: STAT_CALL

      CALL StartupHydro(modelname)                                                               ! in ModuleGlobalData :: just open out files for error and usedkeywords

      !Gets the actual machine time
      call date_and_time(Values = F95Time)                                                       ! intrinsic procedure
      call SetDate      (InitialSystemTime, float(F95Time(1)), float(F95Time(2)),      &         ! in ModuleTime
                                            float(F95Time(3)), float(F95Time(5)),      &
                                            float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)

        !Reads the Keywords    TODO
        call ReadKeywords                                                                        ! in this program : CONSTRUCTOR

        !Constructs Time     TODO
        call StartComputeTime (ObjComputeTime, InitialSystemTime, BeginTime, EndTime, DT, VariableDT, MaxDT, STAT = STAT_CALL)

        !Update Current Time
        CurrentTime  = BeginTime

        !Constructs Basin TODO
        !call ConstructBasin   (ObjBasinID = ObjBasin, ObjTime = ObjComputeTime, ModelName = ModelName, STAT = STAT_CALL)

        !Construct Hydropts: points where timeseries of variables are produced as output
        CALL ConstructHydropts (ObjHydroptsID = ObjHydropts, ObjTime = ObjComputeTime, ModelName = ModelName, STAT = STAT_CALL)

        ModelConstructed = .TRUE.


    END SUBROUTINE ConstructHydroLand

    !--------------------------------------------------------------------------

    SUBROUTINE ReadKeywords

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        character(PathLength)                       :: DataFile, TimeFile, HydroptsFile
        integer                                     :: STAT_CALL
        integer                                     :: ObjEnterData = 0
        integer                                     :: iflag
        character(PathLength)                       :: WatchFile, OutputFile, DTLogFile

        !Monitor Performance of the model execution?
        call ReadFileName('OUTWATCH', WatchFile, Message = 'Start Watch File', STAT = STAT_CALL) ! Reads input from FilesName at path dsn
        if (STAT_CALL == SUCCESS_) then
            STAT_CALL = CreateWatchGroup (WatchFile)
            MonitorPerformance  = .true.
        else
            MonitorPerformance  = .false.
        endif


        !Monitor Performance of the model execution?
        call ReadFileName('OUTPUT_RESULT', OutputFile, Message = 'Start Result File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            SaveOutput  = .true.
        else
            SaveOutput  = .false.
        endif

        call ReadFileName('DT_LOG', DTLogFile, Message = 'Start DTLog File', STAT = STAT_CALL)
        if (STAT_CALL == SUCCESS_) then
            MonitorDT  = .true.
        else
            MonitorDT  = .false.
        endif

        if (MonitorDT) then
            call UnitsManager (UnitDT, OPEN_FILE)
            open(UNIT   = UnitDT, FILE   = DTLogFile, STATUS  = "UNKNOWN", IOSTAT  = STAT_CALL)
            if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR01'
            write (UnitDT, '(A25, A10,A12,A12,A13,A13,A13,A13,A26)') &
                "ModuleName", "iter", "DT", "DNet", "RunOff", "PorousMedia", "Atmosphere", "DTNextEv", "NextTime"
        end if

        call ReadFileName('IN_MODEL', DataFile, "HydroLand Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR02'

        call ReadFileName('IN_TIME', TimeFile, "HydroLand Time Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR03'

        call ReadFileName('IN_HYDROPTS', HydroptsFile, "HydroLand Hydropts Data File", STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR04'


        call ConstructEnterData (ObjEnterData, DataFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR05'

        !Model Name
        call GetData(ModelName,                                                          &
                     ObjEnterData, iflag,                                                &
                     SearchType   = FromFile,                                            &
                     keyword      = 'MODEL_NAME',                                        &
                     default      = 'HYDRO Land Model',                                  &
                     ClientModule = 'HYDROLand',                                         &
                     STAT         = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR06'

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR07'

        call ConstructEnterData (ObjEnterData, TimeFile, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR08'

        call ReadTimeKeyWords   (ObjEnterData, FromFile, BeginTime, EndTime, DT,         &
                                 VariableDT, "Hydro Land", MaxDT, GmtReference,          &
                                 DTPredictionInterval)

        call KillEnterData (ObjEnterData, STAT = STAT_CALL)
        if (STAT_CALL /= SUCCESS_) stop 'ReadKeywords - HydroLand - ERR09'




    end subroutine ReadKeywords

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODIFIER MODI

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR DESTRUCTOR

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    SUBROUTINE KillHydroLand

        !Arguments-------------------------------------------------------------

        !Local-----------------------------------------------------------------
        INTEGER                                     :: STAT_CALL
        REAL                                        :: ElapsedSeconds, TotalCPUTime

        !CALL KillBasin(ObjBasinID = ObjBasin, STAT = STAT_CALL)

        !IF (MonitorPerformance) then
        !    call KillWatchGroup (STAT_CALL)
        !    if (STAT_CALL /= SUCCESS_) stop 'KillMohidLand - MohidLand - ERR01'
        !endif

        !if (MonitorDT) call UnitsManager (UnitDT, CLOSE_FILE)

        CALL date_and_time(Values = F95Time)
        CALL SetDate      (FinalSystemTime,   float(F95Time(1)), float(F95Time(2)),      &
                                              float(F95Time(3)), float(F95Time(5)),      &
                                              float(F95Time(6)), float(F95Time(7))+      &
                                              F95Time(8)/1000.)
        call cpu_time(TotalCPUTime)                                                               ! intrinsic, output (REAL) is current CPU time
        ElapsedSeconds = FinalSystemTime - InitialSystemTime

        !if (SaveOutput) call SaveRunInfo (modelname, ElapsedSeconds, TotalCPUTime, OutputFile) ! in ModuleGlobalData

        call ShutdownHydro (modelname, ElapsedSeconds, TotalCPUTime)

        ModelConstructed = .false.


    END SUBROUTINE KillHydroLand

END PROGRAM hydroLand

!CHARACTER(len=50)  :: string1
!CHARACTER(len=50)  :: string2
!REAL(DP), DIMENSION(3) :: arr1
!REAL(DP), DIMENSION(3) :: arrexp1
!arr1 = [2,3,5]
!arrexp1 = [2,1,3]
!WRITE(*,*) 'arr1**arrexp1', arr1**arrexp1
!REAL(DP), DIMENSION(4) :: arrayfoo
!LOGICAL(LGT), DIMENSION(4) :: arrayfoo2
!string1 = 'esto es una cadena '
!string2 = 'de texto'
!string2 = string1 // string2
!WRITE(*,*) string2
