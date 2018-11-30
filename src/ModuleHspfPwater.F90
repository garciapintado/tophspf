MODULE ModuleHspfPwater
  ! programming conventions
  ! - Temporary arrays in subroutines are created as automatic arrays. Note these cannot be SAVEd.
  ! - Dummy arrays in subroutines are declared as assumed shape dummy arrays (see Chapman, section 6.5).

  USE nrtype
  USE iso_varying_string
  USE strings, only : split
  USE ModuleFloodTools, only : ConstructKin1D, DestructKin1D, kin1D            ! subroutines
  USE ModuleFloodTools, only : ReadTimes2IRTS, WriteIRTS, which                ! subroutines
  USE ModuleFloodTools, only : whichf                                          ! function
  USE ModuleFloodTools, only : T_irts                                          ! data types

  !USE ModuleGlobalData
  !USE ModuleTime
  !USE ModuleRunOff,           only : KineDiffWave

  ! modules for distributed modelling - shallow water equations
  USE topo_module,       ONLY : mxtopo, mytopo     ! just for array DIM checks
  USE amr2_module,       ONLY : amr2               ! distributed routing

!#ifdef _ENABLE_CUDA
!  USE ModuleCuda
!#endif _ENABLE_CUDA

    IMPLICIT NONE

    PRIVATE

    !Subroutines
    !Constructor
    PUBLIC  :: HSPF_Pwater
    PRIVATE ::   read_geometa
    PRIVATE ::   allo_2Darr
    PRIVATE ::   allo_2Darr_spdis
    PRIVATE ::   allo_tsarr_hrus
    PRIVATE ::   deallo_common
    PRIVATE ::   deallo_2Darr
    PRIVATE ::   deallo_2Darr_spdis
    PRIVATE ::   deallo_tsarr_hrus
    PRIVATE ::   read_static_dbl
    PRIVATE ::   read_static_int
    PRIVATE ::   read_logicalASC
    PRIVATE ::   read_intASC
    PRIVATE ::   integrate_to_hru
    PRIVATE ::   read_rain
    PRIVATE ::   read_petinp
    PRIVATE ::   icept
    PRIVATE ::   surfac
    PRIVATE ::   divisn
    PRIVATE ::   dispos
    PRIVATE ::   etbase
    PRIVATE ::   evicep
    PRIVATE ::   etuzon
    PRIVATE ::   etuzs
    PRIVATE ::   etagw
    PRIVATE ::   etlzon
    PRIVATE ::   intflw
    PRIVATE ::   uzone
    PRIVATE ::   uzones
    PRIVATE ::   lzone
    PRIVATE ::   gwater
    PRIVATE ::   proute

    PRIVATE ::   write_ts
    PRIVATE ::   drouteK01
    PRIVATE ::   drouteSWE

    !Types
    !TYPE T_Coupling
    !    LOGICAL                                     :: RunOff               = .FALSE.
    !END TYPE T_Coupling

    !TYPE T_ExtVar
    !    INTEGER, DIMENSION(:),   POINTER            :: HydroPoints          => null()
    !    REAL,    DIMENSION(:,:), POINTER            :: Topography           => null()
    !END TYPE T_ExtVar

    !External variables from Runoff but that will be updated and sent to Runoff
    !TYPE T_ExtUpdate
    !    REAL(DP), dimension(:,:), pointer            :: WaterLevel             => null()
    !    REAL(DP), dimension(:,:), pointer            :: WaterColumn            => null()
    !    REAL(DP), dimension(:,:), pointer            :: WaterColumnOld         => null()
    !END TYPE T_ExtUpdate

    !TYPE T_Files
    !    character(len=PathLength)                   :: TopographicFile        = null_str
    !END TYPE T_Files

    ! Global Module Variables

    ! declare input Variables
!-----------------------
REAL(DP)             :: dx                          !required: geometa data
INTEGER(I4B)         :: ncols, nrows                !required: geometa data [cols = eastings [nx], rows = northings [ny]
REAL(DP)             :: dto, dti                    !required: geometa data
INTEGER(I4B)         :: nto                         !required: geometa data
CHARACTER(len=5)     :: simspdist                   !required: geometa data flag: 'spdis' for spatially distributed simulation, 'smdis' for semidistributed (clasical HSPF) simulation
CHARACTER(len=3)     :: rwave                       !required: geometa data flag: 'k01' [1D kinematic],'k02' [2D kinematic], or 'swe' [shallow water equations]
!INTEGER(I4B)         :: ntob                       !number of base outer timesteps (those with meteorological input; ntob >= nto)
INTEGER(I4B)         :: nx, ny                      ! nx <- cols | ny <- rows

REAL(DP)              :: t0 = 0.0d0                 !model-space time at the start of each interval
TYPE(T_irts), POINTER :: otimes => null()           !model-space time: output times 
INTEGER(I4B), ALLOCATABLE, DIMENSION(:) :: aid     ! general array index vector for sub-array operations

LOGICAL     , ALLOCATABLE, DIMENSION(:,:) :: mask   !required: logical 2D mask .TRUE. for watershed
INTEGER(I4B), ALLOCATABLE, DIMENSION(:,:) :: hrus   !required: 1...nhrus for subwatersheds, 0 for values outside watershed

! PWAT-PARM1
INTEGER(I4B) :: csnofg = 0
INTEGER(I4B) :: rtopfg = 0                          !either 0 (new method) or 1 (clasical) TODO: implement 2 and 3 (both are for hwtfg = 1) rtopfg set to 5 (kinematic) or 6 (diffusion wave) for 'spdis'
INTEGER(I4B) :: uzfg   = 1                          !TODO: implement 0
INTEGER(I4B) :: vcsfg  = 0                          !for cespc. TODO: for v... flags, implement possibly time-varying parameters.
INTEGER(I4B) :: vuzfg  = 0                          ! "  uzsn
INTEGER(I4B) :: vnnfg  = 0                          ! "  nsur
INTEGER(I4B) :: vifwfg = 0                          ! "  intfw
INTEGER(I4B) :: vircfg = 0                          ! "  irc
INTEGER(I4B) :: vlefg  = 0                          ! "  lzetp
INTEGER(I4B) :: iffcfg = 0
INTEGER(I4B) ::  hwtfg = 0
INTEGER(I4B) :: irrgfg = 0
! END PWAT-PARM1
! PWAT-PARM2
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzsn              ! [mm]    lower zone nominal storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infilt            ! [mm/hr] index to the infiltration capacity of the soil
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lsur              ! [m]     length of assumed overland flow plane
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: slsur             ! [0/1]   slope of the overland flow plane
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: kvary             ! [1/mm]  parameter
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwrc             ! [1/day]
! END PWAT-PARM2
! PWAT-PARM3
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infexp            ! [-]     exponent of the infiltration equation                     PWAT-PARM3
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infild            ! [-]     ratio between the maximum and mean infiltration capacities over the PLS                   "
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: deepfr            ! [-]     fraction of groundwater inflow which will enter deep (inactive) groundwater (lost from the system)
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: basetp            ! [-]     fraction of remaining potential ET which can be satisfied from baseflow, if enough is available
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwetp            ! [-]     fraction of remaining potential ET which can be satisfied from active groundwater if available
! END PWAT-PARM3
! PWAT-PARM4
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepsc             ! [mm] maximum interception capacity
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzsn              ! [mm] upper soil zone nominal storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: nsur              ! [cplx]  manning coefficient for overland flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: intfw             ! [-]     interflow inflow parameter
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: irc               ! [1/day] interflow recession parameter (ratio of today's interflow outflow rate to yesterday's rate)
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzetp             ! [-]     lower zone ET parameter. It is an index to the density of deep-rooted vegetation
! END PWAT-PARM4

! TIME SERIES INPUT
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rain       !external: rainfall                      for timestep
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: petinp     !external: potential evapotranspiration  for timestep
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: surli      !external: lateral inflow to surface     for timestep
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzli      !external: lateral inflow to upper zone  for timestep
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwli      !external: lateral inflow to upper zone  for timestep
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzli      !external: lateral inflow to upper zone  for timestep
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwli      !external: lateral inflow to upper zone  for timestep

! INITIALLY COMPUTED VARIABLES
INTEGER(I4B)         :: nti                    !internal number of timesteps. It will be > nto for dti < dto
REAL(DP)             :: stepdiv                     !step divisor.    external flow rates will be evenly divided by this amount for each internal timestep
INTEGER(I4B)         :: nhrus                       !number of hydrological response units
REAL(DP)             :: delt60                      !simulation time interval in hours

! TIME SERIES COMPUTED BY PERLND (allocated just for simspdist='spdis')
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pers       !sum of perlnd storages                 [mm] storage  pers = ceps + surs + ifws + uzs + lzs + agws
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ceps       !actual interception storage            [mm] storage
REAL(DP), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: surs       !actual surface detention storage       [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifws       !actual interflow storage               [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzs        !actual upper zone storage              [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzs        !actual lower zone storage              [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agws       !actual active groundwater storage      [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: gwvs       !groundwater variable storage index     [mm] storage
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: petadj     !adjustment factor for potential ET     [-]  storage !csnofg == 0

REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: supy       !rainfall+snow                          [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: suro       !surface outflow                        [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwo       !interflow outflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwo       !active groundwater outflow             [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pero       !total outflow from PSL                 [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: igwi       !inflow to inactive (deep) GW           [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pet        !potential ET, adjusted for snow and T  [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepe       !actual evapot. from ceps               [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzet       !ET from upper zone                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzet       !ET from lower zone                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwet      !ET from active groundwater storage     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: baset      !ET from active groundwater outflow     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: taet       !total simulated ET                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwi       !interflow inflow (excluded lateral)    [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzi        !upper zone inflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infil      !infiltration to the soil               [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: perc       !percolation from upper to lower zone   [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzi        !lower zone inflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwi       !active groundwater inflow (excl. lat.) [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: suri       !surface inflow (included lateral)      [mm/intv] flow

! TIME SERIES COMPUTED BY PERLND NOT INCLUDED IN THE CATALOG (allocated just for simspdist='spdis')
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepo       !output from ceps                       [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rempet     !remaining potential evapotranspiration [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: iperc      !percolation from upper to lower zone   [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: watin      !total input of water to the pervious land segment [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: watdif     !net input of water to the pervious land segment [mm/intv] flow

!arrays that overwrite themselves
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: dec        !routing parameter
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: src        !routing parameter

!DEFINE SEMIDISTRIBUTED VARIABLES
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: mask_hru

! PWAT-PARM2
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzsn_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infilt_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lsur_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: slsur_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: kvary_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwrc_hru
! END PWAT-PARM2
! PWAT-PARM3
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infexp_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infild_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: deepfr_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: basetp_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwetp_hru
! END PWAT-PARM3
! PWAT-PARM4
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepsc_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzsn_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: nsur_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: intfw_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: irc_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzetp_hru
! END PWAT-PARM4

! TIME SERIES INPUT (always allocated)
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rain_hru        ! arrays (nti,nhrus+1) to store hru-lumped time series
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: petinp_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: surli_hru
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzli_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwli_hru
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzli_hru
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwli_hru

! TIME SERIES COMPUTED BY PERLND (always allocated)
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pers_hru       !sum of perlnd storages                 [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ceps_hru       !actual interception storage            [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: surs_hru       !actual surface detention storage       [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifws_hru       !actual interflow storage               [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzs_hru        !actual upper zone storage              [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzs_hru        !actual lower zone storage              [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agws_hru       !actual active groundwater storage      [mm] storage
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: gwvs_hru       !groundwater variable storage index     [mm] storage
!REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: petadj_hru     !adjustment factor for potential ET    [-]  storage !csnofg == 0

REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: supy_hru       !rainfall+snow                          [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: suro_hru       !surface outflow                        [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwo_hru       !interflow outflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwo_hru       !active groundwater outflow             [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pero_hru       !total outflow from PSL                 [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: igwi_hru       !inflow to inactive (deep) GW           [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: pet_hru        !potential ET, adjunted for snow and Tª [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepe_hru       !actual evapot. from ceps               [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzet_hru       !ET from upper zone                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzet_hru       !ET from lower zone                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwet_hru      !ET from active groundwater storage     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: baset_hru      !ET from active groundwater outflow     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: taet_hru       !total simulated ET                     [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwi_hru       !interflow inflow (excluded lateral)    [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: uzi_hru        !upper zone inflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: infil_hru      !infiltration to the soil               [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: perc_hru       !percolation from upper to lower zone   [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: lzi_hru        !lower zone inflow                      [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: agwi_hru       !active groundwater inflow (excl. lat.) [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: suri_hru       !surface inflow (included lateral)      [mm/intv] flow

! TIME SERIES COMPUTED BY PERLND NOT INCLUDED IN THE CATALOG
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: cepo_hru       !output from ceps                       [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: rempet_hru     !remaining potential evapotranspiration [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: iperc_hru      !percolation from upper to lower zone   [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: watin_hru      !total input of water to the pervious land segment [mm/intv] flow
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: watdif_hru     !net input of water to the pervious land segment [mm/intv] flow

!arrays that overwrite themselves
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: dec_hru        !routing parameter
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: src_hru        !routing parameter

!arrays adaptable to spdis/smdis analysis
LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: maskloc
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: ifwk1, ifwk2, lzrat, rlzrat, lzfrac, gwi, rparm   !overwrite themselves
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: kgw                                               !calculated parameter
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: msupy      !msupy is used in PWATER and again in subroutine group surfac
REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: surss      !surss is the value of surs at the start of the ivl - it is used by module section MSTLAY
INTEGER(I4B), ALLOCATABLE, DIMENSION(:,:) :: smsfg  !flag state variable
INTEGER(I4B), ALLOCATABLE, DIMENSION(:,:) :: fsmsfg !flag state variable


! +++ Global flags and parameters incorporated just to make code closer to HSPF +++
INTEGER(I4B) :: slifp = 0                                  !flag for inflow timeseries| 0 for no lateral inflow, > 0 for lateral inflow
INTEGER(I4B) :: ilifp = 0                                  !flag for interflow lateral inflow timeseries| 0 for no lateral inflow, > 0 for lateral inflow
INTEGER(I4B) :: alifp = 0                                  !flag for groundwater lateral inflow timeseries| 0 for no lateral inflow, > 0 for lateral inflow
REAL(DP)     :: inffac = 1.0_dp                            !factor to reduce infiltration and percolation to account for frozen ground
INTEGER(I4B) :: verbolev = 0                               !control the amount of output information during simulation

CHARACTER(len=PathLength) :: dsno
CHARACTER(len=10) :: t0str, tostr                          ! actual interval time, and closest output time

CONTAINS

SUBROUTINE HSPF_Pwater (modelfilep)
!
! Purpose:
!  distributed version on hydrological part of HSPF model 12.0
!
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: modelfilep             !pathed main input file to topHSPF
 INTEGER(I4B) :: error
 INTEGER(I4B) :: status
 CHARACTER(len=175) :: fgeometap     !file with global data
 CHARACTER(len=175) :: fmaskp        !input file for LOGICAL mask .TRUE. for watershed
 CHARACTER(len=PathLength) :: fhrusp !input file for INTEGER mask 1...nhrus for watershed, 0 outside watershed
 CHARACTER(len=PathLength) :: fk01p  !main input to kinematic wave
 CHARACTER(len=175) :: frainp        !input file for rainfall
 CHARACTER(len=175) :: rain_std      !input type for rainfall
 !TYPE(varying_string) :: rain_std       !input file for rainfall
 CHARACTER(len=175) :: fpetinpp      !input file for potential evapotranspiration
 CHARACTER(len=175)   :: petinp_std  !input type for potential evapotranspiration

 CHARACTER(len=175) :: flzsnp       !PWAT-PARM2 filenames
 CHARACTER(len=175) :: finfiltp     !
 CHARACTER(len=175) :: flsurp       !
 CHARACTER(len=175) :: fslsurp      !
 CHARACTER(len=175) :: fkvaryp      !
 CHARACTER(len=175) :: fagwrcp       !
 CHARACTER(len=175) :: finfexpp     !PWAT-PARM3    "
 CHARACTER(len=175) :: finfildp     !
 CHARACTER(len=175) :: fdeepfrp     !
 CHARACTER(len=175) :: fbasetpp     !
 CHARACTER(len=175) :: fagwetpp     !
 CHARACTER(len=175) :: fcepscp      !PWAT-PARM4    "
 CHARACTER(len=175) :: fuzsnp       !
 CHARACTER(len=175) :: fnsurp       !
 CHARACTER(len=175) :: fintfwp      !
 CHARACTER(len=175) :: fircp        !
 CHARACTER(len=175) :: flzetpp      !
 CHARACTER(len=175) :: fcepsp       !PWAT-STATE1    "
 CHARACTER(len=175) :: fsursp       !
 CHARACTER(len=175) :: fifwsp       !
 CHARACTER(len=175) :: fuzsp        !
 CHARACTER(len=175) :: flzsp        !
 CHARACTER(len=175) :: fagwsp       !
 CHARACTER(len=175) :: fgwvsp       !
 CHARACTER(len=175)   :: lzsn_std     !PWAT-PARM2 spatiotemporal input specifications
 CHARACTER(len=175)   :: infilt_std   !
 CHARACTER(len=175)   :: lsur_std     !
 CHARACTER(len=175)   :: slsur_std    !
 CHARACTER(len=175)   :: kvary_std    !
 CHARACTER(len=175)   :: agwrc_std    !
 CHARACTER(len=175)   :: infexp_std   !PWAT-PARM3    "
 CHARACTER(len=175)   :: infild_std   !
 CHARACTER(len=175)   :: deepfr_std   !
 CHARACTER(len=175)   :: basetp_std   !
 CHARACTER(len=175)   :: agwetp_std   !
 CHARACTER(len=175)   :: cepsc_std    !PWAT-PARM4    "
 CHARACTER(len=175)   :: uzsn_std     !
 CHARACTER(len=175)   :: nsur_std     !
 CHARACTER(len=175)   :: intfw_std    !
 CHARACTER(len=175)   :: irc_std      !
 CHARACTER(len=175)   :: lzetp_std    !
 CHARACTER(len=175)   :: ceps_std     !PWAT-STATE1    "
 CHARACTER(len=175)   :: surs_std     !
 CHARACTER(len=175)   :: ifws_std     !
 CHARACTER(len=175)   :: uzs_std      !
 CHARACTER(len=175)   :: lzs_std      !
 CHARACTER(len=175)   :: agws_std     !
 CHARACTER(len=175)   :: gwvs_std     !
 CHARACTER(len=175)   :: outfld_ts    !
 INTEGER(I4B) :: ii, it, iit         !ii is a loop index, it: external timestep, iit: internal timestep
 INTEGER(I4B) :: tsnum
 CHARACTER(len=PathLength) :: foptimep ! overpass times [i.e. gridded data output]

 REAL(DP), DIMENSION(10) :: uzra   ! currently unused
 REAL(DP), DIMENSION(10) :: intgrl ! currently unused

 ! distributed output filenames
 CHARACTER(len=20) :: rain_fname  = 'fort.rain_xxxxxxxxxx'
 CHARACTER(len=20) :: msupy_fname = 'fort.msupyxxxxxxxxxx'
 CHARACTER(len=20) :: cepo_fname  = 'fort.cepo_xxxxxxxxxx'
 CHARACTER(len=20) :: ifwi_fname  = 'fort.ifwi_xxxxxxxxxx'
 CHARACTER(len=20) :: uzi_fname   = 'fort.uzi__xxxxxxxxxx'
 CHARACTER(len=20) :: infil_fname = 'fort.infilxxxxxxxxxx'
 CHARACTER(len=20) :: suro_fname  = 'fort.suro_xxxxxxxxxx'    ! surs pre-route
 CHARACTER(len=20) :: ceps_fname  = 'fort.ceps_xxxxxxxxxx'    ! -- initial conditions --
 CHARACTER(len=20) :: sursb_fname = 'fort.sursbxxxxxxxxxx'    ! surs pre-route
 CHARACTER(len=20) :: surs_fname  = 'fort.surs_xxxxxxxxxx'    ! surs pos-route
 CHARACTER(len=20) :: ifws_fname  = 'fort.ifws_xxxxxxxxxx'
 CHARACTER(len=20) :: uzs_fname   = 'fort.uzs__xxxxxxxxxx'
 CHARACTER(len=20) :: lzs_fname   = 'fort.lzs__xxxxxxxxxx'
 CHARACTER(len=20) :: agws_fname  = 'fort.agws_xxxxxxxxxx'
 CHARACTER(len=20) :: gwvs_fname  = 'fort.gwvs_xxxxxxxxxx'

 uzra(1)   = 0.0_dp                   !table of coordinates for function used to evaluate upper zone behavior
 uzra(2)   = 1.25
 uzra(3)   = 1.50
 uzra(4)   = 1.75
 uzra(5)   = 2.00
 uzra(6)   = 2.10
 uzra(7)   = 2.20
 uzra(8)   = 2.25
 uzra(9)   = 2.5
 uzra(10)  = 4.0_dp
 intgrl(1) = 0.0_dp
 intgrl(2) = 1.29
 intgrl(3) = 1.58
 intgrl(4) = 1.92
 intgrl(5) = 2.36
 intgrl(6) = 2.81
 intgrl(7) = 3.41
 intgrl(8) = 3.8
 intgrl(9) = 7.1
 intgrl(10)= 3478.0_dp

 OPEN( UNIT=11, FILE=modelfilep, STATUS='OLD', ACTION='READ', &
       IOSTAT=status)
 IF ( status == 0 ) THEN
   READ(11,'(A)') fgeometap
   READ(11,'(A)') fmaskp
   READ(11,'(A)') fhrusp
   READ(11,'(A)') fk01p                 ! file for 1D kinematic wave parameters
   READ(11,'(A)')   rain_std            ! input file + spatiotemporal distribution
   READ(11,'(A)') petinp_std
   READ(11,'(A)')   lzsn_std            !PWAT-PARM2 filenames & spatiotemporal input specifications
   READ(11,'(A)') infilt_std
   READ(11,'(A)')   lsur_std
   READ(11,'(A)')  slsur_std
   READ(11,'(A)')  kvary_std
   READ(11,'(A)')  agwrc_std
   READ(11,'(A)') infexp_std            !PWAT-PARM3  "
   READ(11,'(A)') infild_std
   READ(11,'(A)') deepfr_std
   READ(11,'(A)') basetp_std
   READ(11,'(A)') agwetp_std
   READ(11,'(A)')  cepsc_std            !PWAT-PARM4  "
   READ(11,'(A)')   uzsn_std
   READ(11,'(A)')   nsur_std
   READ(11,'(A)')  intfw_std
   READ(11,'(A)')    irc_std
   READ(11,'(A)')  lzetp_std
   READ(11,'(A)')   ceps_std            !PWAT-STATE1 "
   READ(11,'(A)')   surs_std
   READ(11,'(A)')   ifws_std
   READ(11,'(A)')    uzs_std
   READ(11,'(A)')    lzs_std
   READ(11,'(A)')   agws_std
   READ(11,'(A)')   gwvs_std
   READ(11,'(A)')  outfld_ts
   READ(11,'(A)')   foptimep

   CLOSE(11)

   CALL split(rain_std,  ' ',frainp)
   CALL split(petinp_std,' ',fpetinpp)
   CALL split(lzsn_std,  ' ',flzsnp)
   CALL split(infilt_std,' ',finfiltp)
   CALL split(lsur_std,  ' ',flsurp)
   CALL split(slsur_std, ' ',fslsurp)
   CALL split(kvary_std, ' ',fkvaryp)
   CALL split(agwrc_std, ' ',fagwrcp)
   CALL split(infexp_std,' ',finfexpp)
   CALL split(infild_std,' ',finfildp)
   CALL split(deepfr_std,' ',fdeepfrp)
   CALL split(basetp_std,' ',fbasetpp)
   CALL split(agwetp_std,' ',fagwetpp)
   CALL split(cepsc_std, ' ',fcepscp)
   CALL split(uzsn_std,  ' ',fuzsnp)
   CALL split(nsur_std,  ' ',fnsurp)
   CALL split(intfw_std, ' ',fintfwp)
   CALL split(irc_std,   ' ',fircp)
   CALL split(lzetp_std, ' ',flzetpp)
   CALL split(ceps_std,  ' ',fcepsp)
   CALL split(surs_std,  ' ',fsursp)
   CALL split(ifws_std,  ' ',fifwsp)
   CALL split(uzs_std,   ' ',fuzsp)
   CALL split(lzs_std,   ' ',flzsp)
   CALL split(agws_std,  ' ',fagwsp)
   CALL split(gwvs_std,  ' ',fgwvsp)

   IF (verbose > 0) THEN
     WRITE (*,*) 'fgeometap: ',fgeometap
     WRITE (*,*) 'fmaskp   : ',fmaskp
     WRITE (*,*) 'fhrusp   : ',fhrusp
     WRITE (*,*) 'frainp   : ',frainp,    '|rain_std   : ',rain_std
     WRITE (*,*) 'fpetinpp : ',fpetinpp,  '|petinp_std : ',petinp_std
     WRITE (*,*) '---'
     WRITE (*,*) 'Parameter files:'
     WRITE (*,*) 'flzsnp   : ',flzsnp,    '|lzsn_std   : ',  lzsn_std
     WRITE (*,*) 'finfiltp : ',finfiltp,  '|infilt_std : ',infilt_std
     WRITE (*,*) 'flsurp   : ',flsurp,    '|lsur_std   : ',  lsur_std
     WRITE (*,*) 'fslsurp  : ',fslsurp,   '|slsur_std  : ', slsur_std
     WRITE (*,*) 'fkvaryp  : ',fkvaryp,   '|kvary_std  : ', kvary_std
     WRITE (*,*) 'fagwrcp  : ',fagwrcp,   '|agwrc_std  : ', agwrc_std
     WRITE (*,*) 'finfexpp : ',finfexpp,  '|infexp_std : ',infexp_std
     WRITE (*,*) 'finfildp : ',finfildp,  '|infild_std : ',infild_std
     WRITE (*,*) 'fdeepfrp : ',fdeepfrp,  '|deepfr_std : ',deepfr_std
     WRITE (*,*) 'fbasetpp : ',fbasetpp,  '|basetp_std : ',basetp_std
     WRITE (*,*) 'fagwetpp : ',fagwetpp,  '|agwetp_std : ',agwetp_std
     WRITE (*,*) 'fcepscp  : ',fcepscp,   '|cepsc_std  : ', cepsc_std
     WRITE (*,*) 'fuzsnp   : ',fuzsnp,    '|uzsn_std   : ',  uzsn_std
     WRITE (*,*) 'fnsurp   : ',fnsurp,    '|nsur_std   : ',  nsur_std
     WRITE (*,*) 'fintfwp  : ',fintfwp,   '|intfw_std  : ', intfw_std
     WRITE (*,*) 'fircp    : ',fircp,     '|irc_std    : ',   irc_std
     WRITE (*,*) 'flzetpp  : ',flzetpp,   '|lzetp_std  : ', lzetp_std
     WRITE (*,*) '---'
     WRITE (*,*) 'Initial Conditions files:'
     WRITE (*,*) 'fcepsp   : ',fcepsp,    '|ceps_std  : ',ceps_std
     WRITE (*,*) 'fsursp   : ',fsursp,    '|surs_std  : ',surs_std
     WRITE (*,*) 'fifwsp   : ',fifwsp,    '|ifws_std  : ',ifws_std
     WRITE (*,*) 'fuzsp    : ',fuzsp,     '|uzs_std  : ',  uzs_std
     WRITE (*,*) 'flzsp    : ',flzsp,     '|lzs_std  : ',  lzs_std
     WRITE (*,*) 'fagwsp   : ',fagwsp,    '|agws_std  : ',agws_std
     WRITE (*,*) 'fgwvsp   : ',fgwvsp,    '|gwvs_std  : ',gwvs_std
   END IF
 ELSE
   WRITE (*,1000) status
    1000 FORMAT (1X,'open mainfile failed--status = ', I6)
   STOP
 END IF

 CALL read_geometa (fgeometap, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'read_geometa error: ', error
 END IF

 IF (dti >= dto) THEN
   dti = dto
   stepdiv = 1.0_dp
   nti = nto
 ELSE
   stepdiv = dto/dti
   nti = nto * INT(stepdiv)
 END IF
 delt60 = dti / 3600.00_dp
 IF (verbose > 0) THEN
   WRITE(*,'(1X,A20,F8.2)') 'MOD(dto, dti) = ',MOD(dto, dti)
   WRITE(*,'(1X,A20,2I10)') 'nto, nti = ', nto, nti
 END IF
 IF (MOD(dto, dti) > 0.0_dp) THEN
   WRITE(*,*) 'error! Please adapt dti [s] to be an exact divisor of dto [s]'
 END IF

 nx = ncols                            ! eastings  into global variable
 ny = nrows                            ! northings  "     "       "

 CALL ReadTimes2IRTS(foptimep, otimes, unit=11, STAT=error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'ModuleHspfPWater:: ReadTimes2IRTS error: ', error
 END IF
 IF (verbose > 0) THEN
   WRITE(*,*) 'otimes:'
   CALL WriteIRTS(otimes)
 END IF


 CALL allo_2Darr (nx, ny, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'global allocation error: ', error
 END IF
 SELECT CASE(simspdist)
   CASE('smdis')
     IF (verbose > 0) THEN
       WRITE(*,'(1X,A11,A5,A)') 'simspdist: ',simspdist,'; semidistributed HSPF run!!!'
     END IF
   CASE('spdis')
     IF (verbose > 0) THEN
       WRITE(*,'(1X,A11,A5,A57)') 'simspdist: ',simspdist,'; distributed HSPF run!!!'
     END IF
     CALL allo_2Darr_spdis (nx, ny, error)
     IF (error > 0) THEN
       WRITE(*,'(1X,A30,I4)') 'allocation for spdis 2Darr error: ', error
     END IF
   CASE DEFAULT
     WRITE(*,'(1X,A11,A5,A)') 'simspdist: ',simspdist,' is not available, nothing will be done!!!'
     RETURN
 END SELECT

 CALL read_logicalASC (fmaskp, mask, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'mask input error: ', error
 END IF
 IF (verbose > 0) THEN
   WRITE(*,*) 'COUNT(mask)=',COUNT(mask)
   WRITE(*,*) 'SHAPE(mask)=',SHAPE(mask)
   WRITE(*,*) 'SIZE(mask)=',SIZE(mask)
   WRITE(*,*) 'SIZE(mask,1)=',SIZE(mask,1)
   WRITE(*,*) 'SIZE(mask,2)=',SIZE(mask,2)
 END IF

 CALL read_intASC (fhrusp, hrus, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'hrus input error: ', error
     !OPEN( UNIT=11, FILE='hrus_caquotas.bin', STATUS='REPLACE', ACTION='WRITE', &                 !just to check back proper input
     !      FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     !WRITE(11) hrus
     !CLOSE(11)
 END IF
 nhrus = MAXVAL(hrus)
 IF (verbose > 0) THEN
   WRITE(*,*) 'nhrus:',nhrus
 END IF
 IF (nhrus == 1) THEN
  tsnum = nhrus
 ELSE
  tsnum = nhrus + 1
 END IF

 CALL allo_tsarr_hrus (nti, tsnum, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'hrus time series alloc error: ', error
 END IF
 ! +++ Read input for time-constant parameters +++
 CALL read_static_dbl (flzsnp, lzsn_std, hrus, nhrus, lzsn, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'lzsn input error: ', error
 END IF

 CALL read_static_dbl (finfiltp, infilt_std, hrus, nhrus, infilt, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'infilt input error: ', error
 END IF

 CALL read_static_dbl (flsurp, lsur_std, hrus, nhrus, lsur, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'lsur input error: ', error
 END IF

 CALL read_static_dbl (fslsurp, slsur_std, hrus, nhrus, slsur, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'slsur input error: ', error
 END IF

 CALL read_static_dbl (fkvaryp, kvary_std, hrus, nhrus, kvary, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'kvary input error: ', error
 END IF

 CALL read_static_dbl (fagwrcp, agwrc_std, hrus, nhrus, agwrc, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'agwrc input error: ', error
 END IF

 CALL read_static_dbl (finfexpp, infexp_std, hrus, nhrus, infexp, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'infexp input error: ', error
 END IF

 CALL read_static_dbl (finfildp, infild_std, hrus, nhrus, infild, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'infild input error: ', error
 END IF

 CALL read_static_dbl (fdeepfrp, deepfr_std, hrus, nhrus, deepfr, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'deepfr input error: ', error
 END IF

 CALL read_static_dbl (fbasetpp, basetp_std, hrus, nhrus, basetp, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'basetp input error: ', error
 END IF

 CALL read_static_dbl (fagwetpp, agwetp_std, hrus, nhrus, agwetp, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'agwetp input error: ', error
 END IF

 CALL read_static_dbl (fcepscp, cepsc_std, hrus, nhrus, cepsc, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'cepsc input error: ', error
 END IF

 CALL read_static_dbl (fuzsnp, uzsn_std, hrus, nhrus, uzsn, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'uzsn input error: ', error
 END IF

 CALL read_static_dbl (fnsurp, nsur_std, hrus, nhrus, nsur, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'nsur input error: ', error
 END IF

 CALL read_static_dbl (fintfwp, intfw_std, hrus, nhrus, intfw, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'intfw input error: ', error
 END IF

 CALL read_static_dbl (fircp, irc_std, hrus, nhrus, irc, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'irc input error: ', error
 END IF

 CALL read_static_dbl (flzetpp, lzetp_std, hrus, nhrus, lzetp, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'lzetp input error: ', error
 END IF

 CALL integrate_to_hru (lzsn, hrus, lzsn_hru)          !obtains mean values for each parameter; more complex scaling could be chosen
 CALL integrate_to_hru (infilt, hrus, infilt_hru)      !for specific parameters
 CALL integrate_to_hru (lsur, hrus, lsur_hru)
 CALL integrate_to_hru (slsur, hrus, slsur_hru)
 CALL integrate_to_hru (kvary, hrus, kvary_hru)
 CALL integrate_to_hru (agwrc, hrus, agwrc_hru)
 CALL integrate_to_hru (infexp, hrus, infexp_hru)
 CALL integrate_to_hru (infild, hrus, infild_hru)
 CALL integrate_to_hru (deepfr, hrus, deepfr_hru)
 CALL integrate_to_hru (basetp, hrus, basetp_hru)
 CALL integrate_to_hru (agwetp, hrus, agwetp_hru)
 CALL integrate_to_hru (cepsc, hrus, cepsc_hru)
 CALL integrate_to_hru (uzsn, hrus, uzsn_hru)
 CALL integrate_to_hru (nsur, hrus, nsur_hru)
 CALL integrate_to_hru (intfw, hrus, intfw_hru)
 CALL integrate_to_hru (irc, hrus, irc_hru)
 CALL integrate_to_hru (lzetp, hrus, lzetp_hru)

 mask_hru = .TRUE.

 IF (verbose > 1) THEN
   WRITE(*,*) 'hrus lzsn:',lzsn_hru
   WRITE(*,*) 'hrus cepsc:',cepsc_hru
   WRITE(*,*) 'hrus uzsn:',uzsn_hru
   WRITE(*,*) 'hrus nsur:',nsur_hru
   WRITE(*,*) 'hrus intfw:',intfw_hru
   WRITE(*,*) 'hrus irc:',irc_hru
   WRITE(*,*) 'hrus lzetp:',lzetp_hru
   WRITE(*,*) 'hrus mask:',mask_hru
 END IF

 ! +++ Read initial conditions to stores +++
 CALL read_static_dbl (fcepsp, ceps_std, hrus, nhrus, ceps, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial ceps input error: ', error
 END IF

 CALL read_static_dbl (fsursp, surs_std, hrus, nhrus, surs, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial surs input error: ', error
 END IF

 CALL read_static_dbl (fifwsp, ifws_std, hrus, nhrus, ifws, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial ifws input error: ', error
 END IF

 CALL read_static_dbl (fuzsp, uzs_std, hrus, nhrus, uzs, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial uzs input error: ', error
 END IF

 CALL read_static_dbl (flzsp, lzs_std, hrus, nhrus, lzs, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial lzs input error: ', error
 END IF

 CALL read_static_dbl (fagwsp, agws_std, hrus, nhrus, agws, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial agws input error: ', error
 END IF

 CALL read_static_dbl (fgwvsp, gwvs_std, hrus, nhrus, gwvs, error)
 IF (error > 0) THEN
   WRITE(*,'(1X,A30,I4)') 'initial gwvs input error: ', error
 END IF

 IF (simspdist  == 'spdis') THEN
   infilt = infilt * delt60
   kgw = 1.0_dp - agwrc**(delt60/24.0_dp)
 ELSE
   CALL integrate_to_hru (ceps, hrus, ceps_hru)       ! just as mean values for each hru
   CALL integrate_to_hru (surs, hrus, surs_hru)
   CALL integrate_to_hru (ifws, hrus, ifws_hru)
   CALL integrate_to_hru (uzs, hrus, uzs_hru)
   CALL integrate_to_hru (lzs, hrus, lzs_hru)
   CALL integrate_to_hru (agws, hrus, agws_hru)
   CALL integrate_to_hru (gwvs, hrus, gwvs_hru)

   infilt_hru = infilt_hru * delt60
   kgw = 1.0_dp - agwrc_hru**(delt60/24.0_dp)

   CALL deallo_2Darr (error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A20,I2)') 'deallo_2Darr error: ',error
   END IF
 END IF

 ! CONSTRUCT requested routing modules
 IF (simspdist == 'spdis') THEN
    WRITE(*,*)'initial surs mass:',SUM(surs)*dx*dx/1000.0d0     ! integrate surs [mm] in the domain -> [m3]
    WHERE(.NOT. mask)
      surs = 0.0d0
    END WHERE
    WRITE(*,*)'masked initial surs mass:',SUM(surs)*dx*dx/1000.0d0
    IF (rwave == 'k01') THEN
      CALL ConstructKin1D(fk01p, mask, error)                                ! allocate & init matrices, load data
    ELSE IF (rwave == 'swe') THEN
      CALL amr2(t0, dti, construct=.TRUE.)
      WRITE(*,*) 'mxtopo(1): ',mxtopo(1)
      WRITE(*,*) 'mytopo(1): ',mytopo(1)
      WRITE(*,*) 'nx: ', nx
      WRITE(*,*) 'ny: ', ny
      IF (mxtopo(1) /= nx .OR. mytopo(1) /= ny) THEN
        STOP 'ModuleHspfPwater - HSPF_Pwater - ERR01'
      END IF
    ELSE
      STOP 'ModuleHspfPwater - HSPF_Pwater - ERR02'
    END IF
 END IF


 ! initialize FLAGS
 smsfg = 0
 fsmsfg = 0

 ! initialize values to cummulative flows  and other variables for adaptable arrays +++
 lzrat  = 0.0_dp                                                   !just to initialize masked areas
 rlzrat = -1.0E30
 lzfrac = -1.0E30

 IF (verbose > 0) THEN
   WRITE(*,*) 'frainp:',frainp
   WRITE(*,*) 'rain_std:',rain_std
   WRITE(*,*) 'nto', nto
   WRITE(*,*) 'nx',  nx
   WRITE(*,*) 'ny',  ny
   WRITE(*,*) 'error:', error
 END IF

 iit  = 0
 WRITE(t0str,'(I10.10)') NINT(t0)
 outer: DO it = 1, nto
 !outer: DO it = 1, 5
   IF (verbose > 0) THEN
     WRITE(*,*) 'starting it', it,' | t0:',t0str
   END IF

   CALL read_rain (frainp, rain_std, it, nto, nx, ny, rain, error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A17,I6,A9,I2)') 'read_rain tstep: ',it,'| error: ',error
   END IF

   CALL read_petinp (fpetinpp, petinp_std, it, nto, nx, ny, petinp, error)
   IF (error > 0) THEN
     WRITE(*,'(1X, A17,I6,A9,I2)') 'read_petinp tstep: ',it,'| error: ',error
   END IF

   IF (verbose > 2) THEN
     WRITE(*,*) 'mean rainfall [mm] for tstep ',it,'|',SUM(rain,MASK=mask)/COUNT(mask)
     WRITE(*,*) 'mean EVTP [mm] for tstep ',it,'|',SUM(petinp,MASK=mask)/COUNT(mask)
   END IF

     !IF(it == 2) THEN
     !  supy = supy * 2.5_dp
     !  OPEN( UNIT=11, FILE='supy_foo.bin', STATUS='REPLACE', ACTION='WRITE', &
     !        FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     !  WRITE(11) supy
     !  CLOSE(11)
     !  WHERE (mask)
     !    WHERE (supy >= petinp)
     !     petinp = petinp * 2.0_dp
     !    ELSEWHERE                                      !note that without the WHERE(Mask) surrounding condition, this would involve supy == NA areas
     !     petinp = petinp + REAL(it)*2.0_dp
     !    END WHERE
     !  ELSEWHERE
     !    petinp = 0.0_dp
     !  END WHERE
     !  OPEN( UNIT=11, FILE='petinp_foo.bin', STATUS='REPLACE', ACTION='WRITE', &
     !        FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     !  WRITE(11) petinp
     !  CLOSE(11)
     !  WRITE(*,*) 'petinp:'
     !  DO j=1, ny
     !    WRITE(*,*) petinp(:,j)
     !  END DO
     !END IF

   rain = rain / stepdiv
   petinp = petinp / stepdiv

   !ii = 0
   inner: DO ii = 1, INT(stepdiv)
     IF (verbose > 0) THEN
       WRITE(*,*) NEW_LINE('A'),'HSPF_Pwater :: integrating t -> t+dt : ',t0,'->',t0+dti
     END IF
     iit = iit + 1
     t0 = t0 + dti                                                     ! t0 now refers to the end time of the current integration loop
     WRITE(t0str,'(I10.10)') NINT(t0)

     IF (verbose > 2) THEN
       WRITE(*,*) 'executing iit:', iit
     END IF
     CALL integrate_to_hru (rain, hrus, rain_hru(iit:iit,:))           ! iit:iit avoids rank dropping. For smdis, this is the mean ts over the watershed


     CALL integrate_to_hru (petinp, hrus, petinp_hru(iit:iit,:))
     IF (verbose > 0) THEN
       WRITE(*,*) 'rain_hru(iit:iit,:):',  rain_hru(iit:iit,:)
       WRITE(*,*) 'rain mass [m3] : ',rain_hru(iit,nhrus+1)*COUNT(mask) * dx**2 / 1000.0D0
       !WRITE(*,*) 'petinp_hru(iit:iit,:):',petinp_hru(iit:iit,:)
     END IF
     SELECT CASE(simspdist)
     CASE('smdis')

       IF (iit > 1) THEN                                                                  !state variables
         ceps_hru(iit,:) = ceps_hru(iit-1,:)
         surs_hru(iit,:) = surs_hru(iit-1,:)
         uzs_hru(iit,:)  = uzs_hru(iit-1,:)
         ifws_hru(iit,:) = ifws_hru(iit-1,:)
         lzs_hru(iit,:)  = lzs_hru(iit-1,:)
         agws_hru(iit,:) = agws_hru(iit-1,:)
         gwvs_hru(iit,:) = gwvs_hru(iit-1,:)                               !gwvs is not an storage but an INOUT variable for every step
       END IF

       supy_hru(iit,:) = rain_hru(iit,:)               !csnfg == 0; i.e., snow melting is not considered (inffac = 1.0_dp (scalar REAL(DP); hard-coded input)
       !WRITE(*,*) 'supy_hru(iit:iit,:):',  supy_hru(iit:iit,:)
       pet_hru(iit,:)  = petinp_hru(iit,:)
       !simulate interception
       CALL icept(supy_hru(iit:iit,:), cepsc_hru, &
                  ceps_hru(iit:iit,:), &
                  cepo_hru(iit:iit,:), mask_hru)
       !WRITE(*,*) 'interception storage|      ceps_hru(iit:iit,:):', ceps_hru(iit:iit,:)
       IF (slifp > 0) THEN
         WRITE(*,*) 'Error! surli input not implemented'            !surli input would be required
       ELSE
         surli_hru(iit,:) = 0.0_dp                                  !surface lateral inflow is not considered
       END IF
       suri_hru(iit,:) = cepo_hru(iit,:) + surli_hru(iit,:)         !surface inflow is the sum of interception outflow and surface lateral inflow (if any)
       msupy = suri_hru(iit:iit,:) + surs_hru(iit:iit,:)            !note this is surs at the start of this ivl (calculated at the end of previous step)
       surss = surs_hru(iit:iit,:)                                  ! "    "        "       "             "        "      "

       WHERE (msupy > 0.0_dp .AND. mask_hru)
         WHERE (smsfg == 0)                                         ! this is the first interval with surface moisture supply after one or more intervals with none
           fsmsfg = 1
         ELSEWHERE                                                  ! this is not the first wet interval
           fsmsfg = 0
         END WHERE
         smsfg = 1                                                  ! there is surface moisture supply
       ELSEWHERE                                                    ! there is no surface moisture supply
         smsfg = 0
         fsmsfg = 0
       END WHERE
       WHERE (mask_hru)
         lzrat = lzs_hru(iit:iit,:) / lzsn_hru
       END WHERE

       IF ( iit == 1) THEN
         WHERE (mask_hru)
           dec_hru = 0.00982*(nsur_hru*lsur_hru/SQRT(slsur_hru))**0.6              ! should not be needed for distributed kinematic routing
           src_hru = 1020.0*(SQRT(slsur_hru)/(nsur_hru*lsur_hru))                  !      "          "             "
         END WHERE
       END IF

       ! simulate behavior of water on the land surface - surface detention, infiltration, interflow input, surface outflow
       CALL surfac (lzrat, infilt_hru, infexp_hru, inffac, infild_hru, fsmsfg, lsur_hru, slsur_hru, vifwfg, &
                   msupy, uzsn_hru, uzs_hru(iit:iit,:), delt60, &
                   uzra, intgrl, rtopfg, uzfg, nsur_hru, intfw_hru, dec_hru, src_hru, &
                   surs_hru(iit:iit,:), infil_hru(iit:iit,:), uzi_hru(iit:iit,:), &
                   ifwi_hru(iit:iit,:),  suro_hru(iit:iit,:), mask_hru)

       ! simulate interflow
       IF (ilifp > 0) THEN
         WRITE(*,*) 'Error! ifwli input not implemented'            !ifwli input would be required
       ELSE
         ifwli_hru(iit,:) = 0.0_dp                              !surface lateral inflow is not considered
       END IF

       CALL intflw (delt60, ifwi_hru(iit:iit,:), ifwli_hru(iit:iit,:), &
                    irc_hru, ifwk1, ifwk2, ifws_hru(iit:iit,:), uzs_hru(iit:iit,:), &
                    ifwo_hru(iit:iit,:), mask_hru)
       IF (verbose > 2) THEN
         WRITE(*,*) 'intflw executed'
         WRITE(*,*) 'interflow storage|         ifws_hru(iit,:):', ifws_hru(iit,:)
       END IF

       ! simulate upper zone behavior
       CALL uzone (uzsn_hru, uzi_hru(iit:iit,:), infilt_hru, inffac, lzrat, &
                   uzs_hru(iit:iit,:), perc_hru(iit:iit,:), mask_hru)
       iperc_hru(iit,:) = perc_hru(iit,:) + infil_hru(iit,:)                   ! collects inflows to lower zone and groundwater

       ! simulate lower zone behavior
       CALL lzone (iperc_hru(iit:iit,:), lzrat, &
                   lzfrac, lzs_hru(iit:iit,:), rlzrat, &
                   lzi_hru(iit:iit,:), &
                   mask_hru)

       ! simulated groundwater behavior
       gwi = iperc_hru(iit:iit,:) - lzi_hru(iit:iit,:)
       IF (alifp > 0) THEN
         WRITE(*,*) 'Error! agwli input not implemented'            !agwli input would be required
       ELSE
         agwli_hru(iit,:) = 0.0_dp                              !groundwater lateral inflow is not considered
       END IF
       CALL gwater (deepfr_hru, gwi, kvary_hru, kgw, agwli_hru(iit:iit,:), &
                    agws_hru(iit:iit,:), gwvs_hru(iit:iit,:), &
                    igwi_hru(iit:iit,:), agwi_hru(iit:iit,:), agwo_hru(iit:iit,:), &
                    mask_hru)
       ! simulate ET
       CALL evapt (pet_hru(iit:iit,:), basetp_hru, uzsn_hru, agwetp_hru, kvary_hru, lzsn_hru, &
                   delt60, agwo_hru(iit:iit,:), ceps_hru(iit:iit,:), uzs_hru(iit:iit,:), agws_hru(iit:iit,:), &
                   gwvs_hru(iit:iit,:), lzetp_hru, rparm, lzs_hru(iit:iit,:), &
                   rempet_hru(iit:iit,:), taet_hru(iit:iit,:), baset_hru(iit:iit,:), cepe_hru(iit:iit,:), &
                   uzet_hru(iit:iit,:), agwet_hru(iit:iit,:), lzet_hru(iit:iit,:), mask_hru)

       pero_hru(iit,:)   = suro_hru(iit,:) + ifwo_hru(iit,:) + agwo_hru(iit,:) ! find total outflow
       watin_hru(iit,:)  = supy_hru(iit,:) + surli_hru(iit,:) + &                  ! total input of water to the pervious land segment
                           ifwli_hru(iit,:)+ agwli_hru(iit,:)
       watdif_hru(iit,:) = watin_hru(iit,:) - &                                        ! net input of water to the pervious land segment
                           (pero_hru(iit,:) + igwi_hru(iit,:) + taet_hru(iit,:))

       pers_hru(iit,:)   = ceps_hru(iit,:) + surs_hru(iit,:) + &                   ! total moisture storage
                           ifws_hru(iit,:) + uzs_hru(iit,:) + &
                           lzs_hru(iit,:) + agws_hru(iit,:)

       !WRITE(*,*) 'semidistributed data for iit:',iit
       !WRITE(*,'(/1X,A,I6)') 'semidistributed data for iit:',iit
       !WRITE(*,*) '-----------------------------'
       !WRITE(*,*) 'rain_hru(iit:iit,:):',  rain_hru(iit:iit,:)
       !WRITE(*,*) 'petinp_hru(iit:iit,:):',petinp_hru(iit:iit,:)
       !WRITE(*,*) 'supy_hru(iit:iit,:):',  supy_hru(iit:iit,:)
       !WRITE(*,'(/1X,A,I6)') 'storages [mm]::'
       !WRITE(*,*) '---'
       !WRITE(*,*) 'interception storage|      ceps_hru(iit:iit,:):', ceps_hru(iit:iit,:)
       !WRITE(*,*) 'surface detention storage| surs_hru(iit:iit,:):', surs_hru(iit:iit,:)
       !WRITE(*,*) 'interflow storage|         ifws_hru(iit:iit,:):', ifws_hru(iit:iit,:)
       !WRITE(*,*) 'upper zone storage|        uzs_hru(iit:iit,:) :', uzs_hru(iit:iit,:)
       !WRITE(*,*) 'lower zone storage|        lzs_hru(iit:iit,:) :', lzs_hru(iit:iit,:)
       !WRITE(*,*) 'active groundwater storage|agws_hru(iit:iit,:):', agws_hru(iit:iit,:)
       !WRITE(*,'(/1X,A)') 'evapotranspiration flows [mm/ivl]::'
       !WRITE(*,*) '---'
       !WRITE(*,*) 'baset_hru(iit:iit,:):',  baset_hru(iit:iit,:)
       !WRITE(*,*) 'cepe_hru(iit:iit,:) :',  cepe_hru(iit:iit,:)
       !WRITE(*,*) 'uzet_hru(iit:iit,:) :',  uzet_hru(iit:iit,:)
       !WRITE(*,*) 'agwet_hru(iit:iit,:):',  agwet_hru(iit:iit,:)
       !WRITE(*,*) 'lzet_hru(iit:iit,:) :',  lzet_hru(iit:iit,:)

 !WRITE(*,*) 'executed iit:', iit

     CASE('spdis')
       supy = rain                                                                       ! snow currently not taken into account
       pet = petinp

       !simulate interception
       CALL icept(supy, cepsc, ceps, cepo, mask)
       IF (verbose > 0) THEN
         WRITE(*,*) 'cepo mass: [m3]',SUM(cepo, MASK=mask) * dx**2 / 1000.0D0
       END IF

       IF (slifp > 0) THEN
         WRITE(*,*) 'Error! surli input not implemented'            !surli input would be required
       ELSE
         surli = 0.0_dp                                             !surface lateral inflow is not considered
       END IF
       suri = cepo + surli                                          !surface inflow is the sum of interception outflow and surface lateral inflow (if any)
       msupy = suri + surs                                          !note this is surs at the start of this ivl (calculated at the end of previous step)
       surss = surs                                                 ! "    "        "       "             "        "      "
       IF (verbose > 0) THEN
         WRITE(*,*) 'surss mass: [m3]',SUM(surss, MASK=mask) * dx**2 / 1000.0D0
         WRITE(*,*) 'msupy mass: [m3]',SUM(msupy, MASK=mask) * dx**2 / 1000.0D0
       END IF
       WHERE (msupy > 0.0_dp .AND. mask)
         WHERE (smsfg == 0)                                         ! this is the first interval with surface moisture supply after one or more intervals with none
           fsmsfg = 1
         ELSEWHERE                                                  ! this is not the first wet interval
           fsmsfg = 0
         END WHERE
         smsfg = 1                                                  ! there is surface moisture supply
       ELSEWHERE                                                    ! there is no surface moisture supply
         smsfg = 0
         fsmsfg = 0
       END WHERE
       WHERE (mask)
         lzrat = lzs / lzsn
       END WHERE

       IF (iit == 1) THEN
         WHERE (mask)
           dec = 0.00982*(nsur*lsur/SQRT(slsur))**0.6              ! should not be needed for distributed kinematic routing
           src = 1020.0*(SQRT(slsur)/(nsur*lsur))                  !      "          "             "
         END WHERE
       END IF

       ! simulate behavior of water on the land surface - surface detention, infiltration, interflow input, surface outflow
       CALL surfac (lzrat, infilt, infexp, inffac, infild, fsmsfg, lsur, slsur, vifwfg, &
                   msupy, uzsn, uzs, delt60, &
                   uzra, intgrl, rtopfg, uzfg, nsur, intfw, dec, src, &
                   surs, infil, uzi, &
                   ifwi,  suro, mask)

       ! simulate interflow
       IF (ilifp > 0) THEN
         WRITE(*,*) 'Error! ifwli input not implemented'            !ifwli input would be required
       ELSE
         ifwli = 0.0_dp                              !surface lateral inflow is not considered
       END IF
       CALL intflw (delt60, ifwi, ifwli, &
                    irc, ifwk1, ifwk2, ifws, uzs, &
                    ifwo, mask)
       IF (verbose > 2) THEN
         WRITE(*,*) 'intflw executed'
         WRITE(*,*) 'interflow storage|         ifws_hru(iit:iit,:):', ifws_hru(iit:iit,:)
       END IF

       ! simulate upper zone behavior
       CALL uzone (uzsn, uzi, infilt, inffac, lzrat, &
                   uzs, perc, mask)
       iperc = perc + infil                   ! collects inflows to lower zone and groundwater

       ! simulate lower zone behavior
       CALL lzone (iperc, lzrat, &
                   lzfrac, lzs, rlzrat, &
                   lzi, &
                   mask)

       ! simulated groundwater behavior
       gwi = iperc - lzi
       IF (alifp > 0) THEN
         WRITE(*,*) 'Error! agwli input not implemented'            !agwli input would be required
       ELSE
         agwli = 0.0_dp                              !groundwater lateral inflow is not considered
       END IF
       CALL gwater (deepfr, gwi, kvary, kgw, agwli, &
                    agws, gwvs, &
                    igwi, agwi, agwo, &
                    mask)

       ! simulate ET
       CALL evapt (pet, basetp, uzsn, agwetp, kvary, lzsn, &
                   delt60, agwo, ceps, uzs, agws, &
                   gwvs, lzetp, rparm, lzs, &
                   rempet, taet, baset, cepe, &
                   uzet, agwet, lzet, mask)

       pero   = suro + ifwo + agwo ! find total outflow
       watin  = supy  + surli + &                  ! total input of water to the pervious land segment
                ifwli + agwli
       watdif = watin - &                                        ! net input of water to the pervious land segment
               (pero + igwi + taet)

       pers   = ceps + surs + &                   ! total moisture storage
                ifws + uzs + &
                lzs + agws

       !  TODO: CALL ...something to write flag=T 2Darrays

       ! integrate 2Darrays into temporal timeseries
       CALL integrate_to_hru (ceps, hrus, ceps_hru(iit:iit,:))
       CALL integrate_to_hru (surs, hrus, surs_hru(iit:iit,:))
       CALL integrate_to_hru (ifws, hrus, ifws_hru(iit:iit,:))
       CALL integrate_to_hru (uzs,  hrus, uzs_hru(iit:iit,:))
       CALL integrate_to_hru (lzs,  hrus, lzs_hru(iit:iit,:))
       CALL integrate_to_hru (agws, hrus, agws_hru(iit:iit,:))
       CALL integrate_to_hru (supy, hrus, supy_hru(iit:iit,:))
       CALL integrate_to_hru (suro, hrus, suro_hru(iit:iit,:))
       CALL integrate_to_hru (ifwo, hrus, ifwo_hru(iit:iit,:))
       CALL integrate_to_hru (agwo, hrus, agwo_hru(iit:iit,:))
       CALL integrate_to_hru (pero, hrus, pero_hru(iit:iit,:))
       CALL integrate_to_hru (igwi, hrus, igwi_hru(iit:iit,:))
       CALL integrate_to_hru (cepe, hrus, cepe_hru(iit:iit,:))
       CALL integrate_to_hru (uzet, hrus, uzet_hru(iit:iit,:))
       CALL integrate_to_hru (lzet, hrus, lzet_hru(iit:iit,:))
       CALL integrate_to_hru (agwet, hrus, agwet_hru(iit:iit,:))
       CALL integrate_to_hru (baset, hrus, baset_hru(iit:iit,:))
       CALL integrate_to_hru (taet, hrus, taet_hru(iit:iit,:))
       CALL integrate_to_hru (ifwi, hrus, ifwi_hru(iit:iit,:))
       CALL integrate_to_hru (uzi, hrus, uzi_hru(iit:iit,:))
       CALL integrate_to_hru (infil, hrus, infil_hru(iit:iit,:))
       CALL integrate_to_hru (perc, hrus, perc_hru(iit:iit,:))
       CALL integrate_to_hru (lzi, hrus, lzi_hru(iit:iit,:))
       CALL integrate_to_hru (agwi, hrus, agwi_hru(iit:iit,:))
       CALL integrate_to_hru (suri, hrus, suri_hru(iit:iit,:))
       CALL integrate_to_hru (rempet, hrus, rempet_hru(iit:iit,:))
       CALL integrate_to_hru (gwvs, hrus, gwvs_hru(iit:iit,:))
       CALL integrate_to_hru (watin, hrus, watin_hru(iit:iit,:))
       CALL integrate_to_hru (watdif, hrus, watdif_hru(iit:iit,:))
       CALL integrate_to_hru (pers, hrus, pers_hru(iit:iit,:))

    END SELECT
   END DO inner
   ! distributed output

   dsno = trim(dsnsim) // '/output/maps/' // trim(pid) // '/'
   !WRITE(*,*) 'dsno: ',trim(dsno)
   CALL which(otimes%value(:,1) == 0.D0 .AND. otimes%time <= t0, aid)              ! pending output
   ifaid: IF (aid(1) > 0) THEN
     WRITE(tostr,'(I10.10)') NINT(otimes%time(MAXVAL(aid)))

     rain_fname(11:20)  = tostr
     cepo_fname(11:20)  = tostr
     msupy_fname(11:20) = tostr
     ifwi_fname(11:20)  = tostr
     uzi_fname(11:20)   = tostr
     infil_fname(11:20) = tostr
     suro_fname(11:20)  = tostr

     ceps_fname(11:20)  = tostr
     sursb_fname(11:20) = tostr
     surs_fname(11:20)  = tostr
     ifws_fname(11:20)  = tostr
     uzs_fname(11:20)   = tostr
     lzs_fname(11:20)   = tostr
     agws_fname(11:20)  = tostr
     gwvs_fname(11:20)  = tostr

     OPEN(UNIT=15, FILE=trim(dsno) // rain_fname,  &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) rain
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // cepo_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) cepo
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // msupy_fname, &
         STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) msupy
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // ifwi_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) ifwi
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // uzi_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) uzi
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // infil_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) infil
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // suro_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) suro
     CLOSE(15)
     ! storages - initial conditions
     OPEN(UNIT=15, FILE=trim(dsno) // ceps_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) ceps
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // surs_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) surs
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // ifws_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) ifws
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // uzs_fname,  & 
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) uzs
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // lzs_fname,  &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) lzs
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // agws_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) agws
     CLOSE(15)
     OPEN(UNIT=15, FILE=trim(dsno) // gwvs_fname, &
          STATUS='unknown', ACTION='WRITE', FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
     WRITE(15) gwvs
     CLOSE(15)

     otimes%value(aid,1) = 1.D0                                                    ! aproximation: set to written
     !CALL WriteIRTS(otimes)
     !STOP 'forced stop'
   END IF ifaid
 END DO outer
 CALL write_ts (outfld_ts, &
                     rain_hru, petinp_hru, &
                     ceps_hru, surs_hru, ifws_hru, uzs_hru, lzs_hru, agws_hru, &
                     supy_hru, suro_hru, ifwo_hru, agwo_hru, pero_hru, igwi_hru, &
                     cepe_hru, uzet_hru, lzet_hru, agwet_hru, baset_hru, taet_hru, &
                     ifwi_hru, uzi_hru, infil_hru, perc_hru, lzi_hru, agwi_hru, &
                     suri_hru, rempet_hru, gwvs_hru, watin_hru, watdif_hru, &
                     pers_hru)

 ! DESTRUCTOR
 IF (simspdist == 'spdis') THEN
   IF (rwave == 'k01') THEN
     CALL DestructKin1D(error)
     IF (error > 0) THEN
       WRITE(*,'(/1X,A23,I2)') 'DestructKin1D | error: ',error
     END IF
   ELSE ! 'swe'
     CALL amr2(t0,dti, destruct=.TRUE.)
   END IF
 END IF
 CALL deallo_common (error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A23,I2)') 'deallo_common | error: ',error
   END IF
 CALL deallo_tsarr_hrus (error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A27,I2)') 'deallo_tsarr_hrus | error: ',error
   END IF
 IF (simspdist == 'spdis') THEN
   CALL deallo_2Darr (error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A20,I2)') 'deallo_2Darr error: ',error
   END IF
   CALL deallo_2Darr_spdis (error)
   IF (error > 0) THEN
     WRITE(*,'(/1X,A28,I2)') 'deallo_tsarr_spdis | error: ',error
   END IF
 END IF

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/topHSPF.err', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='FORMATTED', ACCESS='SEQUENTIAL', IOSTAT=status)
 WRITE(15,*) status
 CLOSE(15)

END SUBROUTINE HSPF_Pwater

SUBROUTINE read_geometa (fgeometap, error)
 ! global metadata
 IMPLICIT NONE

 CHARACTER(len=*), INTENT(IN) :: fgeometap         !file with global data
 INTEGER(I4B), INTENT(OUT) :: error
 INTEGER(I4B) :: status
 OPEN (UNIT=11, FILE=fgeometap, STATUS='OLD', ACTION='READ', &
       IOSTAT=status)
 IF ( status == 0 ) THEN
   READ(11,*)
   READ(11,*) dx, ncols, nrows, t0, dto, dti, nto, simspdist, rwave
   CLOSE(11)
   error = 0
   IF (verbose > 0) THEN
     WRITE (*,1002)  dx, ncols, nrows, t0, dto, dti, nto, simspdist, rwave
     1002 FORMAT (/1X,'simulation geometadata:',/, &
          1X,'dx         = ',F11.2/,&
          1X,'ncols[w-e] = ',I8/,   &
          1X,'nrows[n-s] = ',I8/,   &
          1X,'t0         = ',F14.0/,&
          1X,'dto        = ',F9.0/, &
          1X,'dti     = ',F9.0/, &
          1X,'nto     = ',I8/,   &
          1X,'simspdist  = ',A5/,   &
          1X,'rwave      = ',A3)
   END IF
 ELSE
   WRITE (*,1000) status
    1000 FORMAT (1X,'open fgeometa failed--status = ', I6)
   error = 1
 END IF
END SUBROUTINE read_geometa

SUBROUTINE allo_2Darr (nx, ny, error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN)    :: nx, ny
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B), DIMENSION(31) :: status

 ALLOCATE(lzsn(nx,ny),     STAT=status(1))          !PWAT-PARM2
 ALLOCATE(infilt(nx,ny),   STAT=status(2))          !PWAT-PARM2
 ALLOCATE(lsur(nx,ny),     STAT=status(3))          !PWAT-PARM2
 ALLOCATE(slsur(nx,ny),    STAT=status(4))          !PWAT-PARM2
 ALLOCATE(kvary(nx,ny),    STAT=status(5))          !PWAT-PARM2
 ALLOCATE(agwrc(nx,ny),    STAT=status(6))          !PWAT-PARM2
 ALLOCATE(infexp(nx,ny),   STAT=status(7))          !PWAT-PARM3
 ALLOCATE(infild(nx,ny),   STAT=status(8))          !PWAT-PARM3
 ALLOCATE(deepfr(nx,ny),   STAT=status(9))          !PWAT-PARM3
 ALLOCATE(basetp(nx,ny),   STAT=status(10))         !PWAT-PARM3
 ALLOCATE(agwetp(nx,ny),   STAT=status(11))         !PWAT-PARM3
 ALLOCATE(cepsc(nx,ny),    STAT=status(12))         !PWAT-PARM4
 ALLOCATE(uzsn(nx,ny),     STAT=status(13))         !PWAT-PARM4
 ALLOCATE(nsur(nx,ny),     STAT=status(14))         !PWAT-PARM4
 ALLOCATE(intfw(nx,ny),    STAT=status(15))         !PWAT-PARM4
 ALLOCATE(irc(nx,ny),      STAT=status(16))         !PWAT-PARM4
 ALLOCATE(lzetp(nx,ny),    STAT=status(17))         !PWAT-PARM4
 ALLOCATE(ceps(nx,ny),     STAT=status(18))         !PWAT-STATE1 & PWAT internally calculated time series
 ALLOCATE(surs(nx,ny),     STAT=status(19))         ! "
 ALLOCATE(ifws(nx,ny),     STAT=status(20))         ! "
 ALLOCATE(uzs(nx,ny),      STAT=status(21))         ! "
 ALLOCATE(lzs(nx,ny),      STAT=status(22))         ! "
 ALLOCATE(agws(nx,ny),     STAT=status(23))         ! "
 ALLOCATE(gwvs(nx,ny),     STAT=status(24))         ! "

 ALLOCATE(mask(nx,ny),     STAT=status(25))
 ALLOCATE(hrus(nx,ny),     STAT=status(26))
 ALLOCATE(rain(nx,ny),     STAT=status(27))
 ALLOCATE(petinp(nx,ny),   STAT=status(28))
 ALLOCATE(surli(nx,ny),    STAT=status(29))
 !ALLOCATE(uzli(nx,ny),     STAT=status(23))
 ALLOCATE(ifwli(nx,ny),    STAT=status(30))
 !ALLOCATE(lzli(nx,ny),     STAT=status(25))
 ALLOCATE(agwli(nx,ny),    STAT=status(31))

 IF (sum(status) == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
 !WRITE(*,'(A,I6)') 'Allocation sum(status): ', sum(status)
END SUBROUTINE allo_2Darr

SUBROUTINE allo_2Darr_spdis (nx, ny, error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN)    :: nx, ny
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B), DIMENSION(41) :: status

  status = 0
  ALLOCATE(pers(nx,ny),   STAT=status(1))
  !ALLOCATE(petadj(nx,ny), STAT=status(9))

  ALLOCATE(supy(nx,ny),   STAT=status(2))
  ALLOCATE(suro(nx,ny),   STAT=status(3))
  ALLOCATE(ifwo(nx,ny),   STAT=status(4))
  ALLOCATE(agwo(nx,ny),   STAT=status(5))
  ALLOCATE(pero(nx,ny),   STAT=status(6))
  ALLOCATE(igwi(nx,ny),   STAT=status(7))
  ALLOCATE(pet(nx,ny),    STAT=status(8))
  ALLOCATE(cepe(nx,ny),   STAT=status(9))
  ALLOCATE(uzet(nx,ny),   STAT=status(10))
  ALLOCATE(lzet(nx,ny),   STAT=status(11))
  ALLOCATE(agwet(nx,ny),  STAT=status(12))
  ALLOCATE(baset(nx,ny),  STAT=status(13))
  ALLOCATE(taet(nx,ny),   STAT=status(14))
  ALLOCATE(ifwi(nx,ny),   STAT=status(15))
  ALLOCATE(uzi(nx,ny),    STAT=status(16))
  ALLOCATE(infil(nx,ny),  STAT=status(17))
  ALLOCATE(perc(nx,ny),   STAT=status(18))
  ALLOCATE(lzi(nx,ny),    STAT=status(19))
  ALLOCATE(agwi(nx,ny),   STAT=status(20))
  ALLOCATE(suri(nx,ny),   STAT=status(21))

  ALLOCATE(cepo(nx,ny),   STAT=status(22))
  ALLOCATE(rempet(nx,ny), STAT=status(23))
  ALLOCATE(iperc(nx,ny),  STAT=status(24))
  ALLOCATE(dec(nx,ny),    STAT=status(25))
  ALLOCATE(src(nx,ny),    STAT=status(26))
  ALLOCATE(watin(nx,ny),  STAT=status(27))
  ALLOCATE(watdif(nx,ny), STAT=status(28))

  !arrays adaptable to smdis or spdis simulations
  ALLOCATE(maskloc(nx,ny),STAT=status(29))   ! deallocated in deallo_common
  ALLOCATE(ifwk1(nx,ny),  STAT=status(30))
  ALLOCATE(ifwk2(nx,ny),  STAT=status(31))
  ALLOCATE(lzrat(nx,ny),  STAT=status(32))
  ALLOCATE(rlzrat(nx,ny), STAT=status(33))
  ALLOCATE(lzfrac(nx,ny), STAT=status(34))
  ALLOCATE(gwi(nx,ny),    STAT=status(35))
  ALLOCATE(kgw(nx,ny),    STAT=status(36))
  ALLOCATE(rparm(nx,ny),  STAT=status(37))
  ALLOCATE(msupy(nx,ny),  STAT=status(38))
  ALLOCATE(surss(nx,ny),  STAT=status(39))
  ALLOCATE(smsfg(nx,ny),  STAT=status(40))
  ALLOCATE(fsmsfg(nx,ny), STAT=status(41))

  IF (verbose > 0) THEN
    WRITE(*,*) 'allo_2Darr_spdis status:', status
  END IF
  IF (sum(status) == 0) THEN
    error = 0
  ELSE
    error = 1
  END IF
END SUBROUTINE allo_2Darr_spdis

SUBROUTINE allo_tsarr_hrus (nti, tsnum, error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(IN)    :: nti, tsnum
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B), DIMENSION(74) :: status

 status = 0

 ALLOCATE(lzsn_hru(1, tsnum),     STAT=status(1))          !PWAT-PARM2
 ALLOCATE(infilt_hru(1, tsnum),   STAT=status(2))          !PWAT-PARM2
 ALLOCATE(lsur_hru(1, tsnum),     STAT=status(3))          !PWAT-PARM2
 ALLOCATE(slsur_hru(1, tsnum),    STAT=status(4))          !PWAT-PARM2
 ALLOCATE(kvary_hru(1, tsnum),    STAT=status(5))          !PWAT-PARM2
 ALLOCATE(agwrc_hru(1, tsnum),    STAT=status(6))          !PWAT-PARM2
 ALLOCATE(infexp_hru(1, tsnum),   STAT=status(7))          !PWAT-PARM3
 ALLOCATE(infild_hru(1, tsnum),   STAT=status(8))          !PWAT-PARM3
 ALLOCATE(deepfr_hru(1, tsnum),   STAT=status(9))          !PWAT-PARM3
 ALLOCATE(basetp_hru(1, tsnum),   STAT=status(10))         !PWAT-PARM3
 ALLOCATE(agwetp_hru(1, tsnum),   STAT=status(11))         !PWAT-PARM3
 ALLOCATE(cepsc_hru(1, tsnum),    STAT=status(12))         !PWAT-PARM4
 ALLOCATE(uzsn_hru(1, tsnum),     STAT=status(13))         !PWAT-PARM4
 ALLOCATE(nsur_hru(1, tsnum),     STAT=status(14))         !PWAT-PARM4
 ALLOCATE(intfw_hru(1, tsnum),    STAT=status(15))         !PWAT-PARM4
 ALLOCATE(irc_hru(1, tsnum),      STAT=status(16))         !PWAT-PARM4
 ALLOCATE(lzetp_hru(1, tsnum),    STAT=status(17))         !PWAT-PARM4

 ALLOCATE(mask_hru(1, tsnum),           STAT=status(18))

 ALLOCATE(rain_hru(nti, tsnum),    STAT=status(19))                        ! TIME SERIES INPUT
 ALLOCATE(petinp_hru(nti, tsnum),  STAT=status(20))                        !  "
 ALLOCATE(surli_hru(nti, tsnum),   STAT=status(21))                        !  "
 status(22) = 0 !ALLOCATE(uzli_hru(nti, tsnum),    STAT=status(22))
 ALLOCATE(ifwli_hru(nti, tsnum),   STAT=status(23))                        !  "
 status(24) = 0 !ALLOCATE(lzli_hru(nti, tsnum),    STAT=status(24))
 ALLOCATE(agwli_hru(nti, tsnum),   STAT=status(25))                        !  "

 ALLOCATE(pers_hru(nti, tsnum),     STAT=status(26))                       !  TIME SERIES COMPUTED BY PERLND
 ALLOCATE(ceps_hru(nti, tsnum),     STAT=status(27))                       !  "
 ALLOCATE(surs_hru(nti, tsnum),     STAT=status(28))                       !  "
 ALLOCATE(ifws_hru(nti, tsnum),     STAT=status(29))                       !  "
 ALLOCATE(uzs_hru(nti, tsnum),      STAT=status(30))                       !  "
 ALLOCATE(lzs_hru(nti, tsnum),      STAT=status(31))                       !  "
 ALLOCATE(agws_hru(nti, tsnum),     STAT=status(32))                       !  "
 ALLOCATE(gwvs_hru(nti, tsnum),     STAT=status(33))                       !  "
 !ALLOCATE(petadj_hru(nti, tsnum),   STAT=status(34))

 ALLOCATE(supy_hru(nti, tsnum),     STAT=status(35))                       !  "
 ALLOCATE(suro_hru(nti, tsnum),     STAT=status(36))                       !  "
 ALLOCATE(ifwo_hru(nti, tsnum),     STAT=status(37))                       !  "
 ALLOCATE(agwo_hru(nti, tsnum),     STAT=status(38))                       !  "
 ALLOCATE(pero_hru(nti, tsnum),     STAT=status(39))                       !  "
 ALLOCATE(igwi_hru(nti, tsnum),     STAT=status(40))                       !  "
 ALLOCATE(pet_hru(nti, tsnum),      STAT=status(41))                       !  "
 ALLOCATE(cepe_hru(nti, tsnum),     STAT=status(42))                       !  "
 ALLOCATE(uzet_hru(nti, tsnum),     STAT=status(43))                       !  "
 ALLOCATE(lzet_hru(nti, tsnum),     STAT=status(44))                       !  "
 ALLOCATE(agwet_hru(nti, tsnum),    STAT=status(45))                       !  "
 ALLOCATE(baset_hru(nti, tsnum),    STAT=status(46))                       !  "
 ALLOCATE(taet_hru(nti, tsnum),     STAT=status(47))                       !  "
 ALLOCATE(ifwi_hru(nti, tsnum),     STAT=status(48))                       !  "
 ALLOCATE(uzi_hru(nti, tsnum),      STAT=status(49))                       !  "
 ALLOCATE(infil_hru(nti, tsnum),    STAT=status(50))                       !  "
 ALLOCATE(perc_hru(nti, tsnum),     STAT=status(51))                       !  "
 ALLOCATE(lzi_hru(nti, tsnum),      STAT=status(52))                       !  "
 ALLOCATE(agwi_hru(nti, tsnum),     STAT=status(53))                       !  "
 ALLOCATE(suri_hru(nti, tsnum),     STAT=status(54))                       !  "

 ALLOCATE(cepo_hru(nti, tsnum),     STAT=status(55))                       !  TIME SERIES COMPUTED BY PERLND NOT INCLUDED IN THE CATALOG
 ALLOCATE(rempet_hru(nti, tsnum),   STAT=status(56))                       !  "
 ALLOCATE(iperc_hru(nti, tsnum),    STAT=status(57))                       !  "
 ALLOCATE(dec_hru(1, tsnum),             STAT=status(58))                       !  "
 ALLOCATE(src_hru(1, tsnum),             STAT=status(59))                       !  "
 ALLOCATE(watin_hru(nti, tsnum),    STAT=status(60))                       !  "
 ALLOCATE(watdif_hru(nti, tsnum),   STAT=status(61))                       !  "

 IF(.NOT. ALLOCATED(maskloc)) THEN
   ALLOCATE(maskloc(1, tsnum),STAT=status(62))
 END IF
 IF(.NOT. ALLOCATED(ifwk1)) THEN
   ALLOCATE(ifwk1(1, tsnum),STAT=status(63))
 END IF
 IF(.NOT. ALLOCATED(ifwk2)) THEN
   ALLOCATE(ifwk2(1, tsnum),STAT=status(64))
 END IF
 IF(.NOT. ALLOCATED(lzrat)) THEN
   ALLOCATE(lzrat(1, tsnum),STAT=status(65))
 END IF
 IF(.NOT. ALLOCATED(rlzrat)) THEN
   ALLOCATE(rlzrat(1, tsnum),STAT=status(66))
 END IF
 IF(.NOT. ALLOCATED(lzfrac)) THEN
   ALLOCATE(lzfrac(1, tsnum),STAT=status(67))
 END IF
 IF(.NOT. ALLOCATED(gwi)) THEN
   ALLOCATE(gwi(1, tsnum),STAT=status(68))
 END IF
 IF(.NOT. ALLOCATED(kgw)) THEN
   ALLOCATE(kgw(1, tsnum),STAT=status(69))
 END IF
 IF(.NOT. ALLOCATED(rparm)) THEN
   ALLOCATE(rparm(1, tsnum),STAT=status(70))
 END IF
 IF(.NOT. ALLOCATED(msupy)) THEN
   ALLOCATE(msupy(1, tsnum),STAT=status(71))
 END IF
 IF(.NOT. ALLOCATED(surss)) THEN
   ALLOCATE(surss(1, tsnum),STAT=status(72))
 END IF
 IF(.NOT. ALLOCATED(smsfg)) THEN
   ALLOCATE(smsfg(1, tsnum),STAT=status(73))
 END IF
 IF(.NOT. ALLOCATED(fsmsfg)) THEN
   ALLOCATE(fsmsfg(1, tsnum),STAT=status(74))
 END IF

 IF (SUM(status) == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
END SUBROUTINE allo_tsarr_hrus

SUBROUTINE deallo_2Darr (error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B), DIMENSION(4) :: status

 status = 0

 DEALLOCATE(lzsn, infilt, lsur, slsur, kvary, agwrc, STAT=status(1)) !PWAT-PARM2
 DEALLOCATE(infexp, infild, deepfr, basetp, agwetp,  STAT=status(2))  !PWAT-PARM3
 DEALLOCATE(cepsc, uzsn, nsur , intfw, irc, lzetp,   STAT=status(3))  !PWAT-PARM4
 DEALLOCATE(ceps, surs, ifws , uzs, lzs, agws, gwvs, STAT=status(4))  !PWAT-STATE1 & internally recalculated

 ! note: mask, hrus, rain, petinp, surli, ifwli, agwli ALLOCATED in allo_2Darr
 !       but DEALLOCATED in deallo_common
 IF (sum(status) == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
END SUBROUTINE deallo_2Darr

SUBROUTINE deallo_common (error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B) :: status

 DEALLOCATE(maskloc, ifwk1, ifwk2, lzrat, rlzrat, lzfrac, &      ! allocated in allo_2Darr_spdis or allo_tsarr_spdis
            gwi, kgw, rparm, msupy, surss, smsfg, fsmsfg, &      !    "               "                    "
            mask, hrus, rain, petinp, surli, ifwli, agwli, &     ! allocated in allo_2Darr
            STAT=status)

 IF (status == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
END SUBROUTINE deallo_common

SUBROUTINE deallo_2Darr_spdis (error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B) :: status

 DEALLOCATE(pers, &
            supy, suro, ifwo, agwo, pero, igwi, &
            pet, cepe, uzet, lzet, agwet, baset, taet, &
            ifwi, uzi, infil, perc, lzi, agwi, suri, &
            cepo, rempet, iperc, dec, src, watin, watdif, &
            STAT=status)
 ! maskloc,...,fsmsfg allocated in allo_2Darr_spdis is deallocated in deallo_common
 IF (status == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
END SUBROUTINE deallo_2Darr_spdis

SUBROUTINE deallo_tsarr_hrus (error)
 IMPLICIT NONE
 INTEGER(I4B), INTENT(OUT)   :: error
 INTEGER(I4B) :: status

 DEALLOCATE(lzsn_hru, infilt_hru, lsur_hru, slsur_hru, &
            kvary_hru, agwrc_hru, infexp_hru, &
            infild_hru, deepfr_hru, basetp_hru, &
            agwetp_hru, cepsc_hru, uzsn_hru, &
            nsur_hru, intfw_hru, irc_hru, lzetp_hru, &
            mask_hru, rain_hru, petinp_hru, surli_hru, &
            ifwli_hru, agwli_hru, pers_hru, &
            ceps_hru, surs_hru, ifws_hru, uzs_hru, &
            lzs_hru, agws_hru, gwvs_hru, &
            supy_hru, suro_hru, ifwo_hru, agwo_hru, pero_hru, &
            igwi_hru, pet_hru, cepe_hru, uzet_hru, &
            lzet_hru, agwet_hru, baset_hru, taet_hru, &
            ifwi_hru, uzi_hru, infil_hru, perc_hru, &
            lzi_hru, agwi_hru, suri_hru, &
            cepo_hru, rempet_hru, iperc_hru, dec_hru, src_hru, &
            watin_hru, watdif_hru, STAT=status)

 IF (status == 0) THEN
  error = 0
 ELSE
  error = 1
 END IF
END SUBROUTINE deallo_tsarr_hrus


SUBROUTINE read_static_dbl (fvarp, std, hrus, nhrus, retmap, error)           !return an external DOUBLE map
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN)  :: fvarp
 CHARACTER(len=*), INTENT(IN)  :: std                                                ! 'dc','lc','sc'| 'dc' <- binary input; 'lc' & 'sc' <- ASCII
 INTEGER(I4B), DIMENSION(:,:), INTENT(IN)  :: hrus                                    !allocatable & previously allocated array
 INTEGER(I4B),     INTENT(IN)  :: nhrus
 REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: retmap                                     !allocatable & previously allocated array
 INTEGER(I4B),     INTENT(OUT) :: error
 INTEGER(I4B) :: status, i
 REAL(DP), DIMENSION(nhrus)           :: tmpin                                       !automatic array

 error = 0

 SELECT CASE (std)
 CASE('dc')                                                                          !distributed constant
   OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &                         !binary
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
   IF ( status == 0 ) THEN
     READ(11) retmap
     CLOSE(11)
   ELSE
     error = error + 1
   END IF
 CASE('lc')                                                                          !lumped constant
   OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &                         !ASCII
         IOSTAT=status)
   IF ( status == 0 ) THEN
     READ(11,*) retmap(1,1)
     retmap = retmap(1,1)
     CLOSE(11)
   ELSE
     error = error + 1
   END IF
 CASE ('sc')                                                                         !semidistributed constant
   OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &                         !ASCII
         IOSTAT=status)
   IF ( status == 0 ) THEN
     READ(11,*) tmpin
     CLOSE(11)
     retmap = 0.0_dp
     DO i=1,nhrus
       WHERE (hrus == i)
         retmap = tmpin(i)
       END WHERE
     END DO
   ELSE
     error = error + 1
   END IF
 CASE DEFAULT
     error = error +1
 END SELECT
END SUBROUTINE read_static_dbl

SUBROUTINE read_static_int (fvarp, std, retmap, error)        !returns an external INTEGER map
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN)  :: fvarp
 CHARACTER(len=*), INTENT(IN)  :: std
 INTEGER(I4B), DIMENSION(:,:), INTENT(OUT) :: retmap
 INTEGER(I4B),     INTENT(OUT) :: error
 INTEGER(I4B) :: status

 error = 0
 SELECT CASE (std)
 CASE('dc')                                                                          !distributed constant
   OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
   IF ( status == 0 ) THEN
     READ(11) retmap
     CLOSE(11)
   ELSE
     error = error + 1
   END IF
 CASE('lc')
   OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         IOSTAT=status)
   IF ( status == 0 ) THEN
     READ(11,*) retmap(1,1)
     retmap = retmap(1,1)
     CLOSE(11)
   ELSE
     error = error + 1
   END IF
  CASE DEFAULT
    error = error +1
  END SELECT
END SUBROUTINE read_static_int

SUBROUTINE read_logicalASC(flogicp, map, error)
 !read an ASCII [0s & 1s] INTEGER map, as returns it as LOGICAL
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN)         :: flogicp
 LOGICAL, DIMENSION(:,:), INTENT(OUT) :: map
 INTEGER(I4B),            INTENT(OUT) :: error

 ! local
 INTEGER(I4B) :: status
 INTEGER(I4B), DIMENSION(SIZE(map,1),SIZE(map,2)) :: tmp2D                         !automatic array

 IF (verbose > 0) THEN
   WRITE(*,*)'read_logicalASC: flogicp:',flogicp
 END IF
 OPEN( UNIT=11, FILE=flogicp, STATUS='OLD', ACTION='READ', &
       IOSTAT=status)
 IF ( status == 0 ) THEN
     READ(11,*) tmp2D
     map = tmp2D == 1
     CLOSE(11)
     error = 0
 ELSE
     error = 1
 END IF
END SUBROUTINE read_logicalASC

SUBROUTINE read_intASC(fvarp, map, error)
 !read an external ASCII INTEGER map
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN)              :: fvarp
 INTEGER(I4B), DIMENSION(:,:), INTENT(OUT) :: map
 INTEGER(I4B),     INTENT(OUT)             :: error

 ! local
 INTEGER(I4B) :: status

 OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
       IOSTAT=status)
 IF ( status == 0 ) THEN
     READ(11,*) map
     CLOSE(11)
     error = 0
 ELSE
     error = 1
 END IF
END SUBROUTINE read_intASC

SUBROUTINE integrate_to_hru (arr2D, hrus, arrhrus)
 IMPLICIT NONE
 REAL(DP),     DIMENSION(:,:), INTENT(IN)  :: arr2D
 INTEGER(I4B), DIMENSION(:,:), INTENT(IN)  :: hrus
 REAL(DP),     DIMENSION(:,:), INTENT(OUT) :: arrhrus
 INTEGER(I4B) :: i, nhrus

 nhrus = MAXVAL(hrus)
 FORALL(i=1:nhrus)
   arrhrus(1,i) =  SUM(arr2D,hrus==i)/COUNT(hrus==i)                  ! mean (lumped) value for each hru
 END FORALL
 IF (nhrus > 1) THEN
   arrhrus(1,nhrus+1) =  SUM(arr2D,hrus /= 0)/COUNT(hrus /= 0)        ! area-weighted mean value for the watershed
 END IF
END SUBROUTINE integrate_to_hru

SUBROUTINE read_rain (fvarp, std, it, nto_s, nx, ny, map, error)
  !rain just allowed to be 'lts' or 'dts'
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: fvarp
 CHARACTER(len=*), INTENT(IN) :: std                                              !spatiotemporal distribution
 INTEGER(I4B),     INTENT(IN) :: it
 INTEGER(I4B),     INTENT(IN) :: nto_s
 INTEGER(I4B),     INTENT(IN) :: nx, ny
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: map
 INTEGER(I4B), INTENT(OUT)    :: error
 !local variables
 INTEGER(I4B) :: status, i
 CHARACTER(len=150), ALLOCATABLE, DIMENSION(:), SAVE :: binmaps       !points to path/names of binary 2D-array stack
 REAL(DP),           ALLOCATABLE, DIMENSION(:), SAVE :: tmp           !allocatable array
 REAL(SP), DIMENSION(nx,ny) :: tmp2D                        !automatic array to read single precision binary rainfall maps

 error = 0

 IF ( it == 1 ) THEN
   SELECT CASE (std)
   CASE ('lts')                   !lumped time series ASCII input. File is data itself
     IF(.NOT. ALLOCATED(tmp)) THEN
       ALLOCATE(tmp(nto_s),STAT=status)
       error = error + status
     END IF
     OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         IOSTAT=status)
     error = error + status
     IF (error == 0) THEN
       READ(11,*) tmp
       CLOSE(11)
     END IF
   CASE ('dts')                   !lumped time series ASCII input. File points to binary 2D data arrays
     IF(.NOT. ALLOCATED(binmaps)) THEN
       ALLOCATE(binmaps(nto_s),STAT=status)
       error = error + status
     END IF
     OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         IOSTAT=status)
     error = error + status
     IF (error == 0) THEN
       READ(11,'(A)') (binmaps(i), i=1,nto_s)
       CLOSE(11)
       IF (verbose > 2) THEN
         WRITE (*,'(A)') (binmaps(i), i=1,nto_s)
       END IF
     END IF
   CASE DEFAULT
     WRITE(*,*) 'Ilegal spatiotemporal distribution specified for rainfall'
     error = 1
   END SELECT
 END IF

 SELECT CASE (std)
 CASE ('lts')
   map = tmp(it)
   IF (verbose > 2) THEN
     WRITE(*,*) 'it=',it,'|..',fvarp((LEN_TRIM(fvarp)-10):LEN_TRIM(fvarp)),tmp(it)
   END IF
 CASE ('dts')
   OPEN( UNIT=11, FILE=binmaps(it), STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
   error = error + status
   IF ( error == 0 ) THEN
     READ(11) tmp2D
     map = DBLE(tmp2D)
     CLOSE(11)
   END IF
 END SELECT

END SUBROUTINE read_rain

SUBROUTINE read_petinp (fvarp, std, it, nto_s, nx, ny, map, error)
 !petinp just allowed to be 'lts' or 'dts'
 IMPLICIT NONE
 CHARACTER(len=*), INTENT(IN) :: fvarp
 CHARACTER(len=*), INTENT(IN) :: std                                              !spatiotemporal distribution
 INTEGER(I4B),     INTENT(IN) :: it
 INTEGER(I4B),     INTENT(IN) :: nto_s
 INTEGER(I4B),     INTENT(IN) :: nx, ny
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: map
 INTEGER(I4B), INTENT(OUT)    :: error
 !local variables
 INTEGER(I4B) :: status, i
 CHARACTER(len=150), ALLOCATABLE, DIMENSION(:), SAVE :: binmaps       !points to path/names of binary 2D-array stack
 REAL(DP),           ALLOCATABLE, DIMENSION(:), SAVE :: tmp           !allocatable array

 error = 0

 IF ( it == 1 ) THEN
   SELECT CASE (std)
   CASE ('lts')                   !lumped time series ASCII input. File is data itself
     IF(.NOT. ALLOCATED(tmp)) THEN
       ALLOCATE(tmp(nto_s),STAT=status)
       error = error + status
     END IF
     OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         IOSTAT=status)
     error = error + status
     IF (error == 0) THEN
       READ(11,*) tmp
       CLOSE(11)
     END IF
   CASE ('dts')                   !lumped time series ASCII input. File points to binary 2D data arrays
     IF(.NOT. ALLOCATED(binmaps)) THEN
       ALLOCATE(binmaps(nto_s),STAT=status)
       error = error + status
     END IF
     OPEN( UNIT=11, FILE=fvarp, STATUS='OLD', ACTION='READ', &
         IOSTAT=status)
     error = error + status
     IF (error == 0) THEN
       READ(11,'(A)') (binmaps(i), i=1,nto_s)
       CLOSE(11)
       IF (verbose > 2) THEN
         WRITE (*,'(A)') (binmaps(i), i=1,nto_s)
       END IF
     END IF
   CASE DEFAULT
     WRITE(*,*) 'Ilegal spatiotemporal distribution specified for potential evapotranspiration'
     error = 1
   END SELECT
 END IF

 SELECT CASE (std)
 CASE ('lts')
   map = tmp(it)
   IF (verbose > 2) THEN
     WRITE(*,*) 'it=',it,'|..',fvarp((LEN_TRIM(fvarp)-10):LEN_TRIM(fvarp)),tmp(it)
   END IF
 CASE ('dts')
   OPEN( UNIT=11, FILE=binmaps(it), STATUS='OLD', ACTION='READ', &
         FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
   error = error + status
   IF ( error == 0 ) THEN
     READ(11) map
     CLOSE(11)
   END IF
 END SELECT

 !IF ( it == nto_s ) THEN
 !  DEALLOCATE(binmaps, tmp, STAT=status)
 !  error = error + status
 !END IF
END SUBROUTINE read_petinp

SUBROUTINE icept (supy, cepsc, &
                  ceps, &
                  cepo, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: supy, cepsc
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: ceps
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: cepo
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask

 WHERE (mask)
   ceps = ceps + supy
   WHERE (ceps > cepsc)
     cepo = ceps - cepsc
     ceps = cepsc
   ELSEWHERE
     cepo = 0.0_dp
   END WHERE
 ELSEWHERE
   ceps = 0.0_dp
   cepo = 0.0_dp
 END WHERE

 !WRITE(*,*) 'ceps in icept after act:', ceps
 !WRITE(*,*) 'cepo in icept after act:', cepo
END SUBROUTINE icept

SUBROUTINE surfac (lzrat, infilt, infexp, inffac, infild, fsmsfg, &
                   lsur, slsur, vifwfg, msupy, uzsn, uzs, delt60, &
                   uzra, intgrl, rtopfg, uzfg, nsur, intfw, &
                   dec, src, surs, &
                   infil, uzi, ifwi, suro, mask)
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: lzrat, infilt, infexp
  REAL(DP), INTENT(IN)                     :: inffac                      !par
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: infild                      !par
  INTEGER(I4B), DIMENSION(:,:), INTENT(IN) :: fsmsfg                      !var FLAG
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: lsur, slsur                 !par
  INTEGER(I4B), INTENT(IN)                 :: vifwfg                      !flag PWAT-PARM1
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: msupy, uzsn, uzs            !var
  REAL(DP), INTENT(IN)                     :: delt60                      !par
  REAL(DP), DIMENSION(:), INTENT(IN)       :: uzra, intgrl                !par  currently unused
  INTEGER(I4B), INTENT(IN)                 :: rtopfg, uzfg                !flag PWAT-PARM1  currenlty just (0, 1)
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: nsur, intfw                 !par
  REAL(DP), DIMENSION(:,:), INTENT(IN)     :: dec, src                    !var
  REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: surs                        !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: infil, uzi, ifwi, suro      !var
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask                        !par
  !local variables
  REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: ibar, imax, imin, ratio

  ! +++ purpose +++
  !Distribute the water available for infiltration and runoff -
  !units of fluxes are mm/ivl

  !     Establish locations of sloping lines on infiltration/inflow/sur
  !     runoff figure.  Prefix "i" refers to "infiltration" line
  !     Ibar is the mean infiltration capacity over the segment
  !     internal units of 'infilt' are mm/ivl
  WHERE (mask .AND. msupy > 0.0_dp)
    ibar = infilt / (lzrat**infexp)
    !WHERE (inffac < 1.0_dp)                                  ! never happens; inffac is set as scalar parameter == 1.0_dp
    !  ibar = ibar * inffac
    !END WHERE
    imax = infild * ibar                                      ! find the maximum and minimum infiltration capacities. infild input parameter = ratio of maximum to mean infiltration capacity
    imin = ibar - (imax - ibar)                               ! for spatially distributed, normally infild == 1
    ratio = intfw * (2.0_dp**lzrat)
    WHERE (ratio <= 1.0_dp)
      ratio = 1.0001_dp
    END WHERE
  END WHERE
  CALL dispos (imin, imax, ratio, msupy, uzsn, uzs, &
               delt60, dec, src, uzra, &
               intgrl, rtopfg, uzfg, &
               surs, &
               infil, uzi, ifwi, suro, mask)
  WHERE (.NOT. (mask .AND. msupy > 0.0_dp))
         surs  = 0.0_dp
         suro  = 0.0_dp
         ifwi  = 0.0_dp
         infil = 0.0_dp
         uzi   = 0.0_dp
  END WHERE
END SUBROUTINE surfac

SUBROUTINE dispos (imin, imax, ratio, msupy, uzsn, uzs, &
                   delt60, dec, src, uzra, &
                   intgrl, rtopfg, uzfg, &
                   surs, &
                   infil, uzi, ifwi, suro, mask)
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: imin, imax, ratio, msupy, uzsn, uzs
  REAL(DP), INTENT(IN) :: delt60                     !par
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: dec        !var
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: src        !var
  REAL(DP), DIMENSION(:), INTENT(IN) :: uzra, intgrl !par  currently unused
  INTEGER(I4B), INTENT(IN) :: rtopfg                 !flag PWAT-PARM1 currently just 0 : new method
  INTEGER(I4B), INTENT(IN) :: uzfg                   !flag PWAT-PARM1 currently just 1, so uzra,intgrl are unused
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: surs    !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: infil     !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: uzi       !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: ifwi      !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: suro      !var
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask   !par
  ! local variables
  REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: pdro, uzfrac, iimin, iimax, psur, &
                                                    undrii, pifwi
  ! + + + PURPOSE + + +
  !     Dispose of moisture supply on either an individual block of the land segment or on an entire land segment.
  !
  !     determine how much of the MSUPY falls above and below the infiltration line in the infiltration/interflow/surface runoff figure
  CALL divisn (imin, imax, msupy, &                                    ! get pdro & infil
               pdro, infil, mask)                                      ! split msupy into pdro & infil (infiltrated moisture), where pdro >= 0.0

  IF (uzfg == 0) THEN
  !  CALL uzinf1(pdro, uzsn, uzs, uzra, intgrl, &
  !             uzi, mask)
  ELSE
    CALL uzinf2(pdro, uzsn, uzs, &                                     ! get uzi
                uzi, mask)
  END IF
  WHERE (mask .AND. msupy > 0.0_dp)
    WHERE (uzi > pdro) uzi = pdro
    uzfrac = uzi/pdro
    iimin = imin*ratio
    iimax = imax*ratio
  END WHERE
  CALL divisn (iimin, iimax, msupy, &
               psur, undrii, mask)                       !get psur: potential surface detention/runoff. Always > 0.0

  IF (verbose > 0) THEN
    WRITE(*,*) 'dispos :: pos-divisn2 msupy  mass [m3]:', SUM(msupy, MASK=mask) * dx**2 / 1000.0D0
    WRITE(*,*) 'dispos :: pos-divisn2 psur   mass [m3]:', SUM(psur,  MASK=mask) * dx**2 / 1000.0D0
    WRITE(*,*) 'dispos :: pos-divisn2 undrii mass [m3]:', SUM(undrii,  MASK=mask) * dx**2 / 1000.0D0
  END IF
  pifwi = pdro - psur
  WHERE (mask .AND. msupy > 0.0_dp)
    ifwi = pifwi*(1.0_dp - uzfrac)
    psur = psur *(1.0_dp - uzfrac)
  END WHERE
  IF (verbose > 0) THEN 
    WRITE(*,*) 'dispos :: pre-route psur mass [m3]:', SUM(psur,  MASK=mask) * dx**2 / 1000.0D0
  END IF
  SELECT CASE(simspdist)
    CASE('smdis')
      CALL proute(psur, rtopfg, delt60, dec, src, &      ! distribute psur between surs storage + suro
              surs, &
              suro, mask)
    CASE('spdis')
      psur = msupy - (ifwi + uzi + infil)
      IF (verbose > 0) THEN
        WRITE(*,*) 'dispos :: pre-drouteK01 psur mass [m3]:', SUM(psur,  MASK=mask) * dx**2 / 1000.0D0
      END IF
      IF (rwave == 'k01') THEN
         CALL drouteK01(t0,dti,psur,surs,suro,mask)    ! connector to 1D cascade kinematic routing
         ! CALL KinematicWave(t0,dti,p
        ! CALL KineDiffWave(psur, dti, dtm, surs, suro, hrus)
      ELSE ! 'swe' - clawpack
        CALL drouteSWE(t0,dti,psur,surs,suro,mask)    ! connector to clawpack SWE
      END IF
  END SELECT

  WHERE (mask .AND. psur <= 0.0_dp)
    surs = 0.0_dp
    suro = 0.0_dp
  END WHERE
  WHERE (mask)
    WHERE (pdro <= 0.0_dp)
      pdro = 0.0_dp
      surs = 0.0_dp
      suro = 0.0_dp
      ifwi = 0.0_dp
      uzi  = 0.0_dp
    END WHERE
  END WHERE
END SUBROUTINE dispos

SUBROUTINE uzinf2 (pdro, uzsn, uzs, uzi, mask)
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: pdro       !var
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: uzsn       !var
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: uzs        !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: uzi       !var
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask   !par
  !local variables
  REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2))     :: uzrat, k1, k2, uzfrac
  LOGICAL, DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: maskloc

  ! + + + PURPOSE + + +
  ! Compute inflow to upper zone during this interval, using "fully forward" type algorithm  as used in HSPX,ARM and NPS.
  ! Note:  although this method should give results closer to those produced by HSPX, etc , its output will be more sensitive to
  ! the timestep than that given by subroutine uzinf1. TODO: perhaps uzinf1 may improve results for DA experiments updating uzs.
  ! Alternatively, sub-hourly time steps could improve the transfer of soil moisture along the profile
  maskloc = mask .AND. pdro > 0.0_dp
  WHERE(maskloc)
   uzrat = uzs/uzsn
   WHERE (uzrat < 2.0_dp)                                             ! piece-wise Gaussian type decay function with more compact support
     k1 = 3.0_dp - uzrat                                              ! e.g. uzfrac=0.5 for uzrat=2, and uzfrac=0 for uzrat=4
     uzfrac = 1.0_dp - (uzrat*0.5_dp)*((1.0_dp/(4.0_dp - uzrat))**k1)
   ELSEWHERE
     k2 = 2.0_dp * uzrat - 3.0_dp
     uzfrac = (1.0_dp/(1.0_dp + k2))**k2
   END WHERE
   uzi = pdro * uzfrac
  ELSEWHERE
    uzi = 0.0_dp
  END WHERE
END SUBROUTINE uzinf2

SUBROUTINE   divisn (min, max, msupy, over, under, mask)
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: min       !var
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: max       !var
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: msupy     !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: over     !var
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: under    !var
  LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask  !par

!     + + + PURPOSE + + +
!     Determine the divisions of the quantities of the moisture supply above and below the "infiltration" line or the "infiltration +
!     interflow" line in the infiltration/interflow/surface runoff figure.  This routine is used either to simulate the behavior of
!     an individual block of the land segment or the entire land segment.
 WHERE (mask)
   WHERE (msupy <= min)                                 ! msupy line is entirely below other line
    under = msupy
    over = 0.0_dp
   ELSEWHERE
     WHERE (msupy > max)                                ! msupy line is entirely above other line
       under = (min + max) * 0.5_dp
       over = msupy - under
     ELSEWHERE                                          ! msupy line crosses the other line
       over = 0.5_dp * ((msupy - min)**2) / (max - min)
       under = msupy - over
     END WHERE
   END WHERE
 ELSEWHERE
   over  = 0.0_dp
   under = 0.0_dp
 END WHERE
END SUBROUTINE divisn

SUBROUTINE etbase (basetp, agwo, rempet, taet, baset, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)    :: basetp           ! [-] fraction of remaining potential ET which can be satisfied from baseflow, if enough is available, parameter PWAT-PARM3
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: agwo             !
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rempet           !
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: taet             !
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: baset            !
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          !
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: baspet
 !     + + + PURPOSE + + +
 !     Simulate ET from baseflow
 WHERE (basetp > 0.0_dp .AND. mask)                          ! there is et from baseflow
   baspet = basetp * rempet
   WHERE (baspet > agwo)                                     ! baseflow et is limited by quantity available
     baset = agwo
     agwo = 0.0_dp
   ELSEWHERE                                                 ! baseflow et will not exhaust storage, so empty at potential
     baset = baspet
     agwo = agwo - baset
   END WHERE
   taet = taet + baset                                       ! update total where etbase is ocurring
   rempet = rempet - baset                                   ! update total where etbase is ocurring
 ELSEWHERE
   baset = 0.0_dp                                            ! no totals (taet, rempet) to update and no et from baseflow
 END WHERE
END SUBROUTINE etbase

SUBROUTINE evicep (ceps, rempet, taet, cepe, mask)           ! evapotranspiration from interception storage
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: ceps             !actual interception storage
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rempet           !remaining potential evapotranspiration
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: taet             !total actual evapotranspiration
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: cepe             !actual evapotranspiration from ceps
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          !mask where evicep is applied
 !     + + + PURPOSE + + +
 !     Simulate ET from interception storage
 WHERE (ceps > 0.0_dp .AND. mask)                             ! there is something in interception storage to evaporate
   WHERE (rempet > ceps)                                      ! evaporation from interception storage is limited by quantity available
     cepe = ceps
     ceps = 0.0_dp
   ELSEWHERE                                                  ! interception evaporation will not exhaust storage, so empty at potential
     cepe = rempet
     ceps = ceps - cepe
   END WHERE
   taet   = taet + cepe                                       ! update total where evicep is ocurring
   rempet = rempet - cepe                                     ! update total where evicep is ocurring
 ELSEWHERE
   cepe = 0.0_dp                                              ! no totals to update, and nothing evaporates from ceps
 END WHERE
END SUBROUTINE evicep

SUBROUTINE etuzon (uzsn, rempet, uzs, taet, uzet, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)    :: uzsn             ! [mm]      upper zone nominal storage, parameter PWAT-PARM2
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rempet           ! [mm/ivl]  rempet is remaining potential ET
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: uzs              ! [mm]      actual upper zone storage
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: taet             ! [mm/ivl]  total simulated ET
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: uzet             ! [mm/ivl]  ET from upper zone
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          ! [bool]    mask where etuzon applies
 !     + + + PURPOSE + + +
 !    Simulate ET from the upper zone
 CALL etuzs (uzsn, rempet, uzs, uzet, mask)
 WHERE (mask)
   taet   = taet + uzet                                  ! update totals
   rempet = rempet - uzet
 END WHERE
END SUBROUTINE etuzon

SUBROUTINE etuzs (uzsn, rempet, uzs, uzet, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)    :: uzsn             ! [mm]      upper zone nominal storage, parameter PWAT-PARM2
 REAL(DP), DIMENSION(:,:), INTENT(IN)    :: rempet           ! [mm/ivl]  rempet is remaining potential ET
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: uzs              ! [mm]      actual upper zone storage
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: uzet             ! [mm/ivl]  ET from upper zone
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          ! [bool]    mask where etuzs applies
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: uzpet, uzrat
 !     + + + PURPOSE + + +
 !    This is a subsidiary subroutine for computing upper zone ET
 WHERE (uzs > 0.0254 .AND. mask)                                   ! there is et from the upper zone, estimate the uzet opportunity | note units are mm
   uzrat = uzs/uzsn
   WHERE (uzrat > 2.0)
     uzpet = rempet
   ELSEWHERE
     uzpet = 0.5*uzrat*rempet
   END WHERE
   WHERE (uzpet > uzs)                                        ! calculate the actual et; upper zone et is limited by quantity available
     uzet = uzs
     uzs = 0.0_dp
   ELSEWHERE                                                  ! upper zone et will not exhaust storage, so empty at potential
     uzet = uzpet
     uzs  = uzs - uzet
   END WHERE
 ELSEWHERE
   uzet = 0.0_dp                                              !uzs is not modified
 END WHERE
END SUBROUTINE etuzs

SUBROUTINE etagw (agwetp, kvary, rempet, agws, taet, gwvs, agwet, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: agwetp              ! fraction of the remaining potential ET than can be sought from active groundwater, parameter PWAT-PARM3
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: kvary               ! [1/mm]  parameter, PWAT-PARM2
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rempet           !
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: agws             !
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: taet             ! total actual evapotranspiration
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: gwvs             !
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: agwet            !
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          ! mask where etagw is applied
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: gwpet
 !     + + + PURPOSE + + +
 !     Simulate ET from active groundwater storage.
 WHERE (agwetp > 0.0_dp .AND. mask)                          ! there is et from groundwater; determine remaining capacity
   gwpet = rempet * agwetp
   WHERE (gwpet > agws)                                      ! groundwater et is limited by quantity available
     agwet = agws
     agws = 0.0_dp
   ELSEWHERE                                                 ! groundwater et will not exhaust storage, so empty at potential
     agwet = gwpet
     agws  = agws - agwet
   END WHERE
   WHERE (ABS(kvary) > 0.0_dp) gwvs = gwvs - agwet            !update variables, where ET occurs from active groundwater
   taet = taet + agwet                                        !  "
   rempet = rempet - agwet                                    !  "
 ELSEWHERE
   agwet = 0.0_dp                                             ! no totals to update, and there is no et from groundwater
 END WHERE
END SUBROUTINE etagw

SUBROUTINE etlzon (lzsn, delt60, lzetp, rempet, rparm, lzs, taet, lzet, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: lzsn                ! [mm]      lower zone nominal storage, parameter PWAT-PARM2
 REAL(DP), INTENT(IN) ::  delt60                             ! [h]       simulation time interval in hours
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: lzetp               ! [-]       lower zone ET parameter, PWAT-PARM4. Note: if monthly-timevarying lzetpm were provided, this would became INOUT
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rempet           ! [mm/ivl]  rempet is remaining potential ET
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: rparm            ! [mm/ivl]  rparm is max ET opportunity (current value of maximum lower zone)
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: lzs              ! [mm]      actual lower zone storage
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: taet             ! [mm/ivl]  total simulated ET
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: lzet             ! [mm/ivl]  ET from lower zone
 LOGICAL, DIMENSION(:,:), INTENT(IN)   :: mask          ! [bool]    mask where etlzon applies
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: lzpet
 !     + + + PURPOSE + + +
 !     Simulate ET from the lower zone.
 WHERE ((ABS(1.0_dp - lzetp) >= 1.0E-5) .AND. mask)
   rparm = 0.25_dp/(1.0_dp - lzetp)*(lzs/lzsn)*delt60/24.0_dp             ![mm/ivl]
 ELSEWHERE
   rparm = 0.0_dp
 END WHERE
 WHERE ((rempet > 0.0_dp) .AND. (lzs > 0.5_dp) .AND. mask)                !there remains an ET demand, and assumes ET can take place (lzs > 0.5 mm)
   WHERE (ABS(1.0_dp - lzetp) < 1.0E-5)                                   !special case - will try to draw ET from whole land segment at remaining potential rate
    lzpet = rempet
   ELSEWHERE                                                              !usual case - desired et will vary over the whole land seg
     WHERE (rempet > rparm)                                               !potential exceeds opportunity
       lzpet = 0.5_dp * rparm
     ELSEWHERE                                                            !potential exceeds opportunity over only part of the land segment
       lzpet = rempet*(1.0 - rempet/(2.0_dp*rparm))
     END WHERE
     WHERE (lzetp < 0.5) lzpet = lzpet*2.0*lzetp
   END WHERE
   WHERE (lzpet < (lzs * 0.95))                                           !lower zone ET will not exhaust storage, so empty at potential
     lzet = lzpet
   ELSEWHERE
     lzet = lzs * 0.95                                                    !lower zone ET is limited by quantity available (assumed as 95% of storage)
   END WHERE
   lzs    = lzs - lzet                                                    !update INOUT variables
   taet   = taet + lzet                                                   ! "
   rempet = rempet - lzet                                                 ! "
 ELSEWHERE
   lzet = 0.0_dp
 END WHERE
END SUBROUTINE etlzon

SUBROUTINE evapt (pet, basetp, uzsn, agwetp, kvary, lzsn, &
                  delt60, agwo, ceps, uzs, agws, &
                  gwvs, lzetp, rparm, lzs, &
                  rempet, taet, baset, cepe, &
                  uzet, agwet, lzet, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: pet, basetp, uzsn, agwetp, kvary, lzsn
 REAL(DP), INTENT(IN) :: delt60
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: lzetp
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: agwo, ceps, uzs, agws, gwvs, rparm, lzs
 REAL(DP), DIMENSION(:,:), INTENT(OUT)   :: rempet, taet, baset, cepe, uzet, agwet, lzet
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 !local variables
 LOGICAL, DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: maskloc1, maskloc2, maskloc3, maskloc4
! + + + PURPOSE + + +
! Simulate evapotranspiration (ET)
 taet = 0.0_dp                                                                  ! [mm/ivl] taet is total actual et
 rempet = pet

 maskloc1 = mask .AND. (rempet > 0.0_dp)
 CALL etbase (basetp, agwo, rempet, taet, baset, maskloc1)                      !etbase will apply just on maskloc = .TRUE.
 maskloc2 = mask .AND. (rempet > 0.0_dp)
 CALL evicep (ceps, rempet, taet, cepe, maskloc2)
 maskloc3 = mask .AND. (rempet > 0.0_dp)
 CALL etuzon (uzsn, rempet, uzs, taet, uzet, maskloc3)
 maskloc4 = mask .AND. (rempet > 0.0_dp)
 CALL etagw (agwetp, kvary, rempet, agws, taet, gwvs, agwet, maskloc4)

 WHERE (mask .AND. (.NOT. maskloc4))
   agwet = 0.0_dp
 END WHERE
 WHERE (mask .AND. (.NOT. maskloc3))
   uzet  = 0.0_dp
   agwet = 0.0_dp
 END WHERE
 WHERE (mask .AND. (.NOT. maskloc2))
   cepe  = 0.0_dp
   uzet  = 0.0_dp
   agwet = 0.0_dp
 END WHERE
 WHERE (mask .AND. (.NOT. maskloc1))
   baset = 0.0_dp
   cepe  = 0.0_dp
   uzet  = 0.0_dp
   agwet = 0.0_dp
 END WHERE

 maskloc4 = mask .AND. (rempet > 0.0_dp)
 CALL etlzon (lzsn, delt60, lzetp, rempet, rparm, lzs, taet, lzet, maskloc4)
 WHERE (mask .AND. (.NOT. maskloc4)) lzet = 0.0_dp

END SUBROUTINE evapt

SUBROUTINE intflw (delt60, ifwi,ifwli, &
                   irc, ifwk1, ifwk2, ifws, uzs, &
                   ifwo, mask)
 IMPLICIT NONE
 REAL(DP), INTENT(IN) :: delt60
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: ifwi, ifwli, irc
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: ifwk1, ifwk2, ifws, uzs
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: ifwo
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: kifw, inflo, value
 ! + + + PURPOSE + + +
 ! Simulate interflow
 WHERE (mask)
  kifw = -LOG(irc) * delt60 / 24.0_dp
  ifwk2 = 1.0 - EXP(-kifw)
  ifwk1 = 1.0 - (ifwk2/kifw)
  inflo = ifwi + ifwli                                                         ! nblks = 1
  value = inflo + ifws
  WHERE (value > 0.005)                                                        ! [mm] there is something worth routing
    ifwo = ifwk1*inflo + ifwk2*ifws
    ifws = value - ifwo
  ELSEWHERE                                                                    !nothing worth routing-dump back to uzs
    ifwo = 0.0_dp
    !ifws = 0.0_dp
    !uzs  = uzs + value
    ifws = value                                                               !line differs original HSPF code
  END WHERE
 uzs = uzs                                                                     !nothing but avoid warning at compiling
 ELSEWHERE
   ifwo = 0.0_dp      !iwfs, uzs untouched
 END WHERE
END SUBROUTINE intflw

SUBROUTINE uzone (uzsn, uzi, infilt, inffac, lzrat, &
                  uzs, &
                  perc, &
                  mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: uzsn, uzi, infilt, lzrat
 REAL(DP), INTENT(IN) :: inffac
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: uzs
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: perc
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 ! + + + PURPOSE + + +
 ! Simulate upper zone behavior
 CALL uzones (uzsn, uzi, infilt, inffac, lzrat, &                                !nblks == 1
              uzs,  &
              perc, &
              mask)
END SUBROUTINE uzone

SUBROUTINE uzones (uzsn, uzi, infilt, inffac, lzrat, &
                  uzs, &
                  perc, &
                  mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: uzsn, uzi, infilt, lzrat
 REAL(DP), INTENT(IN) :: inffac
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: uzs
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: perc
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: uzrat
 ! + + + PURPOSE + + +
 ! This is a subsidiary subroutine for computing upper zone behavior
 WHERE (mask)
   uzrat = uzs/uzsn                        ! percolation will be based on UZRAT at the start of the interval
   uzs = uzs + uzi                         ! add inflow to UZS
   WHERE ((uzrat - lzrat) > (0.01 * 25.4))
     perc = 0.1 * infilt * inffac * uzsn * (uzrat - lzrat)**3  !simulate percolation. Units of perc are mm/ivl
     WHERE (perc > uzs)
       perc = uzs                          !computed value is too high so merely empty storage
       uzs = 0.0_dp
     ELSEWHERE
      uzs = uzs - perc                     !computed value is ok
     END WHERE
   ELSEWHERE
     perc = 0.0_dp                         !assume there is no percolation
   END WHERE
 ELSEWHERE
   perc = 0.0_dp                            !out of watershed
 END WHERE
END SUBROUTINE uzones

SUBROUTINE lzone (iperc, lzrat, &
                  lzfrac, lzs, rlzrat, &
                  lzi, &
                  mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: iperc, lzrat
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: lzfrac, lzs, rlzrat
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: lzi
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: indx
 ! + + + PURPOSE + + +
 ! Simulate lower zone behavior
 WHERE (mask .AND. iperc > 0.0_dp)     !if necessary, recalculate the fraction of infiltration plus percolation which will be taken by lower zone
   WHERE (ABS(lzrat - rlzrat) > 0.02_dp)  !it is time to recalculate
     rlzrat = lzrat
     WHERE (lzrat <= 1.0_dp)
       indx   = 2.5_dp - 1.5_dp*lzrat
       lzfrac = 1.0_dp - lzrat*(1.0_dp/(1.0_dp + indx))**indx
     ELSEWHERE
       indx   = 1.5_dp * lzrat - 0.5_dp
       lzfrac = (1.0_dp/(1.0_dp + indx))**indx
     END WHERE
   END WHERE
   lzi = lzfrac * iperc                ! lower zone inflow
   lzs = lzs + lzi
 ELSEWHERE
   lzi = 0.0_dp        !either no inflow, or out of watershed
 END WHERE
END SUBROUTINE lzone

SUBROUTINE gwater (deepfr, gwi, kvary, kgw, agwli, &
                  agws, gwvs, &
                  igwi, agwi, agwo, &
                  mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN)     :: deepfr, gwi, kvary, kgw, agwli
 REAL(DP), DIMENSION(:,:), INTENT(INOUT)  :: agws, gwvs
 REAL(DP), DIMENSION(:,:), INTENT(OUT)    :: igwi, agwi, agwo
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
 !local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: ainflo
 ! + + + PURPOSE + + +
 ! Simulate groundwater behavior
 WHERE (mask)
   WHERE (gwi > 0.0_dp)                      !groundwater inflow components
     igwi = deepfr * gwi
     agwi = gwi - igwi
   ELSEWHERE
     igwi = 0.0_dp
     agwi = 0.0_dp
   END WHERE
   ainflo = agwi + agwli                     !active groundwater
   agwo = 0.0_dp
   WHERE (ABS(kvary) > 0.0_dp )              !evaluate groundwater recharge parameter
     gwvs = gwvs + ainflo                    !update the index to variable groundwater slope
     WHERE (agws > 1.0E-20) agwo = kgw*(1.0_dp + kvary*gwvs)*agws
   ELSEWHERE                                 !groundwater outflow(baseflow)
     WHERE (agws > 1.0E-20) agwo = kgw * agws
   END WHERE
   agws = agws + (ainflo - agwo)
 ELSEWHERE
   agws = 0.0_dp
   gwvs = 0.0_dp
   igwi = 0.0_dp
   agwi = 0.0_dp
   agwo = 0.0_dp
 END WHERE
END SUBROUTINE gwater

SUBROUTINE proute(psur, rtopfg, delt60, dec, src, &                               !just apply on maskloc2
                surs, &
                suro, mask)
 IMPLICIT NONE
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: psur        ! var potential surface detention/runoff [mm/itvl
 INTEGER(I4B), INTENT(IN) :: rtopfg                  !flag PWAT-PARM1. Currently just 0 (new method for routing)
 REAL(DP), INTENT(IN) :: delt60                      !par
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: dec, src    !var
 REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: surs    !var
 REAL(DP), DIMENSION(:,:), INTENT(OUT) :: suro      !var
 LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask   !par
 ! local variables
 REAL(DP), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: ssupr, sursm, dummy, surse, &
                                                   tsuro, sursnw, change, fact, &
                                                   a1, sterm, fsuro, ffact, dfact, &
                                                   ratio, dfsuro, dterm, dsuro
 INTEGER(I4B), DIMENSION(SIZE(mask,1),SIZE(mask,2)) :: counter
 INTEGER(I4B) :: i, j
 ! + + + PURPOSE + + +
 ! Determine how much potential surface detention (PSUR) runs-off
 ! in one simulation interval.
 SELECT CASE (rtopfg)
   CASE(0)                                              ! do routing the new hspf way: estimate the rate of supply to the overland flow surface
   !msk2_1: FORALL (i=1:SIZE(mask,1), j=1:SIZE(mask,2), (mask(i,j) .AND. psur(i,j) > 0.0002_dp * 25.4_dp))
   northDO: DO j = 1, SIZE(mask,2)
     !$OMP PARALLEL DO PRIVATE(j)
     eastDO: DO i = 1, SIZE(mask,1)
       msk2_1: IF (mask(i,j) .AND. (psur(i,j) > 0.0002_dp * 25.4_dp)) THEN
         ssupr(i,j) = (psur(i,j) - surs(i,j)) /delt60
         surse(i,j) = 0.0_dp                              ! determine equilibrium depth for this supply rate
         IF (ssupr(i,j) > 0.0_dp) THEN
           surse(i,j) = dec(i,j) * ssupr(i,j)**0.6
         END IF
         sursnw(i,j)  = psur(i,j)                         !determine runoff by iteration - Newton's method
         suro(i,j)    = 0.0_dp                            !estimate the new surface storage
         counter(i,j) = 0                                   !do until relative error is small
         newton: DO
           IF (ssupr(i,j) > 0.0_dp) THEN
             ratio(i,j) = sursnw(i,j) / surse(i,j)
             IF (ratio(i,j) <= 1.0_dp) THEN
               fact(i,j) = 1.0_dp + 0.6_dp*ratio(i,j)**3  !flow is increasing
             ELSE
               fact(i,j) = 1.6_dp
             END IF
           ELSE
             ratio(i,j) = 1.0E30
             fact(i,j)  = 1.6_dp
           END IF
           a1(i,j) = delt60 * src(i,j) * fact(i,j)**1.667_dp
           sterm(i,j)  = sursnw(i,j) ** 1.667_dp
           counter(i,j)  = counter(i,j) + 1
           ffact(i,j)  = a1(i,j) * sterm(i,j)
           fsuro(i,j)  = ffact(i,j) - suro(i,j)
           dfact(i,j)  = -1.667_dp * ffact(i,j)
           dfsuro(i,j) = dfact(i,j)/sursnw(i,j) - 1.0_dp
           IF (ratio(i,j) <= 1.0_dp) THEN                      ! additional term required in derivative wrt suro
             dterm(i,j) = dfact(i,j) / (fact(i,j) * surse(i,j)) * 1.8_dp * ratio(i,j)**2
             dfsuro(i,j) = dfsuro(i,j) + dterm(i,j)
           END IF
           dsuro(i,j) = fsuro(i,j)/dfsuro(i,j)
           suro(i,j) = suro(i,j) - dsuro(i,j)
           IF (suro(i,j) <= 1.0E-10) THEN                       !boundary condition - don't let suro go negative
             suro(i,j) = 0.0_dp
           END IF
           sursnw(i,j) = psur(i,j) - suro(i,j)
           change(i,j) = 0.0_dp
           IF (suro(i,j) > 0.0_dp) THEN
             change(i,j) = ABS(dsuro(i,j)/suro(i,j))
           END IF
           IF(change(i,j) < 0.01 .OR. counter(i,j) > 100) EXIT !test loop exit
         END DO newton
         !TODO: checkings about counter array
         surs(i,j) = sursnw(i,j)
       END IF msk2_1
     END DO eastDO
     !$OMP END PARALLEL DO
   END DO northDO
   !END FORALL msk2_1
   WHERE (mask .AND. psur <= 0.0002_dp * 25.4_dp)
       suro = psur
       surs = 0.0_dp
   END WHERE
   WHERE (.NOT. mask)
     suro = 0.0_dp
     surs = 0.0_dp
   END WHERE
   CASE(1)                                      !do routing the way it is done in arm, nps, and hspx: estimate the rate of supply to the overland flow surface
   msk1_1: WHERE (mask)                         ![mm/ivl]
     msk1_2: WHERE (psur > 0.0002_dp * 25.4_dp) !something is worth routing on the surface
       ssupr = psur - surs
       sursm = (surs + psur)*0.5                !estimate the mean detention storage over the interval
       msk1_3: WHERE (ssupr > 0.0_dp)
         dummy = dec * ssupr**0.6               !preliminary estimate of surse
         msk1_4: WHERE (dummy > sursm)
           surse = dummy                        !flow is increasing
           dummy = sursm*(1.0_dp + 0.6_dp*(sursm/surse)**3)
         ELSEWHERE msk1_4
           dummy = sursm*1.6                    !flow on surface is at equilibrium or receding
         END WHERE msk1_4
       ELSEWHERE msk1_3                         !flow on the surface is receding - equilibrium detention is assume equal to actual detention
         dummy = sursm * 1.6
       END WHERE msk1_3
       tsuro = delt60*src*dummy**1.667
       msk1_5: WHERE (tsuro > psur)
         suro = psur                            !too much surface runoff is estimated
         surs = 0.0_dp
       ELSEWHERE msk1_5
         suro = tsuro
         surs = psur - suro
       END WHERE msk1_5
     ELSEWHERE msk1_2                           !send what is on the overland flow plane straight to the channel
       suro = psur
       surs = 0.0_dp
     END WHERE msk1_2
   ELSEWHERE msk1_1
     suro = 0.0_dp
     surs = 0.0_dp
   END WHERE msk1_1
   CASE DEFAULT
     WRITE(*,*) 'Error in rtopfg value: implemented methods are either 1 or 0'
 END SELECT
 WHERE (suro <= 1.0E-10) suro = 0.0_dp
END SUBROUTINE proute

SUBROUTINE write_ts (outfld_ts, &
                     rain_hru, petinp_hru, &
                     ceps_hru, surs_hru, ifws_hru, uzs_hru, lzs_hru, agws_hru, &
                     supy_hru, suro_hru, ifwo_hru, agwo_hru, pero_hru, igwi_hru, &
                     cepe_hru, uzet_hru, lzet_hru, agwet_hru, baset_hru, taet_hru, &
                     ifwi_hru, uzi_hru, infil_hru, perc_hru, lzi_hru, agwi_hru, &
                     suri_hru, rempet_hru, gwvs_hru, watin_hru, watdif_hru, &
                     pers_hru)

 CHARACTER(len=*), INTENT(IN) ::  outfld_ts
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: rain_hru, petinp_hru                                       !climatological input
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: ceps_hru, surs_hru, ifws_hru, uzs_hru, lzs_hru, agws_hru   !storage series
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: supy_hru, suro_hru, ifwo_hru, agwo_hru, pero_hru, igwi_hru
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: cepe_hru, uzet_hru, lzet_hru, agwet_hru, baset_hru, taet_hru
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: ifwi_hru, uzi_hru, infil_hru, perc_hru, lzi_hru, agwi_hru
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: suri_hru, rempet_hru, gwvs_hru, watin_hru, watdif_hru
 REAL(DP), DIMENSION(:,:), INTENT(IN) :: pers_hru
 INTEGER(I4B) :: nti, nseries, status

 nti = SIZE(rain_hru,1)
 nseries  = SIZE(rain_hru,2)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/rain_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening rain_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(rain_hru)
 CLOSE(15)
 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/petinp_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening petinp_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(petinp_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/ceps_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening ceps_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(ceps_hru)
 CLOSE(15)
 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/surs_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening surs_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(surs_hru)
 CLOSE(15)
 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/ifws_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening ifws_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(ifws_hru)
 CLOSE(15)
 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/uzs_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening uzs_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(uzs_hru)
 CLOSE(15)
 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/lzs_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening lzs_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(lzs_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/agws_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening agws_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(agws_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/supy_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening supy_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(supy_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/suro_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening suro_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(suro_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/ifwo_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening ifwo_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(ifwo_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/agwo_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening agwo_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(agwo_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/pero_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening pero_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(pero_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/igwi_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening igwi_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(igwi_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/cepe_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening cepe_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(cepe_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/uzet_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening uzet_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(uzet_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/lzet_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening lzet_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(lzet_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/agwet_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening agwet_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(agwet_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/baset_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening baset_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(baset_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/taet_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening taet_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(taet_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/ifwi_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening ifwi_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(ifwi_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/uzi_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening uzi_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(uzi_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/infil_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening infil_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(infil_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/perc_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening perc_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(perc_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/lzi_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening lzi_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(lzi_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/agwi_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening agwi_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(agwi_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/suri_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening suri_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(suri_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/rempet_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening rempet_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(rempet_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/gwvs_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening gwvs_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(gwvs_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/watin_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening watin_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(watin_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/watdif_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening watdif_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(watdif_hru)
 CLOSE(15)

 OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/pers_hru_bin.ts', &
       STATUS='REPLACE', ACTION='WRITE', &
       FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
 IF (status /= 0) THEN
  WRITE(*,*) 'error opening pers_hru_bin.ts file for output'
 END IF
 WRITE(15) TRANSPOSE(pers_hru)
 CLOSE(15)
 !alternatively, ASCII array could have been written out as:
 !OPEN (UNIT=15, FILE= outfld_ts(1:LEN_TRIM(outfld_ts)) // '/rain.ts', &
 !      STATUS='REPLACE', ACTION='WRITE', IOSTAT=status)
 !IF (status /= 0) THEN
 ! WRITE(*,*) 'error opening rain.ts file for output'
 !END IF
 !DO i=1,nti
 ! WRITE(15,*) rain_hru(i,:)
 !END DO
 !CLOSE(15)

END SUBROUTINE write_ts

SUBROUTINE drouteSWE (t0, dti, psur, surs, suro, mask)
  ! distributed routing via geoclaw SWE
  USE fixedgrids_module, only : fgrid1_h    ! [nx,ny]

  REAL(DP),                     INTENT(IN)    :: t0      ! initial time for SWE
  REAL(DP),                     INTENT(IN)    :: dti     ! integration time for SWE
  REAL(DP), DIMENSION(:,:),     INTENT(IN)    :: psur    ! water depth [L] added to surs prior to routing
  REAL(DP), DIMENSION(:,:),     INTENT(INOUT) :: surs    !
  REAL(DP), DIMENSION(:,:),     INTENT(OUT)   :: suro
  LOGICAL,  DIMENSION(:,:),     INTENT(IN)    :: mask

  ! local
  INTEGER(I4B) :: status

  WHERE (mask)
    surs = psur           ! update state prior to routing map
  ELSEWHERE
    surs = 0.0d0
  END WHERE

  WHERE(surs < 0.0d0)
     surs = 0.0d0
  END WHERE

  !WRITE(*,*) 'setting surs = 0.01 for debugging:'
  !surs = 0.01d0
  !debugging output
  !WRITE(t0str,'(I10.10)') NINT(t0)
  !surs_fname(11:20) = t0str
  !
  !OPEN(UNIT=15, FILE=surs_fname,  STATUS='unknown', ACTION='WRITE', &
  !     FORM='UNFORMATTED', ACCESS='STREAM', IOSTAT=status)
  !WRITE(15) surs
  !CLOSE(15)
  !WRITE(*,*) 'SIZE(surs) : ',SIZE(surs)
  !WRITE(*,*) 'SHAPE(surs): ',SHAPE(surs)
  !WRITE(*,*) 'SIZE(fgrid1_h) : ',SIZE(fgrid1_h)
  !WRITE(*,*) 'SHAPE(fgrid1_h): ',SHAPE(fgrid1_h)

  suro = surs ! [mm]
  surs = surs / 1000.0D0
  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),'drouteSWE :: surs mass pre-route [m3]:',SUM(surs, MASK=mask) * dx**2
  END IF
  CALL amr2(t0, dti, surs=surs)   !route [SWE] surs -> qinit in read_qinit()
  IF (verbose > 0) THEN
    WRITE(*,*)'drouteSWE:: surs mass pos-route:',SUM(fgrid1_h, MASK=mask) * dx**2
  END IF
  WHERE (mask)
    surs = fgrid1_h * 1000.0d0
    suro = suro - surs
  ELSEWHERE
    suro = 0.0d0
    surs = 0.0d0
  END WHERE

END SUBROUTINE drouteSWE

SUBROUTINE drouteK01 (t0, dti, psur, surs, suro, mask)
  ! distributed routing via 1D kinematic wave cascade

  REAL(DP),                     INTENT(IN)    :: t0      ! real time [from model POW]
  REAL(DP),                     INTENT(IN)    :: dti     ! integration time
  REAL(DP), DIMENSION(:,:),     INTENT(IN)    :: psur    ! water depth [L] added to surs prior to routing
  REAL(DP), DIMENSION(:,:),     INTENT(INOUT) :: surs    ! [L] water depth
  REAL(DP), DIMENSION(:,:),     INTENT(OUT)   :: suro    ! cell-based water surface balance
  LOGICAL,  DIMENSION(:,:),     INTENT(IN)    :: mask

  ! local
  INTEGER(I4B) :: status

  !WRITE(*,*) NEW_LINE('A'),'drouteK01 :: surs mass pre-psur [m3]:',SUM(surs, MASK=mask) * dx**2 / 1000.0D0
  !WRITE(*,*) NEW_LINE('A'),'drouteK01 :: psur mass [m3]         :',SUM(psur, MASK=mask) * dx**2 / 1000.0D0

  WHERE (mask)
    surs = psur          ! update state prior to routing map
  ELSEWHERE
    surs = 0.0d0
  END WHERE

  WHERE(surs < 0.0d0)
    surs = 0.0d0
  END WHERE

  suro = surs ! [mm]
  surs = surs / 1000.0D0                          ! to [m] for routing
  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),'drouteK01 :: surs mass pre-route [m3]:',SUM(surs, MASK=mask) * dx**2
  END IF
  CALL kin1D(t0, dti, mask, surs)
  IF (verbose > 0) THEN
    WRITE(*,*) NEW_LINE('A'),'drouteK01 :: surs mass pos-route [m3]:',SUM(surs, MASK=mask) * dx**2
  END IF
  surs = surs * 1000.0D0                          ! back to [mm]
  WHERE (mask)
    suro = suro - surs                      ! > 0 for outflow from cell
  END WHERE

END SUBROUTINE drouteK01

END MODULE ModuleHspfPwater
