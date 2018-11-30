PROGRAM main
  ! simple launcher to test geoclaw as dependent routing software
  ! this program just initializes a simple domain and everytime step conducts a synthetic addition/abstraction of water [h variable] in the domain

  USE nrtype
  USE amr2_module
  USE fixedgrids_module, only : fgrid1_h
  USE ModuleHspfPwater  ! foo module in this folder just for surs as global module variable

  IMPLICIT NONE

  REAL(DP) :: inctime = 208800.0d0 !total simulation time
  REAL(DP) :: inct    = 3600.0d0   ! outer time step
  REAL(DP) :: deltat  = 900.0d0    ! inner time step
  INTEGER  :: nsteps  = 58         ! 58 hours running time [10*inct]
  INTEGER  :: it,iit,j,i
  REAL(DP) :: stepdiv
  INTEGER  :: innsteps
  REAL(DP) :: runtime = 0.0d0
  REAL(DP) :: t       = 0.0d0    ! model time always initialised to 0.0
  CHARACTER(len=20) :: surs_fname = 'fort.surs_xxxxxxxxxx'
  CHARACTER(len=10) :: out_timeStr
  INTEGER  :: ncols,nrows

  !LOGICAL  :: construct, destruct
  ! CONSTRUCTOR
  ! geoclaw routing model
  call amr2(t, deltat, construct=.TRUE.)         !t,deltat not used. Construct call
  call ConstructSURS

  nrows = SIZE(surs,1)
  ncols = SIZE(surs,2)
  WRITE(*,*) 'SIZE(surs) : ',SIZE(surs)
  WRITE(*,*) 'SHAPE(surs): ',SHAPE(surs)
  WRITE(*,*) 'surs(1,:)  : ',surs(1,:)

  stepdiv = inct/deltat
  innsteps = nsteps * INT(stepdiv)

! RUNNER
!  t = 5000.0d0
!  call amr2(t,deltat)
  !WRITE(*,*) 'fgrid1_h(2:,):'

!  t = 6000.0d0
!  surs = surs * 2.0d0
!  call amr2(t,deltat)
!  WRITE(*,*) 'main:SIZE(fgrid1_h,2)',SIZE(fgrid1_h,2)
!  OPEN(95,file='fort.fg01_h',status='unknown',form='formatted')
!  DO j=1,SIZE(fgrid1_h,2)
!    WRITE (95,*) fgrid1_h(:,j)
!  END DO
!  CLOSE(95)

  WRITE(out_timeStr,'(I10.10)') NINT(t)
  surs_fname(11:20) =  out_timeStr     ! JGP output filenames as timestamps
  OPEN(95,file=surs_fname,status='unknown',form='formatted')
    DO j=1,SIZE(surs,1)
      WRITE (95,*) surs(j,:)
    END DO
  CLOSE(95)

  outer: DO it = 1,nsteps
    WRITE(*,*) 'starting it', it
    inner: DO iit = 1, INT(stepdiv)
       WRITE(*,*) 't:   ', t
      ! scenario a:
      surs = surs ! no hydrologic operation on surs
      !scenario b: simulate some precipitation upstream + infiltration
      surs(1:30,:) = surs(1:30,:) + 100.0E-03 / 3600.0d0 * deltat     ! 100 mm / h rainfall
      forall (i=30:50)
        surs(:,i) = surs(:,i) * 0.5d0                  ! foo infiltration in central columns
      end forall

      ! forecast runoff and update surs
      call amr2(t,deltat)
      forall (j=2:nrows-1,i=2:ncols)
        surs(j,i) = fgrid1_h(j,i)
      end forall

      t = t + deltat
      !output
      WRITE(out_timeStr,'(I10.10)') NINT(t)
      surs_fname(11:20) =  out_timeStr     ! JGP output filenames as timestamps
      OPEN(95,file=surs_fname,status='unknown',form='formatted')
      DO j=1,SIZE(surs,1)
        WRITE (95,*) surs(j,:)
      END DO
  CLOSE(95)

    !  OPEN(95,file='fort.fg01_h',status='unknown',form='formatted')
    !  DO j=1,SIZE(fgrid1_h,2)
    !    WRITE (95,*) fgrid1_h(:,j)
    !  END DO
    !  CLOSE(95)

    END DO inner
  END DO outer

  ! DESTRUCTOR
  call amr2(t,deltat, destruct=.TRUE.)
  !stop

END PROGRAM main
