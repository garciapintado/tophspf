PROGRAM test_kin1Dorder
  ! Compile as

  ! gfortran -c nrtype.F90
  ! gfortran -c iso_varying_string.F90
  ! gfortran -c ModuleFloodTools.F90
  ! gfortran -c test_kin1Dorder.F90
  ! gfortran -o test_kin1Dorder ModuleFloodTools.o test_kin1Dorder.o iso_varying_string.o nrtype.o


  USE ModuleFloodTools, only : kin1Dorder

  IMPLICIT NONE

  INTEGER, PARAMETER         :: nx = 4
  INTEGER, PARAMETER         :: ny = 7
  INTEGER,  DIMENSION(nx,ny) :: aflx                   ! SIZE : eastings, northings
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: kinor        ! simulation indexes matrix (2,:)
  REAL, DIMENSION(4,7) :: r
  INTEGER            :: i,j
  INTEGER            :: nd                             ! state vector size, not taking into account 0 values in A (considered as out of the domain)

 ! CALL init_random_seed()
  !nx = SIZE(r,1)
  !ny = SIZE(r,2)
  CALL RANDOM_NUMBER(r)
  aflx = RESHAPE((/(i,i=1,nx*ny) /),SHAPE(aflx))

  WRITE(*,*) 'r:'
  WRITE(*,*) 'nx [eastings]: ',nx
  WRITE(*,*) 'ny:[northings]:',ny

  aflx = INT(r * 30.)
  WHERE(aflx == 0)
     aflx = 1
  END WHERE
  aflx(nx,ny) = 0
  WRITE(*,*) 'simulated accuflux:'
  DO j=1,ny               ! geographical display
  !  WRITE(*,*)  (/(i,i=(j-1)*nx+1,(j-1)*nx+nx )  /),':'
    WRITE(*,*)  aflx(:,j)
  END DO

  nd = COUNT(aflx /= 0)
  ALLOCATE(kinor(nd,2))

  CALL kin1Dorder(aflx, kinor)

  WRITE(*,*)'kinor (i [east] | j [north]):'
  DO i=1,nd
   WRITE(*,*) kinor(i,:)
  END DO
  DEALLOCATE(kinor)
END PROGRAM test_kin1Dorder
