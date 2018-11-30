PROGRAM test_minloc
  ! Compile as
  ! gfortran -o test_minloc test_minloc.F90

  IMPLICIT NONE

  REAL, DIMENSION(7) :: A
  INTEGER            :: i

  A = (/ 2.,5.,4.,4.,1.,2.,1. /)

  WRITE(*,*) 'MINLOC(A):',    MINLOC(A)      ! returns 5
  WRITE(*,*) 'MINLOC(A(2:7):',MINLOC(A(2:7)) ! returns 4

END PROGRAM test_minloc
