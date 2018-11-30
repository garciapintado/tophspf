PROGRAM test_sort
  ! Compile as

  ! gfortran -c ModuleFloodTools.F90
  ! gfortran -c test_sort.F90
  ! gfortran -o test_sort ModuleFloodTools.o test_sort.o


  USE ModuleFloodTools, only : sorti

  IMPLICIT NONE


  INTEGER, DIMENSION(7) :: A
  INTEGER            :: i

  A = (/ 2,5,4,4,1,2,1 /)

  WRITE(*,*) 'A:',            A
  WRITE(*,*) 'MINLOC(A):',    MINLOC(A)      ! returns 5
  WRITE(*,*) 'MINLOC(A(2:7):',MINLOC(A(2:7)) ! returns 4

  CALL sorti(A)
  WRITE(*,*) 'sort(A):',A

END PROGRAM test_sort
