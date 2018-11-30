PROGRAM r2p
 ! test assignation from real to pointers
 ! gfortran -c nrtype.F90 -o nrtype.o
 ! gfortran -c r2p.F90    -o r2p.o
 ! gfortran -o r2p r2p.o nrtype.o
 USE nrtype
 
 REAL(DP), DIMENSION(4), TARGET  :: x
 REAL(DP), DIMENSION(:), POINTER :: p => null()

 x = (/ 2.d0, 5.d0, 2.d0, 1.d0 /)

 ALLOCATE(p(SIZE(x))) 
 p = x
 WRITE(*,*)'associated(p):',associated(p) ! T
 WRITE(*,*)'p:',p

 x(1) = 25.0
 WRITE(*,*)'p:',p

 DEALLOCATE(p) 
 WRITE(*,*)'associated(p):',associated(p) ! F
 p => x ! pointer assignment
 WRITE(*,*)'associated(p):',associated(p) ! T
 WRITE(*,*)'p:',p

 p(1) = 18.0
 WRITE(*,*)'x:',x

END PROGRAM r2p
