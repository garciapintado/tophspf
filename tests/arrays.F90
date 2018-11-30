PROGRAM arrays

  IMPLICIT NONE

  INTEGER :: i,j,n
  REAL, DIMENSION(4,5) :: A
  REAL, DIMENSION(20)  :: B,C

  n=20

  DO j=1,5
    A(:,j) = j
  END DO

  DO j=1,20
     B(j) = j+2.0
  END DO

  DO i=1,4
   WRITE(*,*) A(i,:)
  END DO

  WRITE(*,*) ' '
  WRITE(*,*) B

  A = reshape(B,SHAPE(A))

  WRITE(*,*) 'SHAPE(A):',SHAPE(A)

  DO i=1,4
   WRITE(*,*) A(i,:)
  END DO


  C = reshape(A,SHAPE(C))

  WRITE(*,*) C

  C = reshape(transpose(A),(/n/))

  WRITE(*,*) ' '
  WRITE(*,*) C
END PROGRAM arrays
