PROGRAM appirts
 ! test irts simple operations on T_irts objects

 ! gfortran -c nrtype.F90 -o nrtype.o
 ! gfortran -c iso_varying_string.F90 -o iso_varying_string.o
 ! gfortran -c ModuleFloodTools.F90 -o ModuleFloodTools.o
 ! gfortran -c appirts.F90 -o appirts.o
 ! gfortran -o appirts appirts.o ModuleFloodTools.o iso_varying_string.o nrtype.o

 USE nrtype
 USE ModuleFloodTools
 
 CHARACTER(len=PathLength)  :: fname = 'irts.dat'                         ! irts input file
 CHARACTER(len=StringLength), DIMENSION(:), POINTER  :: vnames => null()  ! field names in input irts
 INTEGER               :: unit = 11
 INTEGER               :: status
 TYPE(T_irts), POINTER :: irts => null()
 INTEGER, DIMENSION(:), ALLOCATABLE :: ids
 TYPE(T_irts), POINTER :: irts0 => null() !for versatile sampling
 INTEGER               :: ic

 REAL(DP), DIMENSION(6) :: to1
 REAL(DP)               :: to2

 !REAL(DP), DIMENSION(:), POINTER :: p => null()

 ALLOCATE(ids(3))
 ids = [2,4,3]                                               ! From Fortran 2003, this equals (/2,4,3/)
 WRITE(*,*)'ids:',ids

 ALLOCATE(vnames(3))
 vnames(1) = 'foo_p'
 vnames(2) = 'foo_west'
 vnames(3) = 'foo_north'
 !vnames = (/ 'foo_p','foo_west','foo_north' /)


 ! read and print IRTS
 CALL readIRTS(fname, vnames, irts, unit, status)
 CALL writeIRTS(irts)
 
 ! timse series aproximation
 to1 = [ -5., 3., 5., 3600.,6000.,15000. ]        ! vector output times to INTERFACE approxIRTS
 CALL approxIRTS(irts, to1, irts0)
 CALL writeIRTS(irts0)
 
 to2 = 5.4                                        ! scalar outout times to INTERFACE approxIRTS
 CALL approxIRTS(irts, to2, irts0)
 CALL writeIRTS(irts0)

 ! simple operations
 CALL which(irts%vname == 'foo_west',ids)
 WRITE(*,*)'ids:',ids
 irts%value(:,ids) = irts%value(:,ids) * 2.0
 CALL writeIRTS(irts)
 irts%value(6:7,1:2) = 5.5d0 
 irts%value(:,[1,3]) = irts%value(:,[1,3]) * 2.0
 CALL writeIRTS(irts)

END PROGRAM appirts
