SUBROUTINE runhydro(runfile, charsfile)
 ! Development: Javier Garcia-Pintado. Last version: 31/04/2010
 !
 ! For compilation as executable see file runhydro.f90, which is an alternative to this SUBROUTINE
 ! Compilation as dynamic library to be loaded into R package hydrosim:
 ! MAKEFLAGS="CFLAGS=-O3" R CMD SHLIB --output=hydrosim.so nrtype.f90 hydrosim_hspf_pwater_mod.f90 runhydro_so.f90
 ! The command will sequentially execute:

 ! gfortran  -fpic -g -O -c  nrtype.f90 -o nrtype.o
 ! gfortran  -fpic -g -O -c  hydrosim_hspf_mod.f90 -o hydrosim_hspf_mod.o
 ! gfortran  -fpic -g -O -c  runhydro_so.f90 -o runhydro_so.o
 ! gfortran -shared -L/usr/local/lib64 -o hydrosim.so nrtype.o hydrosim_hspf_pwater_mod.o runhydro_so.o
 !
 ! valgrind usage:
 ! R -d valgrind --vanilla < hydrosim.Rcheck/hydrosim-Ex.R

 USE nrtype
 USE ModuleHspfPwater,        only : HSPF_Pwater         ! Just required for the corresponding runs


 IMPLICIT NONE

 CHARACTER(len=255) :: runfile         !I run input data file name (e.g. 'runhydro.run')
 INTEGER(I4B) :: charsfile
 CHARACTER(len=charsfile) :: orunfile  !stripped I run input data file name


 CHARACTER(len=25)  :: modelsel        !I model selected to run. Possible values: 'HSPF_PWATER','MARIAM', 'TOPMODEL'
 CHARACTER(len=25)  :: modelfile       !I name of main model input data file (model-specific)
 CHARACTER(len=150) :: datdir_pid      !I path to the watershed folder where 'modelfile' and other direct input resides

 CHARACTER(len=175) :: modelfilep      !M pathed name of main model input data file

 INTEGER(I4B)       :: status          !O

 INTEGER(I4B) :: verbolev = 0          !control output information

 orunfile = ''
 modelfilep = ''

 WRITE(orunfile,'(A)') runfile(1:charsfile)
 IF (verbolev > 0) THEN
   WRITE (*,*) 'orunfile:   ', orunfile
 END IF

 OPEN( UNIT=11, FILE=orunfile, STATUS='OLD', ACTION='READ', &
       IOSTAT=status)
 IF ( status == 0 ) THEN
   READ(11,'(A)') modelsel
   READ(11,'(A)') modelfile
   READ(11,'(A)') datdir_pid
   CLOSE(11)

   WRITE(modelfilep,'(A,A,A)') datdir_pid(1:LEN_TRIM(datdir_pid)),"/",modelfile(1:LEN_TRIM(modelfile))

   IF (verbolev > 0) THEN
     WRITE (*,*) 'selected model:   ', modelsel
     WRITE (*,*) 'model input file: ', modelfile
     WRITE (*,*) 'datdir_pid: ', datdir_pid
     WRITE (*,*) 'pathed model input file: ', modelfilep
   END IF
 ELSE
   WRITE (*,1002) status
    1002 FORMAT (1X,'open runfile failed--status = ', I6)
    STOP
 END IF

SELECT CASE (modelsel)
CASE('HSPF_PWATER')
 IF (verbolev > 0) THEN
   WRITE (*,'(1X,A6,A8,A10/)') 'Model ',modelsel ,' selected!'
 END IF
 CALL HSPF_Pwater(modelfilep)
CASE('topbal')
 IF (verbolev > 0) THEN
   WRITE (*,'(1X,A6,A8,A10/)') 'Model ',modelsel ,' selected!'
 END IF
 WRITE (*,'(A,A)') 'Model not implemented: ',modelsel
 !CALL topbal_main(modelfilep)
CASE('MARIAM')
 IF (verbolev > 0) THEN
   WRITE (*,'(1X,A6,A8,A10/)') 'Model ',modelsel ,' selected!'
 END IF
 WRITE (*,'(A,A)') 'Model not implemented: ',modelsel
 !CALL hydrosim_MARIAM(datdir_pid, modelfilep)
CASE('topmodel')
 IF (verbolev > 0) THEN
   WRITE (*,'(1X,A6,A8,A10/)') 'Model ',modelsel ,' selected!'
 END IF
 WRITE (*,'(A,A)') 'Model not implemented: ',modelsel
CASE DEFAULT
 WRITE (*,'(A,A)') 'Model not implemented: ',modelsel
END SELECT

RETURN

END SUBROUTINE runhydro
