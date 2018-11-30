PROGRAM hspf90
 ! Development: Javier Garcia-Pintado. Last version: 27/04/2015
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
 USE iso_varying_string
 USE ModuleHspfPwater,        only : HSPF_Pwater         ! Just required for the corresponding runs


 IMPLICIT NONE

 CHARACTER(len=PathLength)   :: runfilep        !I run input data file name, possibly with path (e.g. 'runhydro.run')

 CHARACTER(len=StringLength) :: modelfile       !I name of main model input data file (model-specific)

 CHARACTER(len=PathLength) :: modelfilep

 INTEGER(I4B)       :: status
 INTEGER(I4B)       :: verbolev = 1
 CHARACTER(len=PathLength)       :: tmpdir
 LOGICAL ::            fexist
 CHARACTER(len=PathLength)   :: syscmd

 CALL getarg(1, runfilep)
 WRITE(*,*) 'runfilep:',runfilep

 OPEN(UNIT=11, FILE=runfilep, STATUS='OLD', ACTION='READ', &
      IOSTAT=status)
 IF (status == 0) THEN
   READ(11,'(A)') modelfile
   READ(11,'(A)') dsnsim
   READ(11,'(A)') pid
   CLOSE(11)

   IF (verbolev > 0) THEN
     WRITE (*,*) 'input file: ', modelfile
     WRITE (*,*) 'dsnsim: ',     dsnsim
     WRITE(*,*)  'pid:',         pid
   END IF

   modelfilep = trim(dsnsim) // '/input/' // trim(pid) // '/' // trim(modelfile)
  ! WRITE(modelfilep,'(A,A,A,A)') dsnsim(1:LEN_TRIM(dsnsim)),"/input/",pid(1:LEN_TRIM(pid)),modelfile(1:LEN_TRIM(modelfile))
   WRITE (*,*) 'path/model input file: ', modelfilep
 ELSE
   WRITE (*,1002) status
    1002 FORMAT (1X,'open runfile failed--status = ', I6)
    STOP
 END IF

 ! make output directory structure
 tmpdir = trim(dsnsim) // '/output/maps/' // trim(pid)
 INQUIRE(FILE=tmpdir, EXIST=fexist)
 IF (.NOT. fexist) THEN
   syscmd = 'mkdir -p ' // trim(tmpdir)
   !CALL EXECUTE_COMMAND_LINE (syscmd)
   CALL SYSTEM(syscmd, status)
   IF (status /= 0) STOP 'hspf90 --ERR01--' 
 END IF
 tmpdir = trim(dsnsim) // '/output/ts/' // trim(pid)
 INQUIRE(FILE=tmpdir, EXIST=fexist)
 IF (.NOT. fexist) THEN
   syscmd = 'mkdir -p ' // trim(tmpdir)
   !CALL EXECUTE_COMMAND_LINE (syscmd)
   CALL SYSTEM(syscmd, status)
   IF (status /= 0) STOP 'hspf90 --ERR02--'
 END IF


CALL HSPF_Pwater(modelfilep)

END PROGRAM hspf90
