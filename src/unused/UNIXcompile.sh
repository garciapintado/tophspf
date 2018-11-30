#! /bin/sh
#gfortran -c nrtype.f90
# parallel run
#gfortran -fopenmp -c disthspf_mod.f90
#gfortran -fopenmp -o hydrotools.exe hydrotools.f90 disthspf_mod.o

# non-parallel (standard serial) run
#gfortran -c ohydrosim_hspf_mod.f90
#gfortran -c hydrotools.f90
#gfortran -o hydrotools nrtype.o ohydrosim_hspf_mod.o hydrotools.o
gfortran  -fpic -g -O -c  -Wall nrtype.f90 -o nrtype.o
gfortran  -fpic -g -O -c  -Wall ohydrosim_hspf_mod.f90 -o hydrosim_hspf_mod.o
gfortran  -shared -o libdisthspf.so nrtype.o hydrosim_hspf_mod.o
mv libdisthspf.so ~/lib
gfortran -L/home/pt902904/lib -o runhydro runhydro.f90 -ldisthspf
mv runhydro ~/bin

