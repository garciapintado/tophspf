# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/clawutil/src/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

EXE = hspf90_geoclaw            # Executable to create
CLAW_PKG = geoclaw              # Clawpack package to use
OUTDIR  = output                # Directory for output
PLOTDIR = plots                 # Directory for plots
RESTART = False
#SETRUN_FILE = setrun.py        # File containing function to make data
#SETPLOT_FILE = setplot.py      # File containing function to set plots
# SETRUN_FILE = setrun.R
# SETPLOT_FILE = setplot.R
#PREPROCESSOR = R               # Created by Javier: 'R' or 'python' to be used by Makefile.common to create the data processing for Geoclaw
                               # unused by now

# Environment variable FC should be set to fortran compiler, e.g. gfortran

# Compiler flags can be specified here or set as an environment variable
FFLAGS ?= -g -O
# The following flag to the linker seems to solve an issue about x86_64 architecture on MAC laptop, which prevented correct linking
LFLAGS = -lstdc++
# ---------------------------------
# List of sources for this program:
# ---------------------------------
HSPF = $(HOME)/docs/hydrology/tophspf

TOPLIB  = $(HSPF)/src
AMRLIB  = $(CLAW)/amrclaw/src/2d
GEOLIB  = $(CLAW)/geoclaw/src/2d/shallow



# amr2_module.f90 is the dependent version of amr.f90

MODULES = \
 $(TOPLIB)/nrtype.F90 \
 $(TOPLIB)/iso_varying_string.F90 \
 $(TOPLIB)/Strings/stringmod_nr.F90 \
 $(GEOLIB)/utility_module.f90 \
 $(AMRLIB)/amr_module.f90 \
 $(AMRLIB)/regions_module.f90 \
 $(AMRLIB)/gauges_module.f90 \
 $(GEOLIB)/geoclaw_module.f90 \
 $(GEOLIB)/topo_module.f90 \
 $(TOPLIB)/ModuleFloodTools.F90 \
 $(GEOLIB)/surs2qinit_module.F90 \
 $(GEOLIB)/refinement_module.f90 \
 $(GEOLIB)/fixedgrids_module.f90 \
 $(GEOLIB)/fgmax_module.f90 \
 $(GEOLIB)/amr2_module.F90 \
 $(TOPLIB)/ModuleHspfPwater.F90 \

# note: fgmax_interpolate0.f90 :: 0-order interpolation instead of bilinear
# note: stepgrid_fgvar.F90     :: JGP time-varying fixed_grid output as opossed to stepgrid.f
SOURCES = \
  $(TOPLIB)/setprob.F90 \
  $(GEOLIB)/qinit.f90 \
  $(GEOLIB)/topo_update.f90 \
  $(GEOLIB)/cellgridintegrate2.f \
  $(GEOLIB)/topointegral.f \
  $(GEOLIB)/bilinearintegral.f \
  $(GEOLIB)/stepgrid_fgvar.F90 \
  $(CLAW)/riemann/src/rpn2_geoclaw.f \
  $(CLAW)/riemann/src/rpt2_geoclaw.f \
  $(CLAW)/riemann/src/geoclaw_riemann_utils.f \
  $(GEOLIB)/src2.f90 \
  $(GEOLIB)/src1d.f90 \
  $(GEOLIB)/step2.f90 \
  $(GEOLIB)/flux2fw.f \
  $(GEOLIB)/qad.f \
  $(GEOLIB)/valout.f \
  $(GEOLIB)/filval.f90 \
  $(GEOLIB)/filpatch.f90 \
  $(AMRLIB)/prefilp.f90 \
  $(AMRLIB)/trimbd.f \
  $(AMRLIB)/bound.f \
  $(AMRLIB)/intfil.f \
  $(GEOLIB)/bc2amr.f \
  $(GEOLIB)/update.f \
  $(GEOLIB)/setaux.f90 \
  $(GEOLIB)/flag2refine2.f90  \
  $(AMRLIB)/flagregions2.f90  \
  $(GEOLIB)/allowflag.f90  \
  $(GEOLIB)/b4step2.f90 \
  $(GEOLIB)/dumpgauge.f \
  $(GEOLIB)/upbnd.f  \
  $(GEOLIB)/tick.f \
  $(GEOLIB)/setgrd.f \
  $(GEOLIB)/gfixup.f \
  $(GEOLIB)/ginit.F90 \
  $(GEOLIB)/getmaxspeed.f90 \
  $(GEOLIB)/advanc.f \
  $(GEOLIB)/fgmax_read.f90 \
  $(GEOLIB)/fgmax_frompatch.f90 \
  $(GEOLIB)/fgmax_interpolate0.f90 \
  $(GEOLIB)/fgmax_values.f90 \
  $(GEOLIB)/fgmax_finalize.f90 \
  $(AMRLIB)/quick_sort1.f \
  $(AMRLIB)/estdt.f \
  $(AMRLIB)/check4nans.f90 \
  $(AMRLIB)/spest2.f  \
  $(AMRLIB)/init_iflags.f \
  $(AMRLIB)/igetsp.f \
  $(AMRLIB)/reclam.f \
  $(AMRLIB)/birect.f \
  $(AMRLIB)/cleanup.f \
  $(AMRLIB)/colate2.f \
  $(AMRLIB)/errest.f \
  $(AMRLIB)/prepregstep.f \
  $(AMRLIB)/prepbigstep.f \
  $(AMRLIB)/bufnst2.f \
  $(AMRLIB)/flagger.f \
  $(AMRLIB)/flagregions.f \
  $(AMRLIB)/errf1.f \
  $(AMRLIB)/fixcapaq.f \
  $(AMRLIB)/flglvl2.f \
  $(AMRLIB)/fluxad.f \
  $(AMRLIB)/fluxsv.f \
  $(AMRLIB)/grdfit2.f \
  $(AMRLIB)/moment.f \
  $(AMRLIB)/nestck2.f \
  $(AMRLIB)/prepc.f \
  $(AMRLIB)/prepf.f \
  $(AMRLIB)/projec2.f \
  $(AMRLIB)/signs.f \
  $(AMRLIB)/findcut.f \
  $(AMRLIB)/smartbis.f \
  $(AMRLIB)/putnod.f \
  $(AMRLIB)/putsp.f \
  $(AMRLIB)/regrid.f \
  $(AMRLIB)/setuse.f \
  $(AMRLIB)/stst1.f \
  $(AMRLIB)/nodget.F90 \
  $(AMRLIB)/basic.f \
  $(AMRLIB)/outval.f \
  $(AMRLIB)/copysol.f \
  $(AMRLIB)/outvar.f \
  $(AMRLIB)/outmsh.f \
  $(AMRLIB)/outtre.f \
  $(AMRLIB)/domain.f  \
  $(AMRLIB)/cellave.f \
  $(AMRLIB)/fdisc.f \
  $(AMRLIB)/fss.f \
  $(AMRLIB)/zeroin.f \
  $(AMRLIB)/setflags.f \
  $(AMRLIB)/shiftset2.f \
  $(AMRLIB)/conck.f \
  $(AMRLIB)/domshrink.f \
  $(AMRLIB)/domprep.f \
  $(AMRLIB)/domup.f \
  $(AMRLIB)/domcopy.f \
  $(AMRLIB)/setdomflags.f \
  $(AMRLIB)/setIndices.f \
  $(AMRLIB)/coarseGridFlagSet.f \
  $(AMRLIB)/addflags.f \
  $(AMRLIB)/baseCheck.f \
  $(AMRLIB)/drivesort.f \
  $(AMRLIB)/flagcheck.f \
  $(AMRLIB)/domgrid.f \
  $(AMRLIB)/setPhysBndryFlags.f \
  $(AMRLIB)/griddomup.f \
  $(AMRLIB)/griddomcopy.f \
  $(AMRLIB)/griddomshrink.f \
  $(AMRLIB)/coarsen.f \
  $(AMRLIB)/auxcoarsen.f \
  $(AMRLIB)/intcopy.f \
  $(AMRLIB)/preintcopy.f \
  $(AMRLIB)/icall.f \
  $(AMRLIB)/preicall.f \
  $(AMRLIB)/inlinelimiter.f \
  $(AMRLIB)/cstore.f \
  $(AMRLIB)/saveqc.f \
  $(AMRLIB)/check.f \
  $(AMRLIB)/restrt.f \
  $(AMRLIB)/opendatafile.f \
  $(AMRLIB)/init_alloc.f90 \
  $(AMRLIB)/restrt_alloc.f90 \
  $(AMRLIB)/resize_alloc.f90 \
  $(TOPLIB)/runhydro.F90

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)

# Construct the topography data
#.PHONY: topo all
#topo:
#	python maketopo.py

#all:
#	$(MAKE) topo
#	$(MAKE) .plots
#	$(MAKE) .htmls

