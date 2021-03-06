# +++ R interface to write files for a single tophspf F90 model run +++

writeTopHSPF <- function(pid   = 1,
                    dto=NULL, dti=NULL, nto=NULL,
                    hruV    = NULL,
                    rain    = NULL,
                    petinp  = NULL,
                    parlst  = NULL,
                    RF, K01)
{
 # INTEGER :: pid                          !unique identifier for parallel simulations (MonteCarlo analysis...)
 #                                         !each call to the function should have a distinct runpid to avoid srambbled results
 # REAL ::    dx                           !cell resolution [m]
 # INTEGER :: nrows, ncols                 !region nrows and ncols
 # REAL ::    dto                          !outer time step [s] - match rainfall & evpt input
 # REAL ::    dti                          !inner time step [s]
 # INTEGER :: nto                       !number of outer timesteps
 # CHARACTER(len=*) simspdist              !'smdis' or 'spdis' spatial resolution for model run
 # CHARACTER(len=*) simspdist              !'k01' [1D kinematic routing cascade] or 'swe' [clawpack-SWE]
 # LOGICAL :: mask                         !integer (1,0) map. Being 1 for values within watershed
 # INTEGER :: hruV                         !hydrological response units map (region to where time series are aggregated, disregarding simspdist.  0 = out of watershed; 1,...,n pixels withineach n hru
 # REAL ::    dem                          !digital elevation model [m]
 # CHARACTER/REAL :: rain                  !rain, may be an 1D array of CHARACTER variables pointing (whole path) to a stack of 2D binary maps [mm/tstep],
 #                                          or a 1D vector of length nto [mm/tstep]
 # CHARACTER/REAL :: petinp                !potential evapotranspiration, either an 1D array of CHARACTER variables pointing (whole path) to a stack of 2D binary maps,
 #                                          or a 1D vector of length nto [mm/tstep]
 # STRUCT    :: parlst                     !this is a tree containing all parameters, with a structure as defined in the documentation of the package
 # CHARACTER :: datdirp                    !absolute path to datafile parent folder
 # CHARACTER :: region                     !subfolder of outfldp where simulation watershed output data will be put
 # CHARACTER :: event                      !subfolder of region where simulation watershed output data will be put

 mainfile  <- "runhydro.run"
 modelfile <- "hspf_pwater.run"
 dsnsim <- file.path(RF$path$data,RF$region,RF$event,RF$scn)

 oufldp  <- file.path(dsnsim, "output");   if (!file.exists(oufldp))       dir.create(oufldp)
 oufldp.ts    <- file.path(oufldp,"ts");   if (!file.exists(oufldp.ts))    dir.create(oufldp.ts)      # time series
 oufldp.maps  <- file.path(oufldp,"maps"); if (!file.exists(oufldp.maps))  dir.create(oufldp.maps)    # map stacks
 oufldp.stats <- file.path(oufldp,"stats");if (!file.exists(oufldp.stats)) dir.create(oufldp.stats)   # performance statistics

 if(pid > 99999) stop("pid should be <= 99999")
 pidfmt <- formatC(pid,width=5,flag="0",format="d")
 #ipid <- 0; outfldp.pid <- ""
 #while (ipid == 0 || file.exists(outfldp.pid)) {
 #  ipid <- ipid + 1
 #  rpid <- trunc(runif(1,min=1,max=1.0E8))
 #  pidfmt.pid <- paste(pidfmt,rpid,sep="_")
 #  outfldp.pid <- file.path(outfldp,pidfmt.pid)
 #}; dir.create(outfldp.pid); rm(ipid)
 infldp.pid    <- file.path(dsnsim,'input',pidfmt)
   dir.create(infldp.pid, showWarnings=FALSE)
 oufldp.ts.pid <- file.path(oufldp.ts,pidfmt)                         # scratch output folder to store time series (F90 code)
 if (file.exists(oufldp.ts.pid)) {
  system(paste("rm -rf",oufldp.ts.pid))}
 dir.create(oufldp.ts.pid)

 write(c(modelfile, dsnsim, pidfmt),
       file=file.path(infldp.pid,mainfile))

 nhrus <- max(hruV) # from parent environment

 dx    <- RF$G$nsres
 ncols <- RF$G$cols
 nrows <- RF$G$rows

 otimes  <- RF$y$uzs$s$data$time
 OP4this <- otimes >  RF$timenow &
            otimes <= RF$timenex    # index of overpass files within the following time window simulation
 if (sum(OP4this) > 0)
   otimes <- otimes[OP4this]
 else
   otimes <- RF$timenex
 
 # data spatio-temporal distribution: |'lco': lumped constant|'lts': lumped timeseries|'dco': distributed constant|'dts': distributed ts|
 rain.std   <- ifelse (class(rain) == 'character','dts',
                     ifelse(length(rain) == ncols*nrows,'dc',
                            ifelse(length(rain) == nto,'lts','lc')))

 petinp.std <- ifelse (class(petinp) == 'character','dts',
                     ifelse(length(petinp) == ncols*nrows,'dc',
                            ifelse(length(petinp) == nto,'lts','lc')))

 # parfiles input may be introduced on a semidistributed hru basis.
 parfiles.std <- NULL
 for ( i in 1:length(parlst)){
   parfiles.std[i] <- ifelse (class(parlst[[i]]) == 'character','dts',                         # distributed time series
                              ifelse(length(parlst[[i]]) == ncols*nrows,'dc',                  # distributed constant
                                     ifelse(length(parlst[[i]]) == nhrus,'sc','lc')))          # semidistributed constant, lumped constant
   if (!(parfiles.std[i] %in% c('dts','dc','sc','lc')))
     stop('writeTopHSPF:: ERR 001 --') 
 }

   fpars     <-  paste("f_",names(parlst),".map",sep="")

   mainfilep   <- file.path(infldp.pid,mainfile)
   modelfilep  <- file.path(infldp.pid,modelfile)             # path/files to model
   fgeometap   <- file.path(infldp.pid, 'f_geometa.asc')

   fmaskp      <- file.path(infldp.pid, 'mask.asc')
   fhrusp      <- file.path(infldp.pid, 'hrus.asc')
   fdemp       <- file.path(infldp.pid, 'dem.bin')
   fk01p       <- file.path(infldp.pid, 'f_k01.asc')
   frainp      <- file.path(infldp.pid, 'f_rain.asc')
   fpetinpp    <- file.path(infldp.pid, 'f_petinp.asc')
   fparsp      <- file.path(infldp.pid, fpars)
   foptimep    <- file.path(infldp.pid, 'opts.asc')
 
   #WRITE fgeometap
   zo <- file(fgeometap,'w')
   cat('# dx ncols nrows t0 dto dti nto simspdist rwave\n', file=zo)
   cat(dx, ncols, nrows, RF$timenow, dto, dti, nto, RF$simspdist, RF$rwave,'\n', file=zo)
   close(zo)

   #link mask map
   system(paste('ln -sf', file.path(dsnsim,'input','mask.asc'), infldp.pid))

   # link HRU map
   system(paste('ln -sf', file.path(dsnsim,'input','hrus.asc'), infldp.pid))

   # link DTM map
   system(paste('ln -sf', file.path(dsnsim,'input','dem.bin'), infldp.pid))

   # if kinematic wave
   if (RF$simspdist == 'spdis' && RF$rwave == 'k01') {
     system(paste('ln -sf', file.path(dsnsim,'input','ldd.bin'), infldp.pid))
     system(paste('ln -sf', file.path(dsnsim,'input','aflx.bin'), infldp.pid))
     system(paste('ln -sf', file.path(dsnsim,'input','cmsk.bin'), infldp.pid))
     system(paste('ln -sf', file.path(dsnsim,'input','cwid.bin'), infldp.pid))
     system(paste('ln -sf', file.path(dsnsim,'input','cdep.bin'), infldp.pid))
     system(paste('ln -sf', file.path(dsnsim,'input','G.data'), infldp.pid))

     zo <- file(fk01p,'w')
     cat('dem.bin',       '\n', file=zo)
     cat('ldd.bin',       '\n', file=zo)
     cat('aflx.bin',      '\n', file=zo)
     cat('cmsk.bin',      '\n', file=zo)
     cat('cwid.bin',      '\n', file=zo)
     cat('cdep.bin',      '\n', file=zo)
     cat('SGCpar.asc',    '\n', file=zo)
     cat(K01$par$dsl[pid],'\n', file=zo)
     close(zo)

     #zo <- file('SGCpar.asc','w')
     write.table(K01$par$SGC[,,pid], file=file.path(infldp.pid,'SGCpar.asc'),
                 row.names=FALSE, quote=FALSE)


     #close(zo)

  }

   #WRITE frainp   (always ASCII, but may be either CHARACTERs pointing toward binary 2Darrays, or univariate time series)
   write(rain,
         file=frainp)

   #WRITE fpetinpp (always ASCII, but may be either CHARACTERs pointing toward binary 2Darrays, or univariate time series)
   write(petinp,
         file=fpetinpp)

   #WRITE fpars INPUT to model. Currently non time-varying (just 'lc' or 'dc'). However hrus' based ASCII input ('smc' semidistributed constant) is possible, and will be
   #forwarded to the fortran model as a distributed 'dc' array
   for (i in 1:length(parlst)) {
     if (parfiles.std[i] == 'lc') {
       write(parlst[[i]],
             file=file.path(infldp.pid,fpars[i]))
     } else if (parfiles.std[i] == 'dc') {
       writeBin(parlst[[i]],
                file.path(infldp.pid,fpars[i]))
     } else if (parfiles.std[i] == 'sc') {
       parfiles.std[i] <- 'dc'
       parhru <- parlst[[i]]
       pardc  <- as.numeric(hruV)
       for ( j in 1:nhrus) {
         pardc[pardc==j] <- parhru[j]
       }
       writeBin(pardc,
                file.path(infldp.pid,fpars[i]))
     } else {
       stop('hspf_water:: unknown input style')
     }
   }

  # overpass file: time vector [model space] with distributed output
  zo <- file(foptimep,'w')
  #if (length(otimes) == 0) {
    #otimes <- difftime(RF$timenex, RF$timenow,, units='secs') # map into model space
  #} else {
    #otimes <- difftime(otimes, RF$timenow, units='secs')
  #}
  write(otimes, file=zo, ncolumns=1, append=TRUE)
  close(zo)

  #WRITE modelfilep describing subfiles
  zo <- file(file.path(infldp.pid,modelfile), 'w')
  cat(fgeometap,           '\n',file=zo)
  cat(fmaskp,              '\n',file=zo)
  cat(fhrusp,              '\n',file=zo)
  cat(fk01p,               '\n',file=zo)
  cat(frainp,   rain.std,  '\n',file=zo)
  cat(fpetinpp, petinp.std,'\n',file=zo)
  cat(rbind(fparsp,' ',parfiles.std,'\n'),sep='',file=zo)
  cat(oufldp.ts.pid,       '\n',file=zo)
  cat(foptimep,            '\n',file=zo)
  close(zo)

  return(0)
}
