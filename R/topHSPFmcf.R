topHSPFmcf <- function(dsn, dsndat, region, event, scn,
                       prjCRS='+init=epsg:27700') {

 # +++ purpose +++
 #   Execute a Monte Carlo run for the topHSPF hydrologic model.
 #   Parameters for MC analysis are obtained from parameter files ('headers') for the watershed & event under study.
 #
 # Record of revisions:
 #   Date        Programmer         Description of change
 #   ---         -----------------  ---------------------
 #   10/10/2015  J. Garcia-Pintado  Original code
 #
 #  CHARACTER, INTENT(IN) :: dsn              path to tophspf HOME database
 #  CHARACTER, INTENT(IN) :: dsndat           path to tophspf PERM database
 #  CHARACTER, INTENT(IN) :: region           folder, within dsn-s,  for the region
 #  CHARACTER, INTENT(IN) :: event            folder, within region, for the event
 #  CHARACTER, INTENT(IN) :: scn              folder, within event,  for the scenario
 #  CHARACTER, INTENT(IN) :: prjCRS           EPGS-mode PROJ4 projection code

 # Notes:
 #   This example is to be run on a typical MPI cluster.
 #   MPI cluster should have been launched normally e.g., through impi or hp-mpi tools before this programme is started

 require(hydrosim)
 require(tseries)
 originTimeStampTZ <- '1970-01-01 00:00:00 UTC'
 
 ## set environment

 dsnscn  <- file.path(dsn,    region, event, scn) # data for analysis (whole path) in HOME volume    [backed-up]
 dsnsim  <- file.path(dsndat, region, event, scn) # data for analysis (whole path) in storage volume [non-backed-up]

 RCLAW <- file.path(CLAW,'Rclaw')
 Rsrcs <- dir(RCLAW, pattern='.R$')
 for (i in 1:length(Rsrcs)) {
   source(file.path(RCLAW,Rsrcs[i]))
 }

 source(file.path(dsnscn,'readRF.R'),        local=TRUE)  # basic environment & forecast information
 source(file.path(dsnscn,'header_mc.R'),     local=TRUE)  # parameters for Monte Carlo analysis & sub-catchment inheritance rules
 source(file.path(dsnscn,'header_event.R'),  local=TRUE)  # parameters for input time series info and pre-processing
 RF$path$dsnscn <- dsnscn
 K01 <- NULL

 RF$inctime <- RF$endT - RF$staT
 RF$path$oufld   <- file.path(dsnsim,'output')
 RF$path$outsfld <- file.path(dsnsim,'output','ts')

 makehpcScript(RF$path$hpcscript, RF$path$data,                  # create High Performance Computing .sh script -- RocksCluster
               RF$region, RF$event, RF$scn,
               'input', 'runhydro.run', RF$modexe, RF$m, RF$np)

 ## read static maps
 maskSP <- readGDAL(file.path(dsnscn,RF$mask))
 maskV    <- as.logical(maskSP@data[,1])                         # logical 1D array [1/0]
 maskV[is.na(maskV)] <- FALSE

 RF$G <- getgmeta6.SP(maskSP, proj=RF$dproj)

 demSP <- readGDAL(file.path(dsnscn,RF$dem))
 proj4string(demSP) <- prjCRS
 demV <- demSP@data[,1]

 hruSP <- readGDAL(file.path(dsnscn,RF$hrus))
 hruV <- as.integer(hruSP@data[,1])

 #rm(maskSP,demSP,hruSP)

 # create gauges for model sampling (SpatialPointsDataFrame class)
 RF$gauges <- read.table(file.path(dsnscn,RF$gaugesf), header=FALSE, as.is=TRUE)
 names(RF$gauges) <- c('east','north','name')
 coordinates(RF$gauges) <- c('east','north')
 proj4string(RF$gauges) <- prjCRS
 RF$gauges@data['z'] <- over(RF$gauges,demSP)
 if (1 > 2) {
   library(colorRamps)
   plotGmeta(layer=hruSP, col=matlab.like(9))
   points(RF$gauges,col='darkred',pch=3)
   dsnGBriv <- file.path(HOME,'Documents/hydrology/shareGeo')
   ogrInfo(dsn='.',layer='gb_rivers')
   GBriv   <- readOGR(dsn=dsnGBriv, layer='gb_rivers')
   plot(GBriv, col='navyblue',lwd=1, add=TRUE)
 }
 if (RF$simspdist == 'spdis') {
   if (RF$rwave == 'k01') {
     cat('loading 1D kinematic wave parameters\n')
     source(file.path(dsnscn,'readK01.R'),       local=TRUE)  # for 1D kinematic routing
   } else if (RF$rwave == 'swe') {
     cat('loading clawpack-SWE parameters\n')
         source(file.path(dsnscn,'readGC.R'),        local=TRUE)  # for clawpack-swe routing
     dir.create(file.path(dsnsim,'input','clawpack'))
     writeGC(GC$rundata, dsn=file.path(dsnsim,'input','clawpack'))
     writeGDAL(demSP, fname=GC$demf@fname, drivername='AAIGrid',
             options=c('NODATA=-9999.00','DECIMAL_PRECISION=4'))   # ArcGIS
   }
 }

 # init parlst, which will later regenerated for each ensemble member
 npar <- length(GP)
 parlst <- GP
 for (i in 1:npar) {
   tmpexpr <- eval(parse(file="",text=paste("GP$",names(parlst)[i],sep="")))
   if (class(eval(tmpexpr)) == "numeric") {
     parlst[[i]]  <- eval(tmpexpr)
   } else {                  # "character"
     parlst[[i]]  <- readRast(RF$G, eval(tmpexpr),   dsn=dsnscn, nhrus=nhrus)
     if (class(parlst[[i]]) == 'SpatialGridDataFrame')
       parlst[[i]] <- parlst[[i]]@data[,1]
   }
 }

 # load observation data [for future DA]
 source(file.path(dsnscn,'readOP.R'), local=TRUE)
 # load assimilation protocol (starttime & launching model times) into RF
 source(file.path(dsnscn,'readRuntimes.R'), local=TRUE)

 RF$timenow <- RF$runtimes[1]
 
 ## Monte Carlo sampling & initial conditions
 # global parameters : always generated. Later overwritten either by prior parameters or newly generated by Monte Carlo
 MCgpar <- matrix(NA, nrow=npar, ncol=RF$m)   # standard DA ensemble deployment
 #rownames(MCgpar) <- c("runpid",names(GP)); MCgpar["runpid",] <- 1:maxruns
 rownames(MCgpar) <- names(GP)

 for (ipar in 1:npar) {
   if (MCflags[ipar]) {
     MCgpar[ipar,] <- runif(RF$m, min(GP[[ipar]]),max(GP[[ipar]]))               # random sampling ~U(min,max)
   } else {
    if (is.numeric(GP[[ipar]][1]))
      MCgpar[ipar,] <- rep(GP[[ipar]][1], RF$m)
   }                                                                          # inequality constraints:
   if (names(GP)[ipar] == 'ceps') {                                           # ceps <= cepsc
     MCgpar[ipar,] <- pmin(MCgpar[ipar,], MCgpar['cepsc',])
   }
   if (names(GP)[ipar] == 'uzs') {                                            # uzs <= uzsn
     MCgpar[ipar,] <- pmin(MCgpar[ipar,], MCgpar['uzsn',])
   }
   if (names(GP)[ipar] == 'lzs') {                                            # lzs <= lzsn
     MCgpar[ipar,] <- pmin(MCgpar[ipar,], MCgpar['lzsn',])
   }
 }

 ## Monte Carlo parameter for K01 wave
 if (RF$simspdist == 'spdis' && RF$rwave == 'k01') {
     K01$par <- list()
     K01$par$SGC <- genSGCgrpar(file.path(dsnscn,K01$SGCmc), RF$m,
                                RF$ana$SGCfr$min)                  # [ngr,7,m] randomise friction coefficients
     K01$par$dsl <- pmax(matrix(rnorm(RF$m, mean=K01$dsl[1], sd=K01$dsl[2]),1,RF$m),   # Monte Carlo: FREE downstream slope [,m]
                         RF$ana$dsl$min)
 }

 # subwatershed-specific parameters
 if (RF$mcmode == 0) {                                                           # mcmode == 0 :: inherit parameters from previous run
   if (!file.exists(RF$inh$prior))
     stop('topHSPFmcf:: prior simulation not found')
   load(file.path(RF$inh$prior,'output','MCwpar.Rsav'))                          # prior MCwpar
   load(file.path(RF$inh$prior,'output','glike.Rsav'))                           # prior glike
   MCwpar  <- selectPriorMC(MCwpar, glike, m, RF$inh$stat, RF$inh$decreasing)
 } else {                                                                        # mcmode == 1 :: new Monte Carlo generation
   MCwpar <- vector('list', nhrus)
   if (length(WP) > 0) {
     for (iw in 1:nhrus) {
       MCwpar[[iw]] <- matrix(NA, nrow=npar, ncol=RF$m)
       rownames(MCwpar[[iw]]) <- names(GP)
       for (ipar in 1:npar) {
         if (MCflags[ipar]) {
           #MCwpar[[iw]][ipar,] <- runif(RF$m,min(WP[[iw]][[ipar]]),max(WP[[iw]][[ipar]]))             # random sampling ~U(min,max)
           MCwpar[[iw]][ipar,] <- MCgpar[ipar,] # simple inheritance from Global parameters
         } else {
           if (is.numeric(WP[[iw]][[ipar]][1]))
             MCwpar[[iw]][ipar,] <- rep(WP[[iw]][[ipar]][1], RF$m)
         }                                                                           # inequality constraints:
         if (names(GP)[ipar] == 'ceps') {                                            # ceps <= cepsc
           MCwpar[[iw]][ipar,] <- pmin(MCwpar[[iw]][ipar,], MCwpar[[iw]]['cepsc',])
         }
         if (names(GP)[ipar] == 'uzs') {                                             # uzs <= uzsn
           MCwpar[[iw]][ipar,] <- pmin(MCwpar[[iw]][ipar,], MCwpar[[iw]]['uzsn',])
         }
         if (names(GP)[ipar] == 'lzs') {                                             # lzs <= lzsn
           MCwpar[[iw]][ipar,] <- pmin(MCwpar[[iw]][ipar,], MCwpar[[iw]]['lzsn',])
         }
       }
     }
   }
 } # endif mcmode

 # model files :: write static distributed maps
 topHSPFwriteStatic(RF, K01)

 nt <- length(RF$runtimes)
 ## start looping
 for (it in nt) {                                                     # forecast step

   itStr <- formatC(it, width=4, flag='0')
   RF$it <- it                                                       # just required for debug mode in analyse
   RF$timenow <- RF$runtimes[it]
   if (it < length(RF$runtimes)) {
     RF$timenex <- RF$runtimes[it+1]
     RF$freerun <- FALSE
   } else {
     RF$timenex <- RF$endT
     RF$freerun <- TRUE
   }
   seqfT <- seq(RF$timenow, RF$timenex, by=dto)     # model input | this forecast time sequence
   if (max(seqfT) < RF$timenex)
     seqfT <- c(seqfT,max(seqfT)+dto)              # assure it covers next interruption point (timenex)
   ntof <- length(seqfT)
   ntif <- as.integer(ntof*dto/dti)
   
   # write [links to] dynamic model input
   seqfTStr <- strptime(RF$seqT,format='%Y-%m-%d %H:%M:%S', tz='UTC')
   rainfs  <- paste(substr(seqfTStr,1,4),substr(seqfTStr,6,7),substr(seqfTStr,9,10),
                 substr(seqfTStr,12,13),substr(seqfTStr,15,16),'.bin',sep='')                 # AAAAMMDDhhmm output gridded QPE file
   rainfs  <- file.path(dsnrain,rainfs)                                                        # pathed names to binary rainfall maps
   write(rainfs, file=rainfp)
 
   petinp <- pe$value[pe$time >= min(seqfT) & pe$time <= max(seqfT)]
   if (length(petinp) != ntof)
     stop("petinp tss length != ntof")

   ## serial writing of model files
   cat('writing topHSPF input files\n')
   parnames <- names(GP)
   MCerr <- NULL
   MCres <- list()
   
   for (im in 1:RF$m) {
     for (ipar in 1:npar) {
       if (!is.na(MCgpar[ipar,im]))
         parlst[[ipar]] <-  as.numeric(MCgpar[ipar,im])
       if (length(WP) > 0) {
         for (iw in 1:max(nhrus)) {
           if (!is.na(MCwpar[[iw]][ipar,im]))
             parlst[[ipar]][iw] <- MCwpar[[iw]][ipar,im]
         }
       }
     }
     MCerr[im]  <- writeTopHSPF(pid = im,
                              dto=dto, dti=dti, nto=ntof,
                              hruV     = hruV,
                              rain     = rainfs,
                              petinp   = petinp,
                              parlst   = parlst,
                              RF, K01)
   } # end for im

   ## parallel model run
   #cat('Running array topHSPF job in RocksCluster...\n')

   if (RF$hpc == 'mac') {
     runtime <- macRunTHjob(RF)
   } else {
     err <- rocksRunTHjob(RF)
     # monitoring ::
     # qstat
   }

   ## serial reading of model output files
   simseqT <- seq(min(seqfT),max(seqfT)+dto,by=dti)[-1] # model output [initial time not included]
   
   cat('reading topHSPF output\n')
   for (im in 1:RF$m) { # times match simseqT
     cat('reading im:',im,'\n')
     MCres[[im]] <- readTStopHSPF(dsn  = file.path(RF$path$outsfld, formatC(im,width=5,flag="0",format="d")),
                                  dto = dto, dti = dti, nto = ntof, nhrus = nhrus)
     # mean dti runoff [m3/s] components and total
     MCres[[im]]$suro.m3s <- MCres[[im]]$suro.mm *
                             matrix(rep(c(swareas,sum(swareas)),each=nti),ncol=nhrus+1) / 1000 / dti
     #MCres[[im]]$suro.m3s <- apply(MCres[[im]]$suro.m3s,2,acc2ins)
     MCres[[im]]$ifwo.m3s <- MCres[[im]]$ifwo.mm *
                             matrix(rep(c(swareas,sum(swareas)),each=nti),ncol=nhrus+1) / 1000 / dti
     #MCres[[im]]$ifwo.m3s <- apply(MCres[[im]]$ifwo.m3s,2,acc2ins)
     MCres[[im]]$agwo.m3s <- MCres[[im]]$agwo.mm *
                             matrix(rep(c(swareas,sum(swareas)),each=nti),ncol=nhrus+1) / 1000 / dti
     #MCres[[im]]$agwo.m3s <- apply(MCres[[im]]$agwo.m3s,2,acc2ins)
     MCres[[im]]$rnof.m3s <- MCres[[im]]$suro.m3s  + MCres[[im]]$ifwo.m3s + MCres[[im]]$agwo.m3s
   }

   # calculate generalized likelihood statistics
   stats.nm <- c("runpid","ME","FLOWREL","MAE","RMSE","E","FV","lag.peak","PTLAG",
                 "dif.peak.abs","PDIFF","dif.peak","lambda","HMLE","trust.HMLE")
   glike <- list() # stack of Generalized Likelihood statistics matrices. One [m x nstats] matrix / HRU

   # RF$qob ~> RF$y$q$g
   for (iw in 1:nhrus) {
     if (!is.null(RF$qob[[iw]])) {
       glike[[iw]] <- matrix(NA, nrow=length(stats.nm), ncol=RF$m)
       rownames(glike[[iw]]) <- stats.nm
       for (im in 1:RF$m) {
         qsim <- irts(simseqT,MCres[[im]]$rnof.m3s[,iw])
         glike[[iw]][,im] <- perfstatf(im, qobs = RF$qob[[iw]]$irtsmo, qsim = qsim, dT = dti)
       }
     }
   }
   names(glike) <- names(RF$qob)[1:length(glike)]

   # scatterplots with performance statistics
   glikpar <- list()
   glikpar$wpar <-  MCwpar                # list of matrices
   glikpar$SGCfr <- K01$par$SGC[,'n',]    # augment list with global kinematic wave parameters
   glikpar$dsl   <- K01$par$dsl
   plotGlike(glikpar, glike, stats=stats.nm[-1], MCflags, K01$par,
             dsnplt = file.path(dsnsim,'output','plots'))

   # example simulations: best statistics
   dir.create(file.path(dsnsim,'output','plots','ts'))
   for (iw in 1:nhrus) {
     Emaxid <- which.max(glike[[iw]]['E',])  # Emaxid <- which.max(MCgpar['agwr',])
     topHSPFplotTS(MCres[[Emaxid]], simseqT,   # plot balances for one ensemble member as a simple indication
                   dsnplt = file.path(dsnsim,'output','plots','ts',paste(names(glike)[iw],'_calib',sep='')),
                   plothrus=TRUE, qobs=RF$qob, hrunames=names(glike))
   }
   # example simulations: specific member [e.g. im=1]
     im <- 1
     topHSPFplotTS(MCres[[im]], simseqT,   # plot balances for one ensemble member as a simple indication
                   dsnplt = file.path(dsnsim,'output','plots','ts',paste('im',im,sep='_')),
                   plothrus=TRUE, qobs=RF$qob, hrunames=names(glike))


   # map runoff into a list readable as boundary conditions by the flooding DA software
   qsim <- ts2DA(simseqT, MCres, glike=glike, mo=NULL)
   topHSPFplotTSrunoff(qsim, RF$qob,
                       dsnplt = file.path(dsnsim,'output','plots','ts'))

   save(RF,     file=file.path(RF$path$oufld, paste('RF_',itStr,'.Rsav',sep='')))
   save(MCgpar, file=file.path(RF$path$oufld, paste('MCgpar',itStr,'.Rsav',sep='')))
   save(MCwpar, file=file.path(RF$path$oufld, paste('MCwpar',itStr,'.Rsav',sep='')))
   save(glike,  file=file.path(RF$path$oufld, paste('glike_',itStr,'.Rsav',sep='')))
   save(MCres,  file=file.path(RF$path$oufld, paste('MCres',itStr,'.Rsav',sep='')))
   save(swareas,file=file.path(RF$path$oufld, paste('swareas_',itStr,'.Rsav',sep='')))
   save(K01,    file=file.path(RF$path$oufld, paste('K01_',itStr,'.Rsav',sep='')))
 
   # save runoff timeseries
   dsnsav <- file.path(RF$path$oufld, 'ts', 'R')
   if (!file.exists(dsnsav))
     dir.create(dsnsav)
   save(qsim,     file=file.path(dsnsav, 'qsim.Rsav'))
   #dsnsav <- file.path(RF$path$oufld, 'ts', 'matlab')
   #irts2matlab(qsim, dsnsav, field='irts', onames=names(qsim))
 } # end for it
# return(0)
#} # end function tophspfMCf

                                        # tmp debugging code
#if (3 > 2) {
 library(colorRamps)

 dsnGBriv <- file.path(PERM,'gis','shareGeo','GB_Rivers')
 GBriv <- readOGR(dsn=dsnGBriv,layer='ls_rivers')

 im <- which.max(glike[[7]]['E',])                 # best at Evesham [river Avon]
 imStr <- formatC(im,width=5,flag="0",format="d")

 pltvars <- c('surs','ifws','uzs','lzs','agws')
 zlims <- matrix(c(0,70,
                   0,100,
                   0,10,
                   0,30,
                   0,250),ncol=2,byrow=TRUE)
 for (iv in 1:length(pltvars)) {
   pltvar <- pltvars[iv]
   zlim   <- zlims[iv,]
   
   varfiles <- dir(file.path(dsnsim,'output','maps',imStr),pattern=paste('fort',pltvar,sep='.'))
   tou <- as.POSIXct(as.numeric(substr(varfiles,11,20)), origin=originTimeStampTZ,
                     tz='UTC') # POSIX double representation
   touStr <- strptime(tou,format='%Y-%m-%d %H:%M:%S', tz='UTC')
   varSGDF <- demSP
   #plotGmeta(layer=demSP, col=topo.colors(100))
   #for (i in 1:length(varfiles)) {

   dsnplt <- file.path(dsnsim,'output','plots','maps',pltvar)
   if (!file.exists(dsnplt))
     dir.create(dsnplt, recursive=TRUE)

   for (i in 100:130) {
     varSGDF@data[,1] <- readBin(file.path(dsnsim,'output','maps',imStr,varfiles[i]), what='double',n=RF$G$cols*RF$G$rows)
     #plotGmeta(layer=varSGDF, col=matlab.like(100), add=TRUE, zlim=zlim)
     #text(x=RF$G$e,y=RF$G$n,label=touStr[i],col='white', adj=c(1.1,1.2))
     #Sys.sleep(2.0)
     fname  <- paste(file.path(dsnplt,substr(varfiles[i],6,20))) # .png extension not required
     printMap(varSGDF, fnamep=fname, main = paste(pltvar,touStr[i],sep=' | '),
              zlim=zlim, detlim=0.0, png.plot=TRUE, gmt.border=FALSE,
              scalar=4, colpal=matlab.like(100), colpow=1,
              leglab=paste(pltvar,'[mm]',sep=''), SLDF=GBriv)
   }
 } # end for iv plot
  return(0)
} # end function tophspfMCf

if (1 > 2) {
 qinit <- DEMsp
 ns <- RF$G$cols * RF$G$rows
 qinit@data[,1] <- readBin(file.path(dsnsim,'fort.qinit'), what='double',n=ns)
 plotGmeta(layer=qinit, col=matlab.like(100))
 digitGmeta(layer=qinit, type='points')

 supy <- DEMsp
 supy@data[,1] <- readBin(file.path(dsnsim,'supy.tmp'), what='double',n=ns)
 plotGmeta(layer=supy, col=matlab.like(100))
 digitGmeta(layer=supy, type='points')

 msupy <- DEMsp
 msupy@data[,1] <- readBin(file.path(dsnsim,'msupy.tmp'), what='double',n=ns)
 plotGmeta(layer=msupy, col=matlab.like(100))
 digitGmeta(layer=msupy, type='points')
}

