staTStr <- "2012-04-01 00:00:00 UTC"    # Apr-May 15-day time window run
endTStr <- "2012-05-10 23:00:00 UTC"    #

rainfp   <- file.path(dsnsim,'rain_pointer.txt')              # ASCII file with pointers to rainfall files [mm/tstep]
dsnrain  <- file.path(glusland,'rain/EA_gridded_SA/2012/rih')

# TODO: change PE to real data when available
dsnpe    <- '/home/pt902904/docs/DEMON/data/pe/morecs/2003_2007'      # these two do no include 2012 but are used here as a proxy
peRsav   <- 'morecs_hourly_2007.Rsav'
petinpfp <- file.path(dsnsim,'petinp.tss')                            # evapotranspiration file [mm/tstep]:  "

qobsf    <- file.path(dsnsim,'gaugeEA.Rsav')
qobsmetf <- file.path(dsnsim,'outgrlets.txt')

# ---
staT <- as.POSIXct(staTStr, origin=originTimeStampTZ, tz='UTC')
endT <- as.POSIXct(endTStr, origin=originTimeStampTZ, tz='UTC')
seqT <- seq(staT,endT,by='hour')
seqTStr <- paste(as.character(seqT),'UTC')

rainfs  <- paste(substr(seqTStr,1,4),substr(seqTStr,6,7),substr(seqTStr,9,10),
                 substr(seqTStr,12,13),substr(seqTStr,15,16),'.bin',sep='')                 # AAAAMMDDhhmm output gridded QPE file
rainfs  <- file.path(dsnrain,rainfs)                                                        # pathed names to binary rainfall maps

nsteps <- length(seqT)                                                                      # number of timesteps INTEGER
inct   <- as.double(seqT[2]) - as.double(seqT[1])         # outer timestep (that of model input) [s]
deltat <- inct                                            # inner timestep (<= inct)             [s]
innsteps <- nsteps * inct / deltat                        # length of time series for internal topHSPF calculation & return

write(rainfs, file=rainfp)

load(file.path(dsnpe,peRsav))
pe <- pehlst[['135']]; rm(pehlst)
initT2007 <- staT-5*365*24*3600
peinit <- which(time(pe) == initT2007)
pe <- pe[peinit:(peinit+nsteps-1),2] # 2 corresponds to potential evapotranspiration according to MORECS
pe <- data.frame(tStr=paste(as.character(time(pe)),'UTC'),mm=value(pe))
write.table(pe, file=petinpfp, col.names=FALSE, row.names=FALSE)

load(qobsf);  qobs <- gaugeEA; rm(gaugeEA,qobsf)
qobsMeta <- read.table(qobsmetf, as.is=TRUE); names(qobsMeta) <- c('east','north','name')
qobsMeta$name[1] <- 'haw_bridge_q'
qid <- match(qobsMeta$name, names(qobs))
qobs <- qobs[qid]                             # NULL for NA timeseries
