RF$staTStr <- "2012-11-01 00:00:00 UTC"    # Nov-Dec 2012 event
RF$endTStr <- "2012-12-10 23:00:00 UTC"    #

rainfp   <- file.path(dsnsim,'input','rain_pointer.txt')                      # ASCII file with pointers to rainfall files [mm/tstep]
dsnrain  <- file.path(PERM,'rain','EA_gridded_SA','2002_2012_01h_1000m','sri')

# TODO: change PE to real data when available
dsnpe    <- file.path(HOME,'docs/DEMON/data/pe/morecs/2003_2007')      # 2007 in here is use as proxi for 2012
peRsav   <- 'morecs_hourly_2007.Rsav'

qhobf    <- file.path(dsnscn,'gaugeEA.Rsav')   # h & q observations
qobmetf  <- file.path(dsnscn,'outgrlets.txt')  # file with q gauge locations

# ---
RF$staT    <- as.POSIXct(RF$staTStr, origin=originTimeStampTZ, tz='UTC')
RF$endT    <- as.POSIXct(RF$endTStr, origin=originTimeStampTZ, tz='UTC')
RF$seqT    <- seqT    <- seq(RF$staT,RF$endT,by='hour')
RF$seqTStr <- seqTStr <- strptime(RF$seqT,format='%Y-%m-%d %H:%M:%S', tz='UTC')

nto <- length(seqT)                                   # INTEGER. integration outer time steps
dto <- as.double(seqT[2]) - as.double(seqT[1])        # outer timestep (that of model input)      [s]
dti <- 900                                            # inner timestep [15'] (<= dto)             [s]
nti <- nto * dto / dti                                # length of time series for internal topHSPF calculation & return

load(file.path(dsnpe,peRsav))
pe <- pehlst[['135']]; rm(pehlst)
pe$time <- pe$time + as.numeric(ISOdate(2012,1,1,0) - ISOdate(2007,1,1,0), units='secs')
pe$value <- round(pmax(pe$value[,2],0.1E-10),10) # [mm/h] 2 corresponds to potential evapotranspiration according to MORECS

load(qhobf); qob <- gaugeEA           # flow time series for calibration/validation
qobMeta <- read.table(qobmetf, as.is=TRUE); names(qobMeta) <- c('east','north','name')
qobMeta$name[1] <- 'haw_bridge_q'
qid <- match(qobMeta$name, names(qob))
RF$qob <- qob[qid]                             # NULL for NA timeseries

rm(seqT,seqTStr,gaugeEA,qhobf)
