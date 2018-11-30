# foo overpass to demonstrate specific-time distributed output from topHSPF
# this example does not contain input data for assimilation
dsny <- file.path('/users','mci','ffir','wsm')
ytbl <- read.table(file.path(dsny,'summary_nosnow.txt'),skip=3, as.is=TRUE, header=FALSE,
                   col.names=c('asar','AD','prlev','polar','nascat'), fill=TRUE)
ynmf <- dir(file.path(dsny,'swi_rescaled'), pattern='_sc.img$')
ynam <- substr(ynmf,5,12)
ytbl <- ytbl[ytbl$asar %in% ynam,]

adt <- rep('22:30:00',nrow(ytbl))                         # is this time about right?
adt[ytbl$AD=='D'] <- '10:30:00'
ytbl$asarTime <- as.POSIXct(strptime(paste(ytbl$asar,adt), "%Y%m%d %H:%M:%S", tz='GMT'))

RF$y$uzs$s$data <- list()
RF$y$uzs$s$data$time   <- ytbl$asarTime              # [POSIXct] required
RF$y$uzs$s$data$names  <- ytbl$asar                  # required slot
RF$y$uzs$s$data$sensor <- rep('asar',nrow(ytbl))     # required slot
for (i in 1:nrow(ytbl)) { # for testing. Should be better loaded/unloaded as just required
  #yy   <- substr(ytbl$asar[i],1,4)
  ymap <- readGDAL(file.path(dsny,'swi_rescaled',ynmf[i]))
  ydf  <- data.frame(coordinates(ymap),swi=ymap@data[,1])
  RF$y$uzs$s$data[[ynam[i]]] <- ydf[ydf[,3] != 0.0,]   # 0.0 are NA in this data
}

#yssm <- list()
#yssm$time    <- seq(RF$staT,RF$endT,by=3600*6)                 # every 6 hours 
#yssm$timeStr <- paste(as.character(RF$seqT),'GMT')
#RF$y$uzs$s$data <- yssm; rm(yssm)
