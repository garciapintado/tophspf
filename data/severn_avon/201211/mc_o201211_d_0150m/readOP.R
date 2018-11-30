# foo overpass to demonstrate specific-time distributed output from topHSPF
# this example does not contain input data for assimilation
yssm <- list()
yssm$time    <- seq(RF$staT,RF$endT,by=3600*6)                 # every 6 hours 
#yssm$timeStr <- paste(as.character(RF$seqT),'UTC')
RF$y$uzs$s$data <- yssm; rm(yssm)

