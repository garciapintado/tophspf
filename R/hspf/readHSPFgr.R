readHSPFgr <- function(RF, vtypes=NULL) {
  # read gridded output from topHSPF
  #
  # output is a list with named elements, each with two S3 slots:
  # $time  :: POSIXct vector of times
  # $stack :: list with length(time) [ng,m] matrices

  originTimeStampTZ <- '1970-01-01 00:00:00 UTC'

  n <- RF$G$cells
  m <- RF$m
 
  xtypes  <- names(RF$ana)[names(RF$ana) %in% vtypes]
  #xupdate <- sapply(RF$ana[xtypes], function(x){x$u})
  
  dsno  <- file.path(RF$path$oufld,'maps','00001')
  fnames  <- dir(dsno, pattern=paste('fort.',xtypes[1],sep=''))
  touDstr <- substr(fnames,11,20)                                            # POSIX double [as CHARACTER] representation
  tou <- as.POSIXct(as.numeric(touDstr), origin=originTimeStampTZ, tz='GMT') # POSIXct: all available output time
  tid <- which(tou > RF$timenow & tou <= RF$timenex)                         # INTEGER
  
  
  # init Els structure
  Els <- vector('list',length(xtypes))  
  for (i in 1:length(Els)) {
    Els[[i]] <- list()
    Els[[i]]$time <- tou[tid]
    Els[[i]]$stack <- list()
    for (it in 1:length(tid)) {
      Els[[i]]$stack[[it]] <- matrix(NA, n, m)
    }
  }
  names(Els) <- xtypes

  # get forecast data
  imStrs <- formatC(1:m,width=5,flag="0",format="d")
  for (iv in 1:length(xtypes)) {
    fname <- 'fort._____xxxxxxxxxx'
    substring(fname,6,5+nchar(xtypes[iv])) <- xtypes[iv]
    for (it in 1:length(tid)) {
      substring(fname,11,20) <- touDstr[tid[it]]
      for (im in 1:m) {
        dsno  <- file.path(RF$path$oufld,'maps',imStrs[im])
        Els[[iv]]$stack[[it]][,im] <- readBin(file.path(dsno, fname),
                                              what='double', n=n)    
      } # end im
    } # end it
  } # end iv
  return(Els)
}
