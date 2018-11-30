ts2DA <- function(xm3s, field='rnof', glike=NULL, stat='E', mo=NULL) {
 # map timeseries to a format directly readable by the floodDA() as water inflows
 # based on glike select the m best simulations for each boundary conditions

 # note: the glike-based selection mechanism is not to be used as real forecast mode,
 #       where we don't know which are the best simulations, but all (m=NULL) are selected
 #       it is just intended to be used for testing an uncalibrated model in one event as linked with the lisflood DA software

  mi <- ncol(xm3s[[field]][[1]]$value)             # input size ensemble
  nu  <- length(glike)
  olst <- vector('list', nu)                       # output list with irregularly sampled time series (irts)

  if (length(xm3s[[field]]) != nu)
    stop('ts2DA:: --ERR 001--')
  
  if (!is.null(glike))
    names(olst) <- names(glike)

  for (iu in 1:nu) {
    olst[[iu]]$irts <- xm3s[[field]][[iu]]              # multivariate irts 
  }

  if (is.null(mo))
    mo <- mi

  if (mo > mi)
    stop('ts2DA:: ERR002 : mo > mi')

  if (!is.null(glike)) {  # select the best m simulations in the ensemble
    for (i in 1:nu) {
      stats <- glike[[i]][stat,]
      imd <- order(stats, decreasing=TRUE)[1:mo]      # warning: assumes higher 'stats' is better
      olst[[i]]$irts[[2]] <- olst[[i]]$irts[[2]][,imd]
      olst[[i]]$imrf      <- imd
      olst[[i]][[stat]]   <- glike[[i]][stat,imd]
    }
  }
  return(olst)
}
