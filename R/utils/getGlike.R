getGlike <- function(nu, m, xls, yls, dti) {
   # get a stack of Generalized Likelihood statistics matrices. One [m x nstats] matrix / HRU
   # list, INTENT(IN) :: xls           length(nu) list with ensemble forecast irts for each iu
   # list, INTENT(IN) :: yls           length(nu) list with observation time series. $irts and $pos slot for each iu


  stats.nm <- c("runpid","ME","FLOWREL","MAE","RMSE","E","FV","lag.peak","PTLAG",
                 "dif.peak.abs","PDIFF","dif.peak","lambda","HMLE","trust.HMLE")
  if (length(xls) != nu)
    stop('getGlike:: forecast time series list length != nu')
  if (length(yls) != nu)
    stop('getGlike:: observation time series list length != nu')

  glike <- vector('list',nu)
  for (iu in 1:nu) {
    glike[[iu]] <- matrix(NA, nrow=length(stats.nm), ncol=m)
    rownames(glike[[iu]]) <- stats.nm
    for (im in 1:m) {
      glike[[iu]][,im] <- perfstatf(im, qobs = yls[[iu]]$irts, qsim = xls[[iu]][,im], dT = dti)
    }
  }
  names(glike) <- names(xls)
  return(glike)
} # end function getGlike()
