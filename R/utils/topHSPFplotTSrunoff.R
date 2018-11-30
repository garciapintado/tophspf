topHSPFplotTSrunoff <- function(qsim, qobs=NULL, glike=NULL, stat='E', mo=NULL, dsnplt=NULL) {
 # +++ purpose +++
 # plot the ensemble of runoff timeseries in the topHSPF model
 #
 #
 # qsim     :: list with runoff ensemble simulations from topHSPF run
 # qobs     :: list with observed runoff, if available
 # qlike    :: list with performance statistics, if available
 # stat     :: statistics name to select the best ensemble members
 # mtop     :: number of best members to identify as better
 # dsnplt   :: output folder

 if (!file.exists(dsnplt))
  dir.create(dsnplt)

 nhrus <- length(qsim)
 paper.width  <- 18      # cm
 paper.height <- 15      # cm

 if (!all.equal(names(qsim),names(qobs)))
   stop('topHSPFplotTSrunoff:: check qsim & qobs timeseries match')

 for (iw in 1:nhrus) {   # one PDF / HRU
   cat('plotting runoff ensemble for iw:',iw,'\n')
   m <- ncol(qsim[[iw]]$irts[[2]])
   #col <- rep('grey',m)
   #if (!is.null(glike))
   fname <- paste(names(qsim[iw]),'_runoff.pdf',sep='')
   pdf(file=file.path(dsnplt,fname), width=paper.width/2.54+0.1, height=paper.height/2.54+0.1)
   par.default <- par(no.readonly = TRUE)
   par(cex=1.0)
   par(mgp=c(2.6,1,0))  # para regular la distancia (2.5) de los titulos de ejes a la figura
   par(mai=c(1.5,2.3,0.8,0.1)/2.54)

   main <- gsub('_',' ',names(qsim[iw]))
   main <- substr(main,1,nchar(main)-2)
   ylim <- range(qsim[[iw]]$irts[[2]], na.rm=TRUE)
   if (!is.null(qobs))
     ylim <- range(ylim,qobs[[iw]]$irts[[2]], na.rm=TRUE)
   tlim <- range(qsim[[iw]]$irts[[1]])

   plot(tlim, c(1,1), xlim=tlim, ylim=ylim, main=main, type='n',
        las=1, xlab='',ylab='Runoff [m3/s]',font.lab=2, family='sans')

   for (im in 1:m) {
     lines(qsim[[iw]]$irts[[1]], qsim[[iw]]$irts[[2]][,im], col='grey')
   }
   lines(qsim[[iw]]$irts[[1]],rowMeans(qsim[[iw]]$irts[[2]]), col='darkred',lty=2)
   lines(qobs[[iw]]$irts[[1]], qobs[[iw]]$irts[[2]], col='navyblue')                 # 'm3.s'
   dev.off()
 }
 return(0)
}
