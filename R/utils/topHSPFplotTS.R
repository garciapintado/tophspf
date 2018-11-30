topHSPFplotTS <- function(x, xm3s, dsnplt, plothrus=FALSE, qobs=NULL, hrunames, field='value', im=NULL) {
 # +++ purpose +++
 # plot mean diagnosis timeseries for each store in the topHSPF model
 # a single simulation serves as input
 #
 # x    :: list with time series from topHSPF forecast ensemble [mm]
 # xm3s :: list with forecast elements $suro $ifwo $agwo $rnof catchment-integrated [m3s] 
 # simseqT     :: POSIXct object to label the x (time) axis
 # plothrus :: flag. TRUE to plot timeseries at each subwatershed,
 #                   FALSE to plot mean watershed timeseries
 # qobs     :: if not NULL, observed hydrographs will be added to the runoff plots
  
 if (!file.exists(dsnplt))
  dir.create(dsnplt, recursive=TRUE)

 nhrus <- length(hrunames)
 if (length(x$rain) - 1 != nhrus)
   stop('forecast should have 1 element more than nhrus')

 times <- x$rain[[1]]$time
 dtih   <- difftime(times[2],times[1], units='h')

 for (iu in 1:nhrus) {   # one PDF / HRU
   if (!plothrus)
     iu <- nhrus + 1
   #iu <- ifelse(nhrus == 0, 1, nhrus+1)  # colindex for mean values in watershed
   cat('plotting iu:',iu,'\n')
   nrows <- 3  # up to 12 plots
   ncols <- 4
   paper.width  <- 7.5*ncols      # cm
   paper.height <-   6*nrows      # cm
   fname <- paste('rhu_',iu,'_',hrunames[iu],'.pdf',sep='')

   pdf(file=file.path(dsnplt,fname), width=paper.width/2.54+0.1, height=paper.height/2.54+0.1)
   mat <- matrix(1:(nrows*ncols), nrow=nrows, ncol=ncols, byrow=TRUE)
   layout(mat=mat,widths=lcm(rep(paper.width/ncols,ncols)), heights=lcm(rep(paper.height/nrows,nrows)), respect=TRUE)
   par.default <- par(no.readonly = TRUE)
   par(cex=0.6)
   par(mgp=c(2.3,1,0))  # para regular la distancia (2.5) de los titulos de ejes a la figura
   par(mai=c(1.1,1.5,0.6,0.1)/2.54)

   # plot 1: rainfall & potential evapotranspiration
   ylim <- c(0,range(c(x$rain[[iu]][[field]],x$petinp[[iu]][[field]]))[2])
   plot(x$rain[[iu]]$time, x$rain[[iu]][[field]][,im], type='l', col='navyblue', main="Rainfall and PEVTP",
        xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
        ylab="Watershed Mean [mm]", ylim=ylim, font.lab=2)
   lines(x$petinp[[iu]]$time, x$petinp[[iu]][[field]][,im], col='green4')
   legend(x="topleft",legend=c("rain","petinp"),col=c("navyblue","green4"),
          lty=1, bty="n")

   # water balances
   # ---
   # plot 2: balance around interception storage (ceps) [mm]
   #cat('plot 2\n')
   supycum <- cumsum(x$supy[[iu]][[field]][,im])                # input to CEPS
   ceps0   <- x$ceps[[iu]][[field]][1,im] - x$ceps[[iu]][[field]][1,im] +
              x$suri[[iu]][[field]][1,im] + x$cepe[[iu]][[field]][1,im]
   ceps    <- x$ceps[[iu]][[field]][,im]                        # CEPS storage  (ceps_ini=0.0 mm. TODO: Input from PWAT-STATE1)
   suricum <- cumsum(x$suri[[iu]][[field]][,im])                       # output toward SURS
   cepecum <- cumsum(x$cepe[[iu]][[field]][,im])                       # output to evapotranspiration

   plot(times,supycum,type='l',main="CEPS balance",
        xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
        ylab="Watershed Mean [mm]",lwd=3, font.lab=2)
   lines(times,ceps,col='blue')
   lines(times,suricum,col='orange')
   lines(times,cepecum,col='cyan')
   lines(times,ceps-ceps0+suricum+cepecum,col='green',lty=3)  #should match supycum
   points(times[1],ceps0,col='blue')
   legend(x="topleft",legend=c("supycum","ceps","suricum","cepecum","balance"),
          col=c("black","blue","orange","cyan","green"),
          lty=c(1,1,1,1,3), lwd=c(3,1,1,1,1),bty="n")

   # plot 3: balance around surface detention storage (surs)
   surs0    <- x$surs[[iu]][[field]][1,im] - x$suri[[iu]][[field]][1,im] +
               x$suro[[iu]][[field]][1,im] + x$ifwi[[iu]][[field]][1,im] +
               x$uzi[[iu]][[field]][1,im]  + x$infil[[iu]][[field]][1,im]
   surs     <- x$surs[[iu]][[field]][,im]                        #SURS storage (surs_ini=0.0 mm. TODO: Input from PWAT-STATE1)
   surocum  <- cumsum(x$suro[[iu]][[field]][,im])                #surface lateral outflow
   ifwicum  <- cumsum(x$ifwi[[iu]][[field]][,im])                #output to IFWS
   uzicum   <- cumsum(x$uzi[[iu]][[field]][,im])                 #output to UZS
   infilcum <- cumsum(x$infil[[iu]][[field]][,im])               #output to other storages

   plot(times,suricum,type='l',main="SURS balance",
              xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
              ylab="Watershed Mean [mm]",lwd=3, font.lab=2)
   lines(times,surs,col='blue')
   lines(times,surocum,col='navyblue')
   lines(times,ifwicum,col='red')
   lines(times,uzicum,col='orange')
   lines(times,infilcum,col='violet')
   lines(times,surs-surs0+surocum+ifwicum+uzicum+infilcum,col='green',lty=3)  #should match suricum input to SURS
   points(times[1],surs0,col='blue')
   legend(x="topleft",legend=c("suricum","surs","surocum","ifwicum","uzicum","infilcum","balance"),
          col=c("black","blue","navyblue","red","orange","violet","green"),
          lty=c(1,1,1,1,1,1,3), lwd=c(3,1,1,1,1,1,1),bty="n")

   # plot 4: balance around interflow storage (ifws)
   #cat('plot 4\n')
   ifws0    <- x$ifws[[iu]][[field]][1,im] - x$ifwi[[iu]][[field]][1,im] + x$ifwo[[iu]][[field]][1,im]
   ifws     <- x$ifws[[iu]][[field]][,im]                        #IFWS storage (ifws_ini=0.0 mm. TODO: Input from PWAT-STATE1)
   ifwocum  <- cumsum(x$ifwo[[iu]][[field]][,im])                #interflow lateral outflow

   plot(times,ifwicum,type='l',main="IFWS balance",
              xlab=paste("Inner timestep (",dtih," [h])",sep=""),
              ylab="Watershed Mean [mm]",lwd=3, font.lab=2)
   lines(times,ifws,col='blue')
   lines(times,ifwocum,col='navyblue')
   lines(times,ifws-ifws0+ifwocum,col='green',lty=3)                    #should match ifwicum input to IFWS
   points(times[1],ifws0,col='blue')
   legend(x="topleft",legend=c("ifwicum","ifws","ifwocum","balance"),
          col=c("black","blue","navyblue","green"),
          lty=c(1,1,1,3), lwd=c(3,1,1,1), bty="n")

   # plot 5: balance around upper zone storage (uzs)
   uzs0 <- x$uzs[[iu]][[field]][1,im] - x$uzi[[iu]][[field]][1,im] +
           x$uzet[[iu]][[field]][1,im] + x$perc[[iu]][[field]][1,im]    # this estimate has to be substituted by the real init conditions
   uzs      <- x$uzs[[iu]][[field]][,im]                         # UZS storage
   uzetcum  <- cumsum(x$uzet[[iu]][[field]][,im])                # output to evapotranspiration
   perccum  <- cumsum(x$perc[[iu]][[field]][,im])                # output to other storages (adds to infilcum)

   plot(times,uzicum,type='l',main="UZS balance",          # input to UZS storage
        xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
        ylab="Watershed Mean [mm]",lwd=3, font.lab=2)
   lines(times,uzs,col='blue')
   lines(times,uzetcum,col='cyan')
   lines(times,perccum,col='orange')
   lines(times,uzs-uzs0+uzetcum+perccum,col='green',lty=3)             #should match uzicum input to UZS
   points(times[1],uzs0,col='blue')
   legend(x="topleft",legend=c("uzicum","uzs","uzetcum","perccum","balance"),
          col=c("black","blue","cyan","orange","green"),
          lty=c(1,1,1,1,3), lwd=c(3,1,1,1,1), bty="n")

  # plot 6: balance around lower zone storage (lzs)
  lzicum   <- cumsum(x$lzi[[iu]][[field]][,im])                 #input to LZS storage
  lzs0     <- x$lzs[[iu]][[field]][1,im] - x$lzi[[iu]][[field]][1,im] + x$lzet[[iu]][[field]][1,im]
  lzs      <- x$lzs[[iu]][[field]][,im]                         #LZS storage
  lzetcum  <- cumsum(x$lzet[[iu]][[field]][,im])                #output to evapotranspiration

  plot(times,lzicum,type='l',main="LZS balance",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Watershed Mean [mm]", font.lab=2,
       lwd=3, ylim=c(0,max(lzicum,lzs0,lzetcum)))
  lines(times,lzs,col='blue')
  lines(times,lzetcum,col='cyan')
  lines(times,lzs-lzs0+lzetcum,col='green',lty=3)                     #should match lzicum input to LZS
  points(times[1],lzs0,col='blue')
  legend(x="topleft",legend=c("lzicum","lzs","lzetcum","balance"),
         col=c("black","blue","cyan","green"),
         lty=c(1,1,1,3), lwd=c(3,1,1,1), bty="n")

  # plot 7:  balance around active groundwater storage (agws)
  agwicum  <- cumsum(x$agwi[[iu]][[field]][,im])                #input to agws
  agws0    <- x$agws[[iu]][[field]][1,im] - x$agwi[[iu]][[field]][1,im] +
              x$agwo[[iu]][[field]][1,im] + x$agwet[[iu]][[field]][1,im]
  agws     <- x$agws[[iu]][[field]][,im]                        #AGWS storage
  agwocum  <- cumsum(x$agwo[[iu]][[field]][,im])                #groundwater lateral outflow
  agwetcum <- cumsum(x$agwet[[iu]][[field]][,im])               #output to evapotranspiration

  plot(times,agwicum,type='l',main="AGWS balance",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Watershed Mean [mm]", font.lab=2,
       lwd=3, ylim=c(0,max(agwicum,agws0,agwocum,agwetcum)))
  lines(times,agws,col='blue')
  lines(times,agwocum,col='navyblue')
  lines(times,agwetcum,col='cyan')
  lines(times,agws-agws0+agwocum+agwetcum,col='green',lty=3)          #should match agwicum input to AGWS
  points(times[1],agws0,col='blue')
  legend(x="topleft",legend=c("agwicum","agws","agwocum","agwetcum","balance"),
         col=c("black","blue","navyblue","cyan","green"),
         lty=c(1,1,1,1,3), lwd=c(3,1,1,1,1), bty="n")

  # plot 8: evapotranspiration balance
  #cat('plot 8\n')
  taetcum <- cumsum(x$taet[[iu]][[field]][,im])                    #total evapotranspitation
  rempetcum <- cumsum(x$rempet[[iu]][[field]][,im])                #remaining non-satisfied evapotranspitation
  petinpcum <- cumsum(x$petinp[[iu]][[field]][,im])                #input potential evapotranspitation

  plot(times,taetcum,type='l',main="Evapotranspiration balance",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Watershed Mean [mm]",
       lwd=3,ylim=c(0,max(petinpcum)), font.lab=2)
  lines(times,cepecum,col='blue')
  lines(times,uzetcum,col='violet')
  lines(times,lzetcum,col='orange')
  lines(times,agwetcum,col='cyan')
  lines(times,cepecum+uzetcum+lzetcum+agwetcum,col='green',lty=3)          #should match taetcum
  lines(times,petinpcum,col='navyblue',lwd=3,lty=2)
  lines(times,taetcum+rempetcum,col='yellow',lty=2)                        #should match petinpcum
  legend(x="topleft",legend=c("taetcum","cepecum","uzetcum","lzetcum","agwetcum","balance1","petinpcum","balance2"),
         col=c("black","blue","violet","orange","cyan","green","navyblue","yellow"),
         lty=c(1,1,1,1,1,3,2,2), lwd=c(3,1,1,1,1,1,3,1), bty="n")

  # plot 9: lzs/agws and lost water balance
  igwicum <- cumsum(x$igwi[[iu]][[field]][,im])

  plot(times,infilcum+perccum,type='l', main="lzs/agws distribution and lost water balance",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""), font.lab=2,
       ylab="Watershed Mean [mm]", lwd=3, ylim=c(0,max(infilcum+perccum)))
  lines(times,infilcum,col='orange')
  lines(times,perccum,col='red')
  lines(times,igwicum,col='blue')
  lines(times,lzicum,col='blue4')
  lines(times,agwicum,col='cyan')
  lines(times,igwicum+lzicum+agwicum,col='green',lty=3)          #should match infilcum+perccum
  legend(x="topleft",legend=c("infilcum+perccum","infilcum","perccum","igwicum","lzicum","agwicum","balance"),
         col=c("black","orange","red","blue","blue4","cyan","green"),
         lty=c(1,1,1,1,1,1,3), lwd=c(3,1,1,1,1,1,1), bty="n")

  # plot 10: global balance and all stores
  #cat('plot 10\n')
  raincum <- cumsum(x$rain[[iu]][[field]][,im])
  plot(times,raincum,type='l',main="Global Balance and Stores",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Watershed Mean [mm]", lwd=3, ylim=c(0,max(raincum)), font.lab=2)
  lines(times,taetcum,col='cyan')
  lines(times,surocum,col='blue3')
  lines(times,ifwocum,col='blue4')
  lines(times,agwocum,col='blue2')
  lines(times,ceps,col='pink')
  lines(times,surs,col='yellow')
  lines(times,ifws,col='brown')
  lines(times,uzs,col='orangered')
  lines(times,lzs,col='lightgoldenrod')
  lines(times,agws,col='orange')
  lines(times,taetcum+surocum+ifwocum+agwocum+
              ceps +surs +ifws +uzs +lzs +agws -
             (ceps0+surs0+ifws0+uzs0+lzs0+agws0),col='green',lty=3)          #should match raincum
  legend(x="topleft",legend=c("raincum","taetcum","surocum","ifwocum","agwocum",
                             "ceps","surs","ifws","uzs","lzs","agws","balance"),
        col=c("black","cyan","blue3","blue4","blue2",
              "pink","yellow","brown","orangered","lightgoldenrod","orange","green"),
        lty=c(rep(1,11),3), lwd=c(3,rep(1,11)), bty="n")
# browser()
  # plot 11: total runoff & components :: note runoff list is assumed to be a formal irts
  y <- xm3s$rnof[[iu]]$value[,im]
  ymax <- max(y)
  if (!is.null(qobs)) {
    if (!is.null(qobs[[iu]])) {ymax <- max(ymax,qobs[[iu]]$irts$value)} }
  plot(times,xm3s$rnof[[iu]]$value[,im],type='l',main="Total Runoff & Components",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Runoff [m3/s]", lwd=3, ylim=c(0,ymax),col='navyblue', font.lab=2)
  lines(times, xm3s$suro[[iu]]$value[,im],col='blue')                       #surface lateral outflow
  lines(times, xm3s$ifwo[[iu]]$value[,im],col='green')                      #interflow lateral outflow
  lines(times, xm3s$agwo[[iu]]$value[,im],col='cyan')                       #groundwater lateral outflow
  if (!is.null(qobs)) {
    if (!is.null(qobs[[iu]])) {
      lines(qobs[[iu]]$irts$time,qobs[[iu]]$irts$value, col='grey')
    }
  }
  legend(x="topleft",legend=c("total","suro","ifwo","agwo"),
         col=c("navyblue","blue","green","cyan"),
         lty=1, lwd=c(3,1,1,1), bty="n")

  # plot 12: total runoff and components
  #cat('plot 12\n')
  y <- cumsum(xm3s$rnof[[iu]]$value[,im])*dtih*3600
  plot(times,y,type='l',main="Total Runoff & Components",
       xlab=paste("Inner Timestep (",dtih," [h])",sep=""),
       ylab="Runoff [m3]", lwd=3, ylim=c(0,max(y)),col='navyblue', font.lab=2)
  lines(times,cumsum(xm3s$suro[[iu]]$value[,im]),col='blue')
  lines(times,cumsum(xm3s$ifwo[[iu]]$value[,im]),col='green')
  lines(times,cumsum(xm3s$agwo[[iu]]$value[,im]),col='cyan')
  legend(x="topleft",legend=c("total","suro","ifwo","agwo"),
         col=c("navyblue","blue","green","cyan"),
         lty=1, lwd=c(3,1,1,1), bty="n")
  #cat('all plotted\n')
  dev.off()
  if (!plothrus)
    break
 } # end for iu
}
