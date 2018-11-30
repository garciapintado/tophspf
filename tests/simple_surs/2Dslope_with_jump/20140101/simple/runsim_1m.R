# PROGRAM runsim_1m
# simple launcher to prepare m=1 simulation - with Geoclaw as test case
# m=1,no assimilation, 1 single inflow time series upstrem
# 2D domain with simple constant slope including a jump + valley + incised river

library(hydrosim)
library(tseries)
library(colorRamps)
options(digits=6)
prj <- '+init=epsg:27700'

HOME  <- Sys.getenv('HOME')
CLAW  <- Sys.getenv('CLAW')
CLAW  <- '/Users/pt902904/Documents/hydrology/clawpack-5.2.0'
PERM  <- Sys.getenv('PERM')

RCLAW <- file.path(CLAW,'Rclaw')
Rsrcs <- dir(RCLAW, pattern='.R$')
for (i in 1:length(Rsrcs)) {
  source(file.path(RCLAW,Rsrcs[i]))
}

#dsndat <- file.path(HOME,'Documents/hydrology/clawpack-5.2.0/geoclaw/examples/flood')
dsndat <- file.path(HOME,'Documents/hydrology/tophspf/tests/simple_surs')
region <- '2Dslope_with_jump'
event  <- '20140101'
scn    <- 'simple'

dsnsim <- file.path(dsndat,region,event,scn)

cwd <- getwd()
setwd(dsnsim)

# prepare foo domain
ncols <-   70                                      # longitude
nrows <-   200                                     # latitude
slope <- 0.002
dy    <- dx <- 50                                 # m
ll    <- c(0,0)                                   # lower-left coordinate of the domain (i.e. corner of the ll pixel)
zmax  <- 100                                      # m (height of the northmost pixel in the domain)
prjCRS   <- '+init=epsg:27700'                    # foo projection

makeValley <- function(xseq,depth) {
  w     <- max(xseq) - min(xseq)
  xbar <- min(xseq) + w/2
  abs(xseq - xbar) * depth / (w/2)
}


xseq    <- ll[1] - 0.5*dx + c(1:ncols)*dx         # create DEM (matrix stored in row order from N to S)
yseq    <- ll[2] - 0.5*dy + c(1:nrows)*dy
ryseq   <- yseq[nrows:1]
zseqY   <- zmax + (ryseq-ryseq[1])*slope          # [200 ... 181], dz=-1
rowbar  <- trunc(nrows/2)
zseqY[rowbar:nrows] <- zseqY[rowbar:nrows] - 5       # create jump
zseqX   <- makeValley(xseq, 5)                    # valley-shape domain
zseqX   <- rep(0,ncols)                           # simple downward-slope domain
zY <- matrix(zseqY,nrow=nrows,ncol=ncols)
zX <- matrix(zseqX,nrow=nrows,ncol=ncols,byrow=TRUE)
z <- zX + zY
bb      <- matrix(c(xseq[1]-dx/2,xseq[ncols]+dx/2,
                    yseq[1]-dy/2,yseq[nrows]+dy/2), nrow=2, byrow=TRUE, dimnames=list(c('x','y'),c('min','max')))
DEMsp   <- createSGDF(bb, res=dx, z=as.numeric(t(z)), proj4=prjCRS)
if (1 > 2) {
    library(colorRamps)
    plotGmeta(layer=DEMsp, col=matlab.like(50))
    digitGmeta(layer=DEMsp, type='line')
}
G <- getgmeta6.SP(DEMsp)

rm(ncols,nrows,dy,ll,xseq,yseq,ryseq,rowbar,z,bb)

source('readFF.R') # read Flood Forecast general parameters

                                        # export DEM to ArcGIS format
writeGDAL(DEMsp, fname=FF$demf, drivername='AAIGrid',
          options=c('NODATA=-9999.00','DECIMAL_PRECISION=4'))   # ArcGIS

# create bci file
nbci <- 2
bci <- matrix(NA,nbci,4)
bci[1,] <- c(G$xseq[15],G$xseq[15],rep(G$ryseq[20],2))
bci[2,] <- c(rep(G$xseq[5],2),rep(G$ryseq[20],2))
bcinames <- c('alongN','centerP')
bci <- data.frame(bci,to=c('s','o'),bcinames,
                  stringsAsFactors=FALSE)                                                # inflow are considered to have a direction towards 'e','s','w','n', respectively coded as 1:4 in geoclaw
write.table(bci, file='bci.data', quote=FALSE, row.names=FALSE, col.names=FALSE)
FF$bci <- data.frame(x=rowMeans(bci[,1:2]),
                     y=rowMeans(bci[,3:4]),bci)
coordinates(FF$bci) <- c('x','y')

# Create bdy file
tboundsPOSIXlt <- strptime(c(FF$starttimeStr, FF$endtimeStr),
                           format='%Y-%m-%d %H:%M:%S', tz='UTC')
tbounds <- as.double(tboundsPOSIXlt)
FF$starttime <- tbounds[1]                                                       # double :: POSIX start of simulation time [s] (from 1970-01-01 00:00:00 UTC)
FF$endtime   <- tbounds[2]                                                       # double :: POSIX end of simulation time [s]   (from 1970-01-01 00:00:00 UTC)
FF$inctime   <- FF$endtime - FF$starttime                                        # double :: total simulation time [s]

tseqPOSIX <- seq(tboundsPOSIXlt[1],tboundsPOSIXlt[2],by=60)                      # 15' foo dataset
nt <- length(tseqPOSIX)                                                          # 235
#v1 <- rep(0,nt)                                                                  # m3/s
#v1[100:nt] <- 200*exp(-(c(100:nt)-150)^2/20^2)
v1 <- 200*exp(-as.numeric(tseqPOSIX-tseqPOSIX[round(0.65*nt)])^2/18000^2)
v2 <- 0                                                                          # no inflow at this point
bdyIRTS <- irts(time=tseqPOSIX,value=cbind(v1,v2))

# create gauges (SpatialPointsDataFrame class)
FF$gauges <- read.table(file.path(dsnsim,FF$gaugesfR), header=FALSE, col.names=c('x','y','name'), as.is=TRUE)
coordinates(FF$gauges) <- c('x','y')
proj4string(FF$gauges) <- prjCRS
FF$gauges@data['z'] <- over(FF$gauges,DEMsp)

# create simple initial conditions:
xy  <- coordinates(DEMsp)
#xyc <- rowMeans(bbox(DEMsp))
xyc <- c(100,300)
hstart <- 3 * exp(- hypot(xy[,1] - xyc[1],xy[,2] - xyc[2])/10^2) +
          5 * exp(- hypot(xy[,1] - 2500,xy[,2] - 6000)/10^2)
hstart[hstart < 1.0E-03] <- 0.0
FF$qinitsp <- DEMsp
FF$hstart  <- DEMsp
FF$hstart@data[,1] <- hstart
FF$qinitsp@data[,1] <- hstart
writeGDAL(FF$qinitsp, fname=file.path(dsnsim,'hini.asc'), drivername='AAIGrid')   # ArcGIS

# plotGmeta(layer=qinitsp, col=matlab.like(50))
# digitGmeta(layer=qinitsp, type='points')

# prepare Geoclaw parameters
source(file.path(dsnsim,'readGC.R'))                                      # read Geoclaw parameters

# write input to Geoclaw
FF$timenow <- FF$starttime

#GC$rundata@claw@order <- as.integer(2)                                   #scn3,5
#GC$rundata@claw@transverse_waves <- as.integer(0)                        #scn5,6
#GC$rundata@claw@use_fwaves <- FALSE                                      #scn6
#GC$rundata@claw@use_fwaves <- TRUE                                       #scn3,5,7
#GC$rundata@claw@transverse_waves <- as.integer(2)                        #scn7

writeGC(GC$rundata, dsn=dsnsim)

#grid2Topotype(file.path(dsnsim,FF$demf),       GC$demf@fname,   topotype=3)
#grid2Topotype(file.path(dsnsim,'habsini.asc'), GC$qinitf@fname, topotype=1)
writeGDAL(FF$qinitsp, fname=GC$qinitf@fname, drivername='XYZ')   # x,y,z values, one per line in standard order from NW corner to SE


GC$bdy <- cbind(as.double(tseqPOSIX-tseqPOSIX[1]),bdyIRTS[[2]])
colnames(GC$bdy) <- c('time',bcinames)
write.table(GC$bdy, file='bdy.data', quote=FALSE, col.names=TRUE, row.names=FALSE)                                        # export for Geoclaw

# integrate Geoclaw
timenow <- Sys.time()
runclaw(xclawcmd=FF$modexe, GC, rundir=dsnsim, outdir=FF$path$dirroot)
runtime <- Sys.time() - timenow


# Get & plot model output
# -----------------------
# plot forecast at gauges
FF$forecast <- readGaugeSim(file.path(dsnsim,'gcresults','fort.gauge'),
                            FF$gauges$name,
                            FF$starttime)
for (ig in 1:length(FF$forecast)) {
 if (!is.null(FF$forecast[[ig]])) {
     X11(); #plot(FF$forecast[[ig]][,c('h','eta')], main=FF$gauges$name[ig])
             plot(FF$forecast[[ig]], main=FF$gauges$name[ig])
 }
}

# observe specific level versus inflow at a gauge
ig <- 2
X11();plot(FF$forecast[[ig]][,'h'], main=FF$gauges$name[ig], ylim=c(0,0.5))
lines(bdyIRTS[[1]],bdyIRTS[[2]][,1]/500, col='blue')      # re-scaled inflow

# Plot sequence of grids showing water depth

plotGmeta(layer=DEMsp, col=terrain.colors(100))
for (i in 1:232) {
    cat('plotting output slide:',i,'\n')
    SGDF <- readGeoclawfor(file.path(dsnsim,'gcresults'), fname=paste('fort.q',formatC(i,format='d',width=4,flag=0),sep=''),
                           G, proj4=prjCRS)
    #plotGmeta(layer=DEMsp, col=terrain.colors(100), add=TRUE)
    plotGmeta(layer=SGDF$h,col=matlab.like(100), zlim=c(0,0.5),add=TRUE)
    points(FF$gauges, pch=21, col='red')
    points(FF$bci,    pch=21, col='green')

    # Sys.sleep(0.05)
}

inith <- DEMsp
inith@data[,1] <- SGDF$habs@data[,1] - DEMsp@data[,1]
plotGmeta(layer=DEMsp)
plotGmeta(layer=inith, col=matlab.like(100), add=TRUE)

bth <- DEMsp
bth@data[,1] <- (SGDF$habs@data[,1] - SGDF$h@data[,1])     # model bathymetry
btherr  <- DEMsp
btherr@data[,1] <- DEMsp@data[,1] - bth@data[,1]           # error in model bathymetry with respect to input bathymetry
plotGmeta(layer=btherr,  col=matlab.like(100))
bte <- digitGmeta(layer=btherr, type='points')


# post analysis with fixed grid output
# read surs

sursfiles <- dir(file.path(dsnsim,'gcresults'),pattern='fort.surs')
sursSGDF <- DEMsp
plotGmeta(layer=DEMsp, col=matlab.like(50))
#for (i in 1:length(sursfiles)) {
for (i in 1:1) {
 cc <-    scan(file.path(dsnsim,'gcresults',sursfiles[i]))
 sursSGDF@data[,1] <- scan(file.path(dsnsim,'gcresults',sursfiles[i]))
 plotGmeta(layer=sursSGDF, col=matlab.like(50), add=TRUE)
 Sys.sleep(1)
}

                                        # read fixgrid


fgfiles <- dir(file.path(dsnsim,'gcresults'),pattern='fort.fg')

fgrid <- readGeoclawfixgrid(dsno =file.path(dsnsim,'gcresults'),
                            fname=fgfiles[1],
                            proj4=prjCRS)

plotGmeta(layer=fgrid$SGDF['h'], col=matlab.like(100))        # water depth
plotGmeta(layer=fgrid$SGDF['hu'], col=matlab.like(100))       # momentum N-S
plotGmeta(layer=fgrid$SGDF['eta'], col=matlab.like(100))      # absolute water height
plotGmeta(layer=fgrid$SGDF['bathy'], col=matlab.like(100))    # absolute water height

fgerr <- fgrid$SGDF['bathy']
fgerr@data[,1] <- fgerr@data[,1] - DEMsp@data[,1]
plotGmeta(layer=fgerr, col=matlab.like(100))

# this is better as actually check the computational grid against the source DTM:
fgerr@data[,1] <- fgrid$SGDF['eta']@data[,1] - fgrid$SGDF['h']@data[,1] - DEMsp@data[,1]
plotGmeta(layer=fgerr, col=matlab.like(100))

fgfiles <- dir(file.path(dsnsim,'gcresults'),pattern='fort.fg')
plotGmeta(layer=DEMsp, col=terrain.colors(100))
wst <- NA
wsv <- NA

for (i in 1:(length(fgfiles))) {
#    for (i in 1:5) {
  fgrid <- readGeoclawfixgrid(dsno =file.path(dsnsim,'gcresults'),
                            fname=fgfiles[i],
                              proj4=prjCRS, mask=TRUE)
  fgrid$time <- as.POSIXct(FF$timenow + fgrid$time, origin='1970-01-01 00:00:00', tz='GMT')

  #plotGmeta(layer=DEMsp, col=terrain.colors(100), add=TRUE)
                                        # plotGmeta(layer=fgrid$SGDF['eta'],col=matlab.like(100), zlim=range(DEMsp@data[,1]),add=TRUE)
  cat('plotting fixed_grid:',i,'\n')
  plotGmeta(layer=fgrid$SGDF['h'],col=matlab.like(100), add=TRUE)
  points(FF$gauges, pch=21, col='red')
  points(FF$bci,    pch=21, col='green')
  #Sys.sleep(1)

  wst[i] <- fgrid$time
  wsv[i] <- sum(fgrid$SGDF['h']@data[,1], na.rm=TRUE) * dx^2
  #dev.off()
}
ws <- irts(time=wst,value=wsv)

# water balance:
# input
nt <- length(bdyIRTS$time)
dt <- as.numeric(difftime(bdyIRTS$time[-1],bdyIRTS$time[-nt], units='secs'))
wi <- rowSums(bdyIRTS[[2]])
wi <- c(0,cumsum((wi[-1] + wi[-nt]) / 2 * dt))
#wi <- wi + sum(hstart)*dx^2 # add initial condition water level integrated over the surface
wi <- irts(time=bdyIRTS$time, value=wi+wsv[1])
X11();plot(wi, ylim=c(0,max(wi$value)))                                                                                 # cumulative [m3 along time] input into the model
lines(ws, col='blue')
abline(h=0, col='grey', lty=3)
                                        # model stored




fgh <- fgrid$SGDF['']
fgh@data[,1]   <- fgh@data[,1] - DEMsp@data[,1]
plotGmeta(layer=fgh,   col=matlab.like(100))

                                        # check diff with source DT



