# temporary file: data preprocess for Severn-Avon 1D kinematic simulation & topHSPF

library(hydrosim)

dsnm <- file.path(dsn,region,'watmaps')


nexb <- readGDAL(file.path(dsnm,'nexburn1000_dtme_crp.tif'))       # nx=191,ny=124
demf <- readGDAL(file.path(dsnm,'demfilled_patch.tif'))       # nx=191,ny=124
mask <- readGDAL(file.path(dsnm,'mask.tif'))       # nx=191,ny=124
ldd  <- readGDAL(file.path(dsnm,'lddfilled.tif'))       # nx=191,ny=124
aflx <- readGDAL(file.path(dsnm,'accufluxfilled.tif'))       # nx=191,ny=124
cmsk <- readGDAL(file.path(dsnm,'cmsk.tif'))       # nx=191,ny=124
cwid <- readGDAL(file.path(dsnm,'cwid.tif'))       # nx=191,ny=124
cdep <- readGDAL(file.path(dsnm,'cdep.tif'))       # nx=191,ny=124
outs <-

#difference betwwen demf and nexb:
difm <- nexb
difm@data[,1] <- demf@data[,1] - nexb@data[,1]
plotGmeta(layer=difm)
summary(difm@data[,1])


plotGmeta(layer=demf)
plotGmeta(layer=demf, zlim='strloc',xlim=c(350000,380000),ylim=c(260000,280000), col=matlab.like(100))


PSP <- tmpSP
PSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','P.for'), n=G$cells, what='double')
plotGmeta(layer=PSP, col=matlab.like(100))

ckSP <- tmpSP
ckSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','ck.for'), n=G$cells, what='double')
plotGmeta(layer=ckSP, col=matlab.like(100))
digitGmeta(layer=ckSP, type='points')

QbSP <- tmpSP
QbSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','Qb.for'), n=G$cells, what='double')
plotGmeta(layer=QbSP, col=matlab.like(100))
digitGmeta(layer=QbSP, type='points')


QaSP <- tmpSP
QaSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','Qa.for'), n=G$cells, what='double')
plotGmeta(layer=QaSP, col=matlab.like(100))
digitGmeta(layer=QaSP, type='points')

hpreSP <- tmpSP
hpreSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','hpre.for'), n=G$cells, what='double')
plotGmeta(layer=hpreSP, col=matlab.like(100))
digitGmeta(layer=hpreSP, type='points')

hposSP <- tmpSP
hposSP@data[,1] <- readBin(file.path(dsnsim,'input','00001','hpos.for'), n=G$cells, what='double')
plotGmeta(layer=hposSP, col=matlab.like(100))
digitGmeta(layer=hposSP, type='points')

dh <- tmpSP
dh@data[,1] = hposSP@data[,1] - hpreSP@data[,1]
plotGmeta(layer=dh, col=matlab.like(100))
digitGmeta(layer=dh, type='points')

tmpStr <- '0000000900'

agws  <- tmpSP
fname <- paste('fort.agws_',tmpStr,sep='')
agws@data[,1] <- readBin(file.path(dsnsim,'input','00001',fname), n=G$cells, what='double')
plotGmeta(layer=agws, col=matlab.like(100))
digitGmeta(layer=agws, type='points')

surs  <- tmpSP
fname <- paste('fort.surs',tmpStr,sep='')
surs@data[,1] <- readBin(file.path(dsnsim,'input','00001',fname), n=G$cells, what='double')
plotGmeta(layer=surs, col=matlab.like(100))
digitGmeta(layer=surs, type='points')

dsno <- file.path(dsnsim,'output','maps','00001')
dt <- 900
t0 <- tnow <- 900
nt <- 500

surs <-  demSP
fnames <- dir(dsno,pattern='fort.sursa')
for (it in 1:length(fnames)) {
 #tmpStr <- formatC(tnow,width=10,flag='0',format='d')
 #fname <- paste('fort.uzs__',tmpStr,sep='')
 fname <- fnames[it]
 surs@data[,1] <- readBin(file.path(dsno,fname), n=RF$G$cells, what='double')
 plotGmeta(layer=surs, col=matlab.like(100), add=TRUE)
 #tnow <- tnow + dt
 Sys.sleep(0.1)
}

dsno <- file.path(dsnsim,'output','maps','00001')
suro <-  demSP
fnames <- dir(dsno,pattern='fort.suro')
for (it in 1:length(fnames)) {
 #tmpStr <- formatC(tnow,width=10,flag='0',format='d')
 #fname <- paste('fort.uzs__',tmpStr,sep='')
 fname <- fnames[it]
 suro@data[,1] <- readBin(file.path(dsno,fname), n=RF$G$cells, what='double')
 plotGmeta(layer=suro, col=matlab.like(100), add=ifelse(it==1,FALSE,TRUE))
 #tnow <- tnow + dt
 Sys.sleep(0.1)
}

uzs <-  demSP
fnames <- dir(dsno,pattern='fort.uzs')
for (it in 1:nt) {
 #tmpStr <- formatC(tnow,width=10,flag='0',format='d')
 #fname <- paste('fort.uzs__',tmpStr,sep='')
 fname <- fnames[it]
 uzs@data[,1] <- readBin(file.path(dsnsim,'input','00001',fname), n=G$cells, what='double')
 plotGmeta(layer=uzs, col=matlab.like(100), add=TRUE)
 #tnow <- tnow + dt
 Sys.sleep(0.1)
}
