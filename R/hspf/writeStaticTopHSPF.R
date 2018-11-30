topHSPFwriteStatic <- function(RF, K01) {

        # RF:  mask, hrus, dem,
        # K01: ldd, aflx, cmsk, cwid, cdep
        # dsnscn [.tif] -> dsnsim [.bin]
  dsnsim <- file.path(RF$path$data,RF$region,RF$event,RF$scn)

  # LOGICAL mask  as integer (1 as .TRUE., 0 as .FALSE| always ASCII)
  maskSP <- readGDAL(file.path(dsnscn,RF$mask))
  write(as.integer(maskSP@data[,1]),
        file=file.path(dsnsim,'input','mask.asc'))

  # hrus map as integer (always ASCII) [0 == out of hydrological domain]
  # writeBin(as.integer(hrus), fhrus)   # cac <- readBin('f_hrus.map',what="integer", n=G$cols*G$rows, size=4)
  hruSP <- readGDAL(file.path(dsnscn,RF$hrus))
  write(as.integer(hruSP@data[,1]),
        file=file.path(dsnsim,'input','hrus.asc'))

  if (!is.null(K01)) {

    # dem map as binary double [NA out of domain]
    demSP <- readGDAL(file.path(dsnscn,K01$dem))
    writeBin(demSP@data[,1],
             file.path(dsnsim,'input','dem.bin')) # K = 8 Bytes * cols * rows / 1024  - Fortran (DP)
    # demV <- demSP@data[,1]
    # demV[is.na(demV)] <- max(demV,na.rm=TRUE) ! in case a SWE is used

    # ldd map integer [NA out of domain]
    lddSP <- readGDAL(file.path(dsnscn,K01$ldd)) # integer input : K =
    writeBin(lddSP@data[,1],
             file.path(dsnsim,'input','ldd.bin')) # K = 4 Bytes * cols * rows / 1024  - Fortran (I4B)

    # accuflux map integer [NA out of domain] [with 4 Bytes we can have until n <= 4.2E09 pixels]
    aflxSP <- readGDAL(file.path(dsnscn,K01$aflx)) # integer input : K =
    mode(aflxSP@data[,1]) <- 'integer'
    writeBin(aflxSP@data[,1],
             file.path(dsnsim,'input','aflx.bin')) # K = 4 Bytes * cols * rows / 1024  - Fortran (I4B)


    # channel mask for friction & subgrid parameters [0 out of domain as weel as floodplain]
    cmskSP <- readGDAL(file.path(dsnscn,K01$cmsk)) # double input
    mode(cmskSP@data[,1]) <- 'integer'
    writeBin(cmskSP@data[,1],
             file.path(dsnsim,'input','cmsk.bin')) # K = 4 Bytes * cols * rows / 1024  - Fortran (I4B)

    # channel width double [0 out of domain as well as floodplain]
    cwidSP <- readGDAL(file.path(dsnscn,K01$cwid))
    writeBin(cwidSP@data[,1],
             file.path(dsnsim,'input','cwid.bin')) # K = 8 Bytes * cols * rows / 1024  - Fortran (DP)

    # channel depth double [0 out of domain as well as floodplain]
    cdepSP <- readGDAL(file.path(dsnscn,K01$cdep))
    writeBin(cdepSP@data[,1],
             file.path(dsnsim,'input','cdep.bin')) # K = 8 Bytes * cols * rows / 1024  - Fortran (DP)

    #cdepSP <- cropSGDF(cdepSP,bbox(maskSP))
    #writeGDAL(cdepSP, fname=file.path(dsn,region,'watmaps',K01$cdep))
    # G geographical metadata file
    writeGmeta(RF$G, file.path(dsnsim,'input','G.data'))

  }
} # end function topHSPFwriteStatic()
