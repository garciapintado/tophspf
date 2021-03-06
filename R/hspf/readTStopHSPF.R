readTStopHSPF <- function(dsn = NULL, dto = NULL, dti = NULL, nto = NULL,
                          nhrus = NULL)
 # read timeseries of HRUs integrations generated by topHSPF
 #                    saveBinfl = 0,
 #                    savehdf5  = FALSE,
 #                    rmScratch = TRUE,
{

  innhrus  <- ifelse(nhrus > 1, nhrus+1, 1)           # number of input timeseries (RHUs plus one lumped output)

  if (dti > dto) {
   stop('topHSPFreadTS - ERR01') }
  if (dti == dto) {                                   # time series length
    stepdiv = 1.0
    nti = nto
  } else {
    stepdiv = dto / dti
    nti = nto * stepdiv
  }

  # binary timeseries generated by topHSPF
  outbints <- c("rain","petinp","ceps","surs","ifws","uzs","lzs","agws",
                "supy","suro","ifwo","agwo","pero","igwi",
                "cepe","uzet","lzet","agwet","baset","taet",
                "ifwi","uzi","infil","perc","lzi","agwi",
                "suri","rempet","gwvs","watin","watdif",
                "pers")
  modres <- list()

  tmptsf <- paste(outbints,"_hru_bin.ts",sep="")
  tmptsl <- paste(outbints,".mm",sep="")
  for(i in 1:length(outbints)){
    modres[[tmptsl[i]]] <- matrix(readBin(file.path(dsn,tmptsf[i]),
                             what="double", n=nti*innhrus),
                             nrow=nti, ncol=innhrus)
  }

  #save(modres, file=paste(outfld.ts,"/spsums_",runpidfmt,".Rdata",sep=""))       # re-export as an R structure

  #if(rmScratch){
  #  file.remove(c(fgeometap, fmaskp, fhrusp, frainp, fpetinpp,
  #               fparsp, modelfilep, mainfilep,
  #               file.path(outfld.ts.pid,tmptsf)
  #              ))
  #  file.remove(outfldp.pid)
  #  file.remove(outfld.ts.pid)
  #}
  return(modres)
}
