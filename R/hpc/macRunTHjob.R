macRunTHjob <- function(RF) {
 ct <- Sys.time()

 dsnsim <- file.path(RF$path$data,RF$region,RF$event,RF$scn)

 cwd <- getwd()
 setwd(dsnsim)

 pauset <- 10       # checks for job completion every pauset processing time [s]

 m  <- RF$m
 np <- RF$np
 if (RF$simspdist == 'smdis') {
   exec <- 'hspf'                # just for lumped model
 } else {
   exec <- 'hspf90_geoclaw'      # includes SWE / K01 for overland flow
 }

 # serial running: TODO adapt F90 to make it able to run with I files from everywhere
 for (im in 1:m) {
   cat('running im: ',im,'\n')
   imStr <- formatC(im,width=5,flag='0',format='d')
   syscmd <- paste(exec,file.path(dsnsim,'input',imStr,'runhydro.run'))
   dsnim <- file.path(dsnsim,'input',imStr)
   setwd(dsnim)
   system(syscmd, wait=TRUE)
 }

 done <- FALSE
 waittime <- 0
 while (!done) {
   cat('waittime =',waittime * pauset / 3600,'hours\n')
   waittime = waittime + 1
   Sys.sleep(pauset)
   done <- TRUE
   for (im in 1:m) {
     imStr <- formatC(im,width=5,flag="0",format="d")
     if (!file.exists(file.path(RF$path$outsfld,imStr,"topHSPF.err")))
       done <- FALSE
   }
 } # end while
 setwd(cwd)
 return(Sys.time() - ct)
}
