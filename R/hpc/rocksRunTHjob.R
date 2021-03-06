rocksRunTHjob <- function(RF, hrt='05') {
 # hrt   : [hh] string with maximum running time for job before being deleted

 pauset <- 60                # checks for job completion every pauset processing time [s]

 m  <- RF$m
 np <- RF$np

 # re-assure previous run output error code file is deleted
 for (im in 1:RF$m) {
   imStr <- formatC(im,width=5,flag="0",format="d")   
   fname <- file.path(RF$path$outsfld,imStr,"topHSPF.err")
   if (file.exists(fname))
     file.remove(fname)
 }
 
 # syscmd <- paste('qsub -cwd -l batch,h_cpu=2:00:0 -t 1-',np,' ', RF$path$hpcscript, sep='')
 syscmd <- paste('qsub -cwd -P sinatra -l project -l h_rt=',hrt,':00:00 -t 1-',np,' ', RF$path$hpcscript, sep='')       # SGE - SINATRA nodes

 system(syscmd)

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
 return(0)
}
