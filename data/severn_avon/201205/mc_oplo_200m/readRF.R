# Severn-Teme-Avon watershed
# run-off forecast ensemble parameters
RF <- list()
RF$m    <- 200                        # ensemble size
RF$np   <- 10                         # number of processors for parallel computing
RF$path <- list()
RF$path$hpcscript <- file.path(dsnsim,'RocksCluster_topHSPF.sh')
RF$path$data <- dsndat
RF$region    <- region
RF$event     <- event
RF$scn       <- scn

RF$mcmode    <- 1                     # 1 for MonteCarlo simulation | 2 for parameter inheritance
RF$mo        <- 60                    # or pre-calibrated models (mcmode==0) this has to NULL

regionou <- region  # this match is not required, just helps match input/output
eventou  <- event
