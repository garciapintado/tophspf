# Severn-Teme-Avon watershed
# run-off forecast ensemble parameters
RF <- list()
RF$m    <- 150                        # ensemble size
RF$np   <- 5                         # number of processors for parallel computing
RF$path <- list()
RF$path$hpcscript <- file.path(dsnsim,'RocksCluster_topHSPF.sh')
RF$path$data <- dsndat
RF$region    <- region
RF$event     <- event
RF$scn       <- scn

RF$mcmode    <- 0                     # 1 for MonteCarlo simulation | 0 for parameter inheritance
RF$mo        <- 150                   # for pre-calibrated models NULL or equal to m [forecasting mode (mcmode==0)]

RF$inh            <- list()                                                                      # non-empty list for mcmode == 0
RF$inh$prior      <- file.path(PERM,'tophspf/data/severn_avon/200707/mc_oplo_500m')              # previous simulation, required if RF$mcmode = 0, useless otherwise
RF$inh$stat       <- 'E'
RF$inh$decreasing <- TRUE


regionou <- region  # this match is not required, just helps match input/output
eventou  <- event
