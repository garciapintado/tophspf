# Severn-Teme-Avon watershed
# run-off forecast ensemble parameters

RF <- list()
RF$path <- list()
RF$path$data <- dsndat
RF$region    <- region
RF$event     <- event
RF$scn       <- scn
RF$path$hpcscript <- file.path(dsndat, region, event, scn, 'input','RocksCluster_topHSPF.sh')
RF$modexe    <- 'hspf90_geoclaw'
RF$staTStr   <- "2007-05-01 00:00:00 GMT"   # May-Aug 2007 event
RF$endTStr   <- "2007-08-10 23:00:00 GMT"   #
RF$m         <- 150                         # ensemble size
RF$npc       <-  15                         # number of processors [coworkers] for HPC forecast
RF$hpc       <- 'essc'                      # %in% ['c2a','cca','essc','mac']
RF$mpi       <- TRUE                        # MPI for R-DA

RF$mcmode         <- 1                # 1 for MonteCarlo simulation | 0 for parameter inheritance
RF$mo             <- RF$m             # for pre-calibrated models (mcmode==0) this has to NULL
RF$inh            <- list()           # non-empty list for mcmode == 0
RF$inh$prior      <- file.path(PERM,'tophspf/data/severn_avon/200707/mc_oplo_500m') # previous simulation just used if RF$mcmode = 0
RF$inh$stat       <- 'E'
RF$inh$decreasing <- TRUE

RF$dto <- 3600                             # forcing data timestep [outer timestep model]
RF$dti <-  900                             # calculation timestep  [inner timestep model]
RF$simspdist <- 'spdis'                    # "spdis" for spatially distributed, "smdis" for spatially semidistributed
RF$rwave     <- 'k01'                      # 'routing wave for distributed model. swe': clawpack-SWE | 'k01': 1D kinematic routing

RF$maskf   <- 'mask.tif'                   # 1/0 for cells within/out of watershed. Geometadata is obtained from this map. path relative to regionp
RF$demf    <- 'demfilled_patch.tif'        # enhanced-connectivity   DEM [m]
RF$hrusf   <- 'hrus.tif'                   # HRU class map [0=out of domain]
RF$dproj   <- CRS('+init=epsg:27700')      # region projection. 27700 :: British National Grid, Ordnance Survery 1936

y <- list()
y$uzs   <- list()                                                   # surface soil moisture  [SWI ~ uzrat=usz/uzsn]
y$uzs$g <- list()                                                   # gauges
y$uzs$g$use <- FALSE
y$uzs$g$r      <- 0.10^2                                            # observation variance. Scalar
y$uzs$g$w      <- FALSE                                             # weighted observations?
y$uzs$g$qc     <- 'none'                                            # QC of SAR-based WLOs: oulier analysis of the innovations [0: noQC, < 0: global z-scores, > 0: radius for moving window z-scores, 'gstat': geostatistically obtained window]
y$uzs$g$twin   <- 60*60                                             # maximum observed time-window
y$uzs$g$nto    <- 5                                                 # maximum time observations
y$uzs$g$fname  <- ''

# note: all SSM datasets should be pre-merged into a common structure for DA, 
# pre-saved in a common .rds/.Rsav binary file, and loaded here as RF$y$uzs$data [see readOP.R]

y$uzs$s <- list()                                                     # EO SWI
y$uzs$s$use <- TRUE                                                   # do not assimilate but consider $data slot [if exists] for model output
y$uzs$s$r     <- 0.52                                                 # 10% grassland | scaled asar variance [note ascat is 0.025]
y$uzs$s$r     <- 0.19                                                 # 50% grassland | scaled asar variance [note ascat is 0.19]
y$uzs$s$w     <- FALSE
y$uzs$s$qc    <- 'none'
y$uzs$s$twin  <- Inf
y$uzs$s$nto   <- 10
y$uzs$s$fname <- ''                                                   # set in readOb.R

y$qou <- list() # global outflow from each catchment
y$qou$g <- list()
y$qou$g$use     <- FALSE
y$qou$g$r       <- 0.2^2                                                 # [m3/s]^2
y$qou$g$w       <- FALSE
y$qou$g$qc      <- 'none'
y$qou$g$twin    <- 60*60
y$qou$g$nto     <- 5
y$qou$g$fname   <- 'gaugeEA.Rsav'  # flow observations
y$qou$g$fmeta   <- 'outgrlets.txt' # metadata data.frame

RF$y <- y; rm(y)

anagen <- list() # generic analysis item
anagen$u     <- TRUE                    # updated by the analysis?
anagen$ll    <- 0.0                     # global filtering
anagen$infac <- 1.0                     # no inflation
anagen$min   <- NA                      # no lower bound
anagen$pos   <- c(0,0)                  # [0,0] => no localisation
anagen$trf   <- 'none'                  # no transformation

ana <- list()

ana$qou        <- anagen                # [m3/s?] catchment-scale outflows - just for augmentation 
ana$qou$ll     <- 20*km
ana$qou$infac  <- 0.0
ana$qou$min    <- 0.0

ana$ceps       <- anagen
ana$ceps$u     <- FALSE

ana$surs       <- anagen
ana$surs$ll    <- 150*km
ana$surs$infac <- 0.0
ana$surs$min   <- 0.0

ana$ifws       <- anagen
ana$ifws$ll    <- 150*km
ana$ifws$infac <- 0.0
ana$ifws$min   <- 0.0

ana$uzs        <- anagen
ana$uzs$ll     <- 150*km
ana$uzs$infac  <- 0.0
ana$uzs$min    <- 0.0
ana$uzs$trf    <- 'linearScale'

ana$lzs        <- anagen
ana$lzs$ll     <- 150*km
ana$lzs$infac  <- 0.0
ana$lzs$min    <- 0.0

ana$agws       <- anagen
ana$agws$ll    <- 150*km
ana$agws$infac <- 0.0
ana$agws$min   <- 0.0

ana$gwvs       <- anagen
ana$gwvs$ll    <- 150*km
ana$gwvs$infac <- 0.0
ana$gwvs$min   <- 0.0

ana$infilt       <- anagen
ana$infilt$infac <- 0.5
ana$infilt$min   <- 0.0025

ana$irc        <- anagen
ana$irc$infac  <- 0.5
ana$irc$min    <- 1.0E-30

ana$SGCfr       <- anagen                                                  # channel Manning's friction
ana$SGCfr$u     <- c(TRUE,TRUE,TRUE)
ana$SGCfr$infac <- 0.5                                                 # not NA to keep a minimum constant variance
ana$SGCfr$min   <- c(0.04,0.03,0.025)

ana$dsl         <- anagen      # downstream slope
ana$dsl$u       <- FALSE
ana$dsl$ll      <- 10*km
ana$dsl$infac   <- 0.0
ana$dsl$min     <- 3.0E-05
ana$dsl$pos     <- c(384537.5, 227887.5)

#for (i in 1:length(ana)) {
#  ana[[i]]$pos <- c(NA,NA)
#}

RF$ana <- ana; rm(ana)
