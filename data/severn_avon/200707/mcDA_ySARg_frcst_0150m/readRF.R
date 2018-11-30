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
y$uzs$g$qc     <- 'none'                                            # QC of SAR-based WLOs: oulier analysis of the innovations [0: noQC, < 0: global z-scores, > 0: radius for moving window z-scores, 'gstat': geostatistically obtained window]
y$uzs$g$twin   <- 60*60                                             # maximum observed time-window
y$uzs$g$nto    <- 5                                                 # maximum time observations
y$uzs$g$fname  <- ''

# note: all SSM datasets should be pre-merged into a common structure for DA, 
# pre-saved in a common .rds/.Rsav binary file, and loaded here as RF$y$uzs$data [see readOP.R]

y$uzs$s <- list()                                                   # EO SWI
y$uzs$s$use <- TRUE                                                 # do not assimilate but consider $data slot [if exists] for model output
y$uzs$s$r     <- 0.52                                               # 10% grassland | scaled asar variance [note ascat is 0.025]
y$uzs$s$r     <- 0.41                                               # 50% grassland | scaled asar variance [note ascat is 0.19]
y$uzs$s$qc    <- 'none'
y$uzs$s$twin  <- Inf
y$uzs$s$nto   <- 10
y$uzs$s$fname <- ''

y$qou <- list() # global outflow from each catchment
y$qou$g <- list()
y$qou$g$use     <- TRUE
y$qou$g$r       <- 0.2^2                                                 # [m3/s]^2
y$qou$g$qc      <- 'none'
y$qou$g$twin    <- 60*60
y$qou$g$nto     <- 5
y$qou$g$fname   <- 'gaugeEA.Rsav'  # flow observations
y$qou$g$fmeta   <- 'outgrlets.txt' # metadata data.frame

RF$y <- y; rm(y)

ana <- list()

ana$qou    <- list() # catchment-scale outflows - just for augmentation 
ana$qou$u      <- TRUE  # dummy
ana$qou$ll     <- 20.E03
ana$qou$infac  <- 0.0
ana$qou$min    <- 0.0

ana$ceps <- list()
ana$ceps$u     <- FALSE
ana$ceps$ll    <- 20.E03
ana$ceps$infac <- 0.5
ana$ceps$min   <- 0.0

ana$surs <- list()
ana$surs$u     <- FALSE
ana$surs$ll    <- 20.E03
ana$surs$infac <- 0.5
ana$surs$min   <- 0.0

ana$ifws <- list()
ana$ifws$u     <- TRUE
ana$ifws$ll    <- 20.E03
ana$ifws$infac <- 0.5
ana$ifws$min   <- 0.0

ana$uzs <- list()
ana$uzs$u      <- TRUE
ana$uzs$ll     <- 20.E03
ana$uzs$infac  <- 0.5
ana$uzs$min    <- 0.0

ana$lzs <- list()
ana$lzs$u     <- TRUE
ana$lzs$ll    <- 20.E03
ana$lzs$infac <- 0.5
ana$lzs$min   <- 0.0

ana$agws <- list()
ana$agws$u     <- TRUE
ana$agws$ll    <- 20.E03
ana$agws$infac <- 0.5
ana$agws$min   <- 0.0

ana$gwvs <- list()
ana$gwvs$u     <- TRUE
ana$gwvs$ll    <- 20.E03
ana$gwvs$infac <- 0.5
ana$gwvs$min   <- 0.0

ana$SGCfr <- list()                                                     # channel Manning's friction
ana$SGCfr$u      <- rep(TRUE,3)
ana$SGCfr$ll     <- 0.0
ana$SGCfr$infac  <- 1.0                                                 # not NA to keep a minimum constant variance
ana$SGCfr$min    <- c(0.04,0.03,0.025)

ana$dsl <- list()      # downstream slope
ana$dsl$u     <- TRUE
ana$dsl$ll    <- 10.E03
ana$dsl$infac <- 0.0
ana$dsl$min   <- 3.0E-05

for (i in 1:length(ana)) {
  ana[[i]]$pos <- c(NA,NA)
  ana[[i]]$trf <- 'none'
}
ana$uzs$trf <- 'linearScale'

ana$dsl$pos     <- c(384537.5, 227887.5)
RF$ana <- ana; rm(ana)
