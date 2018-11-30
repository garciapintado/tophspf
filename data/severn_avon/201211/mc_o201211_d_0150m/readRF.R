# Severn-Teme-Avon watershed
# run-off forecast ensemble parameters
RF <- list()
RF$m      <- 150                         # ensemble size
RF$hpc    <- 'essc'                      # %in% ['c2a','cca','essc','mac']
RF$npc    <-  15                         # number of processors [coworkers] for parallel computing
RF$modexe <- 'hspf90_geoclaw'

#FF$dirroot        <- 'gcresults'                                              # relative path to Lisflood results, within .../scn
#FF$path$dirroot   <- file.path(dsndat, region, event, scn, FF$dirroot)        # absolute path to Lisflood simulation output
#FF$resroot        <- 'simple'                                                 # base naming of model input/output. Is is not required to match FF.region
#FF$srt   <- 3600*5                                                            # [s] maximum time allowed for jobs to run. Members exceeding this runtime will be killed
#FF$obs_variance_op  <- 0.25^2                                                 # observation variance for WL overpass points [m^2], i.e 0.25 m sdev
#FF$obs_variance_gau <- 0.10^2

RF$path <- list()
RF$path$hpcscript <- file.path(dsnsim,'input','RocksCluster_topHSPF.sh')
RF$path$data <- dsndat
RF$region    <- region
RF$event     <- event
RF$scn       <- scn

RF$mcmode         <- 1                # 1 for MonteCarlo simulation | 0 for parameter inheritance
RF$mo             <- RF$m             # for pre-calibrated models (mcmode==0) this has to NULL
RF$inh            <- list()           # non-empty list for mcmode == 0
RF$inh$prior      <- file.path(PERM,'tophspf/data/severn_avon/200707/mc_oplo_500m') # previous simulation just used if RF$mcmode = 0
RF$inh$stat       <- 'E'
RF$inh$decreasing <- TRUE

RF$simspdist <- 'spdis'                    # "spdis" for spatially distributed, "smdis" for spatially semidistributed
RF$rwave     <- 'k01'                      # 'routing wave for distributed model. swe': clawpack-SWE | 'k01': 1D kinematic routing

RF$mask    <- 'mask.tif'                   # 1/0 for cells within/out of watershed. Geometadata is obtained from this map. path relative to regionp
RF$dem     <- 'demfilled_patch.tif'        # enhanced-connectivity   DEM [m]
RF$hrus    <- 'hrus.tif'                   # HRU class map [0=out of domain]
RF$dproj   <- CRS('+init=epsg:27700')      # region projection. 27700 :: British National Grid, Ordnance Survery 1936

RF$qhobf     <- 'gaugeEA.Rsav'             # gauge flow and h observations. '' for NA
RF$gaugesf   <- 'outgrlets.txt'            # locations for timeseries output sampling

y <- list()
y$surs <- list()
y$surs$g <- list()
y$surs$g$use <- FALSE
y$surs$g$r       <- 0.10^2                                            # observation variance. Scalar
y$surs$g$qc      <- 'none'                                            # QC of SAR-based WLOs: oulier analysis of the innovations [0: noQC, < 0: global z-scores, > 0: radius for moving window z-scores, 'gstat': geostatistically obtained window]
y$surs$g$twin    <- 60*60                                             # maximum observed time-window
y$surs$g$nto     <- 5                                                 # maximum time observations
y$surs$g$fname   <- ''

y$uzs <- list()                                                     # mapping of surface soil moisture
y$uzs$g <- list()
y$uzs$g$use <- FALSE
y$uzs$g$r      <- 0.10^2                                            # observation variance. Scalar
y$uzs$g$qc     <- 'none'                                            # QC of SAR-based WLOs: oulier analysis of the innovations [0: noQC, < 0: global z-scores, > 0: radius for moving window z-scores, 'gstat': geostatistically obtained window]
y$uzs$g$twin   <- 60*60                                             # maximum observed time-window
y$uzs$g$nto    <- 5                                                 # maximum time observations
y$uzs$g$fname  <- ''

# note: all SSM datasets should be pre-merged into a common structure for DA, 
# pre-saved in a common .rds/.Rsav binary file, and loaded here as RF$y$uzs$data [see readOP.R]
y$uzs$s <- list()
y$uzs$s$use <- FALSE                                                # do not assimilate but consider $data slot [if exists] for model output
y$uzs$s$r     <- 2.0E-03^2                                          # [m^2] variance [2mm sd]
y$uzs$s$qc    <- 'none'
y$uzs$s$twin  <- Inf
y$uzs$s$nto   <- 10
y$uzs$s$fname <- ''
RF$y <- y; rm(y)

ana <- list()
ana$uzs <- list()
ana$uzs$u  <- TRUE
ana$uzs$ll <- 20.E03
ana$uzs$infac <- 0.5
ana$uzs$min   <- 0.0
ana$lzs <- list()
ana$lzs$u  <- TRUE
ana$lzs$ll <- 20.E03
ana$lzs$infac <- 0.5
ana$lzs$min   <- 0.0
ana$ifws <- list()
ana$ifws$u  <- TRUE
ana$ifws$ll <- 20.E03
ana$ifws$infac <- 0.5
ana$ifws$min   <- 0.0
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
ana$dsl$pos     <- c(384537.5, 227887.5)
RF$ana <- ana; rm(ana)
