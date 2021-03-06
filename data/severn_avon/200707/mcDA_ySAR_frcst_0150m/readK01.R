# parameter for 1D cascade kinematic wave
K01 <- list()

K01$dem    <- 'demfilled.tif'       # K01-compliant DTM
                                    #                                7 8 9
K01$ldd    <- 'lddfilled.tif'       # [-] local drainage directions  4 5 6
                                    #                                1 2 3
K01$aflx   <- 'accufluxfilled.tif'  # [-] accumulated [pixel number] influx into cells
K01$cmsk   <- 'cmsk.tif'      # channel classes [for friction]
K01$cwid   <- 'cwid.tif'      # channel width [raw powerlaw estimate based on accumulated flux [makeChannels.R]
K01$cdep   <- 'cdep.tif'      # channel depth [ "    " ...]

K01$SGCmc  <- 'SA_SGCmc.pram' # friction parameters
K01$dsl    <- c(6E-05, 1.5E-06) # downstream boundary condition - bathymetric slope S_0 (mean, sd)
#K01$slo0   <- 'SA_slo0.bin'   # bathymetric-based slope

