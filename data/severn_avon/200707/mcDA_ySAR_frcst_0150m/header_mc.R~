# parameter & initial condition input for the Severn-Teme-Avon watersheds and topHSPF model
# ------------------------------------------------------------------
#ldd    <- "ldd.map"
#drofac <- 0.2                   # direct runoff factor            [0/1]
#irofac <- 0.5                   # indirect runoff factor          [0/1]
#sodept.mul  <- 1                # multiplier of estimated soil depth
#wsoil.fc.01 <- 0.35             # soil-moisture at field capacity [0/1]
#wsoil.wp.01 <- 0.05             # soil-moisture at wilting point  [0/1]
#wsoil.in.01 <- 0.20             # initial soil-moisture           [0/1]
#ksat.mul    <- 1                # multiplier of ksat.map

# note: lsur is later estimated, for 'smdis' as: lsur <- (sqrt(summary(as.factor(hrus@data[,1]))*1000000)/2)[-1]
# parameters to model (all subject to calibration)

# global parameters for topHSPF
GP <- list()
GP$lzsn   <- c(.25,25.0)         # [mm]    [.25,2500.] lower zone nominal storage                                                           PWAT-PARM2
GP$infilt <- c(0.1,0.5)          # [mm/hr] [.0025,2500.] index to the infiltration capacity of the soil                                           "
GP$lsur   <- 0                   # [m]     [.3,none] length of assumed overland flow plane " (specified below for every subwatershed)
GP$slsur <- "slpfilled.tif"      # [0/1] [.000001,10.] distributed slope of the overland flow plane                . path relative to regionp
GP$kvary  <-   0.0               # [1/mm]  [.0,none] parameter which affects the behavior of groundwater recession flow, enabling
                                 #                   it to be non-exponential in its decay with time                                              "
GP$agwr   <-  c(0.6,0.9)         # [1/day] [.001-0.999] basic groundwater recession rate if KVARY is zero and there is no inflow to groundwater   "
GP$infexp <-  c(0.1,4.0)         # [-]     2.0 is default [.0,10.] exponent of the infiltration equation                                                            PWAT-PARM3
GP$infild <-    2.0              # [-]     2.0 is default [1.,2.]  ratio between the maximum and mean infiltration capacities over the PLS                                "
GP$deepfr <-    0.0              # [-]     0.0 is default [.0,1.] fraction of groundwater inflow which will enter deep (inactive) groundwater, and, thus, be lost from the system as defined in HSPF
GP$basetp <-    0.0              # [-]     0.0 is default [.0,1.] fraction of remaining potential ET which can be satisfied from baseflow (groundwater outflow), if enough is available
GP$agwetp <-    0.0              # [-]     0.0 is default [.0,1.] fraction of remaining potential ET which can be satisfied from active groundwater storage if enough is available
GP$cepsc  <- c(0.0,0.4)          # [mm]    0.0 is default [.0,250.] maximum interception capacity                                PWAT-PARM4
GP$uzsn   <- c(.35,1.0)          # [mm]    none is default [.25,250.] upper soil zone nominal storage (hrus's distributed)             "
GP$nsur   <- c(0.06,0.1)         # [cplx]  0.1 is default [0.001,1.0] manning coefficient for overland flow                            "
GP$intfw  <- c(3.0,8.0)          # [-]     none is default [0.0,none] interflow inflow parameter                                       "
GP$irc    <- c(0.70,0.75)        # [1/day] none is default [1.0E-30,0.999] interflow recession parameter. Under zero inflow,
                                 #         this is the ratio of today's interflow outflow rate to yesterday's rate               PWAT-PARM4
GP$lzetp  <- c(0.0,0.5)          # [-]     0.0 is default [0.0,1.0] lower zone ET parameter;
                                 #         it is an index to the density of deep-rooted vegetation                               PWAT-PARM4
# initial conditions [HSPF: PWAT-STATE1] here embedded in GP
GP$ceps <- c(0.15,0.4)   # [mm] interception storage (in [0,cepsc])          PWAT-STATE1
GP$surs <- c(.025,2.0)   # [mm] surface detention storage                    PWAT-STATE1
GP$ifws <- c(0.0,1.0)    # [mm] interflow storage                            PWAT-STATE1
GP$uzs  <- c(.25,1.0)    # [mm] upper zone storage                           PWAT-STATE1
GP$lzs  <- c(.25,25.0)   # [mm] lower zone storage                           PWAT-STATE1
GP$agws <- c(2.0,5.0)    # [mm] active groundwater storage                   PWAT-STATE1
GP$gwvs <- c(0.0,5.0)    # [mm] groundwater variable storage index           PWAT-STATE1

MCflags <- rep(FALSE,length(GP))
names(MCflags) <- names(GP)
MCflags[c('lzsn','infilt','agwr','infexp','cepsc','uzsn','nsur','intfw','irc','lzetp',
          'ceps','surs','ifws','uzs','lzs','agws','gwvs')] <- TRUE   # parameters & initial conditions to randomize for MC

# subwatershed specific parameters
npar <- length(GP)
hruSP <- readGDAL(file.path(dsnscn,RF$hrus))
nhrus <- max(hruSP@data[,1], na.rm=TRUE)
dx <- hruSP@grid@cellsize[1]

WP <- list()                                        # specific subwatershed data
swareas <- rep(0,nhrus)
for (iw in 1:nhrus) {
  swareas[iw] <- sum(hruSP@data[,1]==iw)*dx^2    # [m2] subwatershed areas
  WP[[iw]] <- list()
  for (ipar in 1:npar) {
    WP[[iw]] <- GP                                  # inherit from global parameters
  }
  WP[[iw]]$lsur <- sqrt(swareas[iw]) / 2            # [m] rough 'half-book' estimate for 'smdis'
}
RF$nhrus <- nhrus
