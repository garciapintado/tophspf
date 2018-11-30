# readprm.R (filter parameters)
#
#readprm <- function(){
prm <- list()
prm$method  <- 'ETKF'                                                 # DA method
prm$rfactor <- 1                                                      # multiple (scaling) for R (observation error covariance)
prm$loc_boo      <- TRUE                                              # whether to conduct localization 0 or 1 [apply to all elements in the augmented state vector]
prm$loc_method   <- 'LA'                                              # localization method. 'CF' (covariance filtering), or 'LA' (local analysis)
prm$loc_function <- 'Gaspari_Cohn'                                    # tag for the localisation function (calc_loccoeffs.m)
prm$loc_distype  <- '2D'                                              # '2D' for Euclidean, or 'SG' for along network distances [apply to all elements in the augmented state vector]
prm$sgf          <- 'sgriv.rds'                                       # SpatialGraph representing the river network. Just required if loc_distype=='SG'
prm$gridSGf      <- 'GgridSG.rds'                                     # precalculated grid-to-SpatialGraph connectivity object [SpatialPointsDataFrame]
prm$rotate       <- FALSE                                             # 0 = do not rotate. 1 = rotate. The code has to prepared to set prm.rotate each some assimilation steps
prm$rotate_ampl  <- 1                                                 # always 1. Code not prepared for other values
prm$HnnDEM       <- FALSE                                             # whether to use the downstream-transform to obtain HE
prm$downstrmapf <- 'downstrfilled.asc'                                # downstream map. Just required if HnnDEM == 1 
#return(prm)
#}
