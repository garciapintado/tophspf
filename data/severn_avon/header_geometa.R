# Severn-Teme-Avon watershed
# CHARACTER      :: hrusf                "NA" for either lumped simulation or completely distributed simulation.
#                                        If available, the input map is a mask map of integers, where 0 is for the area outsize the watershed, and the indices within the watershed must be coded as the sequence 1:nhrus

RF$maskf <- 'mask.tif'                   # 1/0 for cells within/out of watershed. Geometadata is obtained from this map. path relative to regionp
RF$demf  <- 'demfilled_patch.tif'        # enhanced-connectivity   DEM [m]
RF$hrusf <- 'hrus.tif'                   # HRU class map [0=out of domain]
RF$dproj <- CRS('+init=epsg:27700')      # region projection. 27700 :: British National Grid, Ordnance Survery 1936
