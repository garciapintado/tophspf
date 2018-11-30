rainfp   <- file.path(dsnsim,'rain_pointer.txt')                      # ASCII file with pointers to rainfall files [mm/tstep]
dsnrain  <- file.path(PERM,'rain','EA_gridded_SA','2002_2012_01h_1000m','sri')

dsnpe    <- file.path(HOME,'docs/DEMON/data/pe/morecs/2003_2007')
peRsav   <- 'morecs_hourly_2007.Rsav'

# ---

load(file.path(dsnpe,peRsav))
pe <- pehlst[['135']]; rm(pehlst)
pe$value <- round(pmax(pe$value[,2],0.1E-10),10) # [mm/h] 2 corresponds to potential evapotranspiration according to MORECS

