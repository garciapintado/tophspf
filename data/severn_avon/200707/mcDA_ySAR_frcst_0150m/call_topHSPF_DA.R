# call to topHSPFmcf
#
# launcher to conduct a forward MC analysis for calibrating topHSPF at the Severn-Avon
#
# non-interactive execution:
# R CMD BATCH --no-save --no-restore-data --no-readline $HOME/docs/hydrology/tophspf/data/severn_avon/200707/mcDA_ySA/call_topHSPFmcf.R &


 startT <- Sys.time()

 HOME <- Sys.getenv('HOME')
 PERM <- Sys.getenv('NCEO')
 HSPF <- file.path(HOME,'docs','hydrology','tophspf')          # source data & code files
 #HSPF <- Sys.getenv("HSPF")

 dsnDAR <- file.path(HOME,'docs','DEMON','DA','R')
 
 dsn    <- file.path(HSPF,'data')
 dsndat <- file.path(PERM,'tophspf','data')       # PERMANENT [non-backed-up]
 region <- 'severn_avon'
 event  <- '200707'
 scn    <- 'mcDA_ySAR_frcst_0150m'

 dsnqhob <- file.path(HOME,'docs','DEMON','data','flow_stage',event,'flow_level')
 qhobf   <- 'gaugeEA.Rsav'

 dsnscn <- file.path(dsn,    region, event, scn)  # source data in backed-up volume
 dsnsim <- file.path(dsndat, region, event, scn)  # simulation data in non-backed-up volume

 if (file.exists(dsnsim))                        # create IO folders
   system(paste("rm -r -f",dsnsim))
 dir.create(dsnsim, recursive=TRUE)
 dir.create(file.path(dsnsim,'input'))
 dir.create(file.path(dsnsim,'output'))

 # link header / files
 system(paste("ln -sf", file.path(dsn,region,'header_geometa.R'),     dsnscn))
 system(paste("ln -sf", file.path(dsn,region,event,'header_event.R'), dsnscn))
 system(paste("ln -sf", file.path(dsn,region, 'watmaps','*'),         dsnscn))
 system(paste("ln -sf", file.path(dsnqhob,'R',qhobf),                 dsnscn))   # observations file [flow & level time series]
 #system(paste("ln -sf", file.path(dsnscn,'low_sevn_sub.stage'),       dsnscn))

 rivSPf <- list()
 rivSPf$dsn   <- file.path(PERM,'gis','shareGeo','GB_Rivers')
 rivSPf$layer <- 'gb_rivers' # 'ls_rivers' also shows the major rivers for the Severn-Avon
 prjCRS       <- '+init=epsg:27700'

 source(file.path(HSPF, 'R', 'topHSPF_DA.R'))
 #err <- topHSPF_DA(dsn, dsndat, region, event, scn, dsnDAR,
 #                  HSPF, prjCRS, rivSPf)
 #cat('error:',err,'\n')

if (1 > 2) {
  cwd <- getwd()
  setwd(dsnsim)
  load('MCres.Rsav')
  load('MCglopar.Rsav')
  load('glike.Rsav')
  setwd(cwd)
}
