writeGmeta <- function(G, fname) {
  # write G meta to an ASCII file to be imported by Fortran
  zo  <- file(fname,'w')
  cat(G$LOCATION,          '\n', file=zo)
  cat(G$MAPSET,            '\n', file=zo)
  cat(G$proj4@projargs,    '\n', file=zo)
  cat(G$n,G$s,G$w,G$e,     '\n', file=zo)
  cat(G$nsres,G$ewres,     '\n', file=zo)
  cat(G$rows,G$cols,       '\n', file=zo)
  cat(G$cells,             '\n', file=zo)
  cat(G$xlims,             '\n', file=zo)
  cat(G$ylims,             '\n', file=zo)
  cat(G$xseq,              '\n', file=zo)
  cat(G$yseq,              '\n', file=zo)
  cat(G$ryseq,             '\n', file=zo)
  close(zo)
}
