makehpcScript <- function(hpcscript, dsndat, region, event, scn,
                          input, mainfile, modexe, m, np) {
  # +++ purpose +++
  # write the .sh script to run this specific simulation in parallel computing [tophspf]
  # ---------
  tmpf <- hpcscript
  zo <- file(tmpf,'w')
  cat("#!/bin/bash\n",                         file=zo)
  cat("dsndat='",dsndat,        "'\n", sep="", file=zo)
  cat("region='",region,        "'\n", sep="", file=zo)
  cat("event='",event,          "'\n", sep="", file=zo)
  cat("scn='",scn,              "'\n", sep="", file=zo)
  cat("input='",input,          "'\n", sep="", file=zo)
  cat("mainfile='",mainfile,    "'\n", sep="", file=zo)
  cat("\n",                                    file=zo)
  cat("modexe='",modexe,        "'\n", sep="", file=zo)
  cat("\n",                                    file=zo)
  cat("m=",m,                   "\n",  sep="", file=zo)
  cat("np=",np,                 "\n",  sep="", file=zo)
  cat("\n",                                    file=zo)
  cat("ls $dsndat\n",           file=zo)
  cat("ls $dsndat/$region/$event/$scn\n",      file=zo)
  cat('let "slptime = 5+$SGE_TASK_ID"\n',      file=zo)
  cat("sleep $slptime\n",                      file=zo)
  cat("\n",                                    file=zo)
  cat("cwd=$PWD\n",                            file=zo)
  cat("\n",                                    file=zo)
  cat("export PATH=$HOME/bin:$PATH\n",         file=zo)
  cat("export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH\n", file=zo)
  cat("\n",                                    file=zo)
  cat('let "mp = m / np"\n',
      'let "ini = ($SGE_TASK_ID - 1) * mp + 1"\n',
      'let "end = ini + mp - 1"\n\n',
      'for ((PID=$ini; PID<=$end; PID++))\n',
      'do\n',
      '  x=`printf "%05d" $PID`\n',
      '  cd $dsndat/$region/$event/$scn/$input/$x\n',
      '  echo "$modexe $dsndat/$region/$event/$scn/$input/$x/$mainfile"\n',
      '  $modexe $dsndat/$region/$event/$scn/$input/$x/$mainfile\n',
      'done\n', sep='', file=zo)
  cat('cd $cwd\n',                             file=zo)
  close(zo)
} # end function makehpcScript
