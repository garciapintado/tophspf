#!/bin/bash
dsndat='/glusterfs/land/users/jgp/tophspf'
region='severn_avon'
event='201205'
scn='mc_oplo_200m'
input='input'
mainfile='runhydro.run'

modexe='runhydro'

m=200
np=10

export PATH=$HOME/bin:$PATH
export LD_LIBRARY_PATH=$HOME/lib:$LD_LIBRARY_PATH

let "mp = m / np"
let "ini = ($SGE_TASK_ID - 1) * mp + 1"
let "end = ini + mp - 1"

for ((PID=$ini; PID<=$end; PID++))
do
  pidfmt=`printf "%05d" $PID`
 echo "$modexe $dsndat/$region/$event/$scn/$input/$pidfmt/$mainfile"
 $modexe $dsndat/$region/$event/$scn/$input/$pidfmt/$mainfile
done
