#!/usr/bin/env bash
#

#if [ "x" == "x$target" ] ; then
# bail out if NO TARGET !!!
#fi
#if [ "x" == "x$momentum" ] ; then
## bail out if NO MOMENTUM !!!
#fi

source /cvmfs/geant4-ib.opensciencegrid.org/products/setup
setup root v6_24_06b -q e20:p399:prof

# setup G4 datasets
#
source /wclustre/g4p/g4p/download/g4data/g4-datasets-${G4RELEASE}.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${SLURM_SUBMIT_DIR}/..
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"
fi

cd ${G4WORKDIR}

# ---> leftover from using pbs_multiCore_master.sh; 
# ---> now passed in as an explicit env.var. ---> nevents=${1}

# JobID=1
JobID=$((1+${SLURM_PROCID}))
seed=$((1234+${JobID}))

bm=${beam}
if [ "x${beam}" == "xpiplus" ]; then
   bm="pi+"
fi
if [ "x${beam}" == "xpiminus" ]; then
   bm="pi-"
fi

config=${beam}.${target}.${momentum}.HARP.${physlist}

if [ "x" == "x$JobID" ]; then
   echo "Process entire statistics in one job"
else
 config=${beam}.${target}.${momentum}.HARP.${physlist}.${JobID}
fi

if [ -e ${G4WORKDIR}/${config} ]; then
/bin/rm -f ${G4WORKDIR}/${config}
fi

/usr/bin/printf "#verbose \n-1 \n#rad \n" >> ${config}
/usr/bin/printf "#events \n" >> ${config}
/usr/bin/printf "%d\n" ${nevents} >> ${config}
if [ "x" == "x$JobID" ]; then
   echo "skip seed and jobid settings"
else
/usr/bin/printf "#randomSeed\n" >> ${config}
/usr/bin/printf "%d\n" ${seed}  >> ${config}
/usr/bin/printf "#jobID\n" >> ${config}
/usr/bin/printf "%d\n" ${JobID} >> ${config}
fi
/usr/bin/printf "//-------- \n" >> ${config}
/usr/bin/printf "#particle \n" >> ${config}
/usr/bin/printf "%s\n" ${bm} >> ${config}
/usr/bin/printf "#is%s\n" ${EXP} >> ${config}
/usr/bin/printf "#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> ${config}

/usr/bin/printf "//-------- \n" >> ${config}

/usr/bin/printf "#target-geom(mm) \n0.0 3.15 160. G4Tubs \n" >> ${config}
/usr/bin/printf "#material \n" >> ${config}
/usr/bin/printf "%s\n" ${target} >> ${config}

/usr/bin/printf "#momentum(MeV/c) \n" >> ${config}
/usr/bin/printf "%d\n" ${momentum} >> ${config}

/usr/bin/printf "// --- \n#physicslist \n" >> ${config}
/usr/bin/printf "%s\n" ${physlist} >> ${config}
/usr/bin/printf "#run \n" >> ${config}

/usr/bin/printf "#exit\n" >> ${config}


./test23 ${config}

if [ -e ${config} ]; then
/bin/rm -f ${config}
fi
