#!/usr/bin/env bash
#

if [ "x" == "x$G4RELEASE" ]; then
echo " you must specify the Geang4 release/version; exit"
exit 3
fi

source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh

# setup G4 datasets
#
source /g4/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup-${G4RELEASE}.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
# ---> PBS/obsolete ---> G4WORKDIR=${PBS_O_WORKDIR}/.. 
G4WORKDIR=${SLURM_SUBMIT_DIR}/..
else
    echo "Variable says: $G4WORKDIR"
    # ---> PBS/obsolete ---> echo "Variable PBS_O_WORKDIR says: $PBS_O_WORKDIR"
    echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"
fi

cd ${G4WORKDIR}

# ---> leftover from using pbs_multiCore_master.sh 
# ---> now passed in as an explicit env.var. ---> nevents=${1}

# JobID=1
#JobID=${PBS_ARRAYID}
# ---> leftover from using pbs_multiCore_master.sh ---> JobID=${2}
JobID=${SLURM_ARRAY_TASK_ID}
seed=$((1234+${JobID}))

echo " JobID = $JobID"
echo " seed = $seed"
echo " nevents = $nevents"


#config=test47.proton.${target}.pbarprod
config=proton.${target}.${momentum}.pbarprod.${model}
if [ "x" == "x$JobID" ]; then
   echo "Process entire statisticas in one job"
else
config=proton.${target}.${momentum}.pbarprod.${model}.${JobID}
fi

if [ -e ${config} ]; then
/bin/rm -f ${config}
fi

printf "#verbose \n-1 \n#rad \n" >> ${config}
printf "#events \n" >> ${config}
printf "%d\n" ${nevents} >> ${config}

if [ "x" == "x$JobID" ]; then
   echo "skip seed and jobid settings"
else
printf "#randomSeed\n" >> ${config}
printf "%d\n" ${seed}  >> ${config}
printf "#jobID\n" >> ${config}
printf "%d\n" ${JobID} >> ${config}
fi

printf "//--------proton_processes \n" >> ${config}
printf "#particle \nproton \n#isPbarProd \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n" >> ${config}

printf "//-------- \n" >> ${config}

printf "#material \n" >> ${config}
printf "%s\n" ${target} >> ${config}

printf "#momentum(MeV/c) \n" >> ${config}
#printf "9956.0 \n">> ${config}
printf "%d\n" ${momentum} >> ${config}

printf "// --- \n#generator \n" >> ${config}
# printf "%s\n" ftfp >> ${config}
printf "%s\n" ${model} >> ${config}
printf "#run \n" >> ${config}

printf "#exit\n" >> ${config}

./test47 ${config}

if [ -e ${config} ]; then
/bin/rm -f ${config}
fi
