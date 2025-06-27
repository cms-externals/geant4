#!/usr/bin/env bash
#

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

config=muminus.${target}

# JobID=1
#JobID=${PBS_ARRAYID}  # ---> would now be: SLURM_ARRAY_TASK_ID
#seed=$((1234+${JobID}))


printf "#verbose \n-1 \n#rad \n" >> ${config}
printf "#events \n1000000 \n" >> ${config}

#printf "#randomSeed\n" >> ${config}
#printf "%d\n" ${seed}  >> ${config}
#printf "#jobID\n" >> ${config}
#printf "%d\n" ${JobID} >> ${config}

printf "//--------Muminus_processes \n" >> ${config}
printf "#particle \nmu- \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n//-------- \n" >> ${config}

printf "#material \n${target} \n" >> ${config}

printf "#momentum(MeV/c) \n0.0 \n" >> ${config}

printf "// --- \n" >> ${config}

# --> printf "generator \nstopping \n#run \n" >> ${config}
printf "#generator \ncaptureUpdate \n#run \n" >> ${config}

printf "#exit\n" >> ${config}


#./test48 test48.muminus
./test48 ${config}

$ROOTSYS/bin/root -b -p -q MuMinusModels.C\(\"${target}\"\)  
 
if [ "x" == "x$G4RELEASE" ] ; then
    echo "Variable G4RELEASE is not set"
else
. batch/copy_results.sh ${G4RELEASE}
fi

/bin/rm -r ${config}
