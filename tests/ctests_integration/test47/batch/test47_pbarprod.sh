
#!/usr/bin/env bash
#

if [ "x" == "x$G4RELEASE" ]; then
echo " you must specify the Geang4 release/version; exit"
exit 3
fi

#source /products/setup
#setup root v5_34_01 -q "e2:prof"
source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh

# setup G4 datasets
#
#source /home/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup-${G4RELEASE}.sh
source /g4/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup-${G4RELEASE}.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${PBS_O_WORKDIR}/.. 
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable PBS_O_WORKDIR says: $PBS_O_WORKDIR"
fi

cd ${G4WORKDIR}

# JobID=1
JobID=${PBS_ARRAYID}
seed=$((1234+${JobID}))

#config=test47.proton.${target}.pbarprod
config=proton.${target}.${momentum}.pbarprod.${model}
if [ "x" == "x$PBS_ARRAYID" ]; then
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

if [ "x" == "x$PBS_ARRAYID" ]; then
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
