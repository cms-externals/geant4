#!/usr/bin/env bash
#

source /products/setup
setup root v5_34_01 -q "e2:prof"

# setup G4 datasets
#
source /home/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup.sh

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

printf "#verbose\n0\n#rad\n" >> ITEP.proton.${JobID}
printf "#events \n31250 \n" >> ITEP.proton.${JobID}
printf "#randomSeed\n" >> ITEP.proton.${JobID}
printf "%d\n" ${seed}  >> ITEP.proton.${JobID}
printf "#jobID\n" >> ITEP.proton.${JobID}
printf "%d\n" ${JobID} >> ITEP.proton.${JobID}
printf "//--------Proton_processes \n" >> ITEP.proton.${JobID}
printf "#particle\nproton\n#isITEP\n#position(mm)\n0. 0. 0.\n#direction\n0. 0. 1.\n//-------- \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#material \nC \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n7500.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#material \nU \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n1400.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nbinary \n#run \n" >> ITEP.proton.${JobID}

printf "//-------- \n" >> ITEP.proton.${JobID}
printf "#momentum(MeV/c) \n7500.0 \n" >> ITEP.proton.${JobID}
printf "// --- \n#generator \nbertini \n#run \n// --- \n#generator \nftfp \n#run \n" >> ITEP.proton.${JobID}

printf "#exit\n" >> ITEP.proton.${JobID}


${G4EXE}/test47 ITEP.proton.${JobID}
