#!/usr/bin/env bash
#
# Once submitted to batch queue via SLURM, this script will check 
# the number of cores per worker node (e.g. 32 for an amd32 node)
# and will launch equivalent number of parallel jobs, splitting 
# the total statistics evenly among those jobs
#
# Usage example:
# sbatch -p amd32_g4val_slow -N1 -t 24:00:00 \
#        --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=NuBeam,nevtotal=1000000" \
#        -A g4p proton_NA61_SLURM.sh
# 
# Alternatively, one can use amd32, or intel12 queue, 
# or even amd32_g4perf with the --qos=g4perf for fastest processing,
# but if the later, keep in mind 24hrs limit 

g4had_run() {

JobID=$1
nevents=$2

source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh
source /g4/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup-${G4RELEASE}.sh

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${SLURM_SUBMIT_DIR}/..
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"
fi

echo " G4WORKDIR = ${G4WORKDIR}"

cd ${G4WORKDIR}
echo " We are in ${PWD}" 

seed=$((1234+${JobID}))
echo " random seed = ${seed}"

config=proton.${target}.${momentum}.NA61.${physlist}.${JobID}

if [ -e ${G4WORKDIR}/${config} ]; then
/bin/rm -f ${G4WORKDIR}/${config}
fi

/usr/bin/printf "#verbose \n-1 \n#rad \n" >> ${config}
/usr/bin/printf "#events \n" >> ${config}
/usr/bin/printf "%d\n" ${nevents} >> ${config}
/usr/bin/printf "#randomSeed\n" >> ${config}
/usr/bin/printf "%d\n" ${seed}  >> ${config}
/usr/bin/printf "#jobID\n" >> ${config}
/usr/bin/printf "%d\n" ${JobID} >> ${config}
/usr/bin/printf "//--------Proton_processes \n" >> ${config}
/usr/bin/printf "#particle \nproton \n#isNA61 \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> ${config}

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

}

coresPerNode=`cat /proc/cpuinfo | grep -c processor`
nodeName=`uname -n`
echo "number of cores on node ${nodeName} = ${coresPerNode}"

nevt=$(( nevtotal/coresPerNode ))

for JobID in $(seq ${coresPerNode}); do

   g4had_run ${JobID} ${nevt} &

done

#wait until all background jobs complete
wait

exit

