#!/usr/bin/env bash
#
# Once submitted to batch queue via SLURM, this script will check 
# the number of cores per worker node (e.g. 32 for an amd32 node)
# and will launch equivalent number of parallel jobs, splitting 
# the total statistics evenly among those jobs
#
# Usage example:
# sbatch -p amd32_g4val_slow -N1 -t 24:00:00 \
#        --export="G4RELEASE=${g4version},target=C,momentum=3000,physlist=NuBeam,nevents=1000000" \
#        -A g4p test_proton_HARP.sh
# 

g4had_run() {

JobID=$1
nevents=$2

#target=$3
#momentum=$4
#physlist=$5
#G4RELEASE=geant4-10-05

source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh
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

# G4WORKDIR=${PWD}/..
echo " G4WORKDIR = ${G4WORKDIR}"

cd ${G4WORKDIR}
echo " We are in ${PWD}" 

seed=$((1234+${JobID}))
echo " random seed = ${seed}"

config=proton.${target}.${momentum}.HARP.${physlist}.${JobID}

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
/usr/bin/printf "#particle \nproton \n#isHARP \n#position(mm) \n0. 0. -100. \n#direction \n0. 0. 1. \n" >> ${config}

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

# to be passed in later...
#target=Ta
#momentum=8000
#physlist=ftfp_bert
#nevtotal=100000
nevt=$(( nevtotal/coresPerNode ))

for JobID in $(seq ${coresPerNode}); do

#   g4had_run ${JobID} ${nevt} ${target} ${momentum} ${physlist} >& proton.${target}.${momentum}.${physlist}-${JobID}.log &
#   g4had_run ${JobID} ${nevt} >& proton.${target}.${momentum}.${physlist}-${JobID}.log &
   g4had_run ${JobID} ${nevt} &

done

#wait until all background jobs complete
wait

exit

