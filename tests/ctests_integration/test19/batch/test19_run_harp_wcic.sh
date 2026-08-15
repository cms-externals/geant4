#!/usr/bin/env bash
#

# SL7
# source /cvmfs/geant4-ib.opensciencegrid.org/products/setup
# setup root v6_24_06b -q e20:p399:prof

# --> export G4INSTALL=/work1/g4v/yarba_j/geant4-local-builds/gcc-9.3.0/geant4-${G4RELEASE}

# setup G4 datasets
#
# --> no need --> source /wclustre/g4p/g4p/download/g4data/g4-datasets-${G4RELEASE}.sh

echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${SLURM_SUBMIT_DIR}/..
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"
fi

# --->

JobID=$((1+${SLURM_PROCID}))

target_list=( Be C Ta )
mom_list=( 3000 5000 8000 12000 )

ntgts=${#target_list[@]}
nmoms=${#mom_list[@]}

njobsmax=$((${ntgts}*${nmoms}))

# Note: see earlier that JobID = SLURM_PROCID + 1
#
if [ ${JobID} -gt ${njobsmax} ]; then
exit
fi

beam=${1}
itgt=$((${SLURM_PROCID}/${nmoms}))
imom=$((${SLURM_PROCID} % ${nmoms}))
target=${target_list[$itgt]}
momentum=${mom_list[$imom]}

# --->

cd ${G4WORKDIR}

config=test19.harp.${beam}.${target}.${momentum}

if [ -e ${G4WORKDIR}/${config} ]; then
/bin/rm -f ${G4WORKDIR}/${config}
fi

bm=${beam}
if [ "x${beam}" == "xpiplus" ]; then
   bm="pi+"
fi
if [ "x${beam}" == "xpiminus" ]; then
   bm="pi-"
fi

mom=$((${momentum}/1000))

/usr/bin/printf "#verbose \n-1 \n#rad \n" >> ${config}
/usr/bin/printf "#events \n" >> ${config}
/usr/bin/printf "1000000 \n" >> ${config}

/usr/bin/printf "//--------%s_processes \n" ${beam} >> ${config}
/usr/bin/printf "#particle \n" >> ${config}
/usr/bin/printf "%s\n" ${bm} >> ${config}
/usr/bin/printf "#isHARP \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n" >> ${config}

/usr/bin/printf "//-------- \n" >> ${config}

/usr/bin/printf "#material \n" >> ${config}
/usr/bin/printf "%s\n" ${target} >> ${config}

/usr/bin/printf "#momentum(MeV/c) \n" >> ${config}
/usr/bin/printf "%d\n" ${momentum} >> ${config}

/usr/bin/printf "// --- \n#generator \n" >> ${config}
/usr/bin/printf "bertini \n" >> ${config}
/usr/bin/printf "#run \n" >> ${config}
/usr/bin/printf "// --- \n#generator \n" >> ${config}
/usr/bin/printf "ftfp \n" >> ${config}
/usr/bin/printf "#run \n" >> ${config}
/usr/bin/printf "// --- \n#generator \n" >> ${config}
/usr/bin/printf "ftfp_tune3 \n" >> ${config}
/usr/bin/printf "#run \n" >> ${config}
#
# TMP deactivate validation of inclxx because of a random crashes that
#     happen for 8-12 GeV/c proton beam incident on light target (Be, C)
#
#/usr/bin/printf "// --- \n#generator \n" >> ${config}
#/usr/bin/printf "inclxx \n" >> ${config}
#/usr/bin/printf "#run \n" >> ${config}

/usr/bin/printf "#exit\n" >> ${config}

# EL8

module load gcc/11.4.0

export HOME=/work1/g4p/g4p/products-el8
export SPACK_ROOT=/work1/g4p/g4p/products-el8/spack
. $SPACK_ROOT/share/spack/setup-env.sh
spack load root@6.26.06

./test19 ${config}

#if [ -e ${config} ]; then
#/bin/rm -f ${config}
#fi


