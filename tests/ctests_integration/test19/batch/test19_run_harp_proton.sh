#!/usr/bin/env bash
#

source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh

# setup G4 datasets
#
source /g4/g4p/pbs/g4-had-validation/env-setup/g4-datasets-setup-${G4RELEASE}.sh

export G4INSTALL=/g4/g4p/pbs/g4-had-validation/g4-releases/${G4RELEASE}

echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"

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

config=test19.harp.proton.${target}.${momentum}

if [ -e ${G4WORKDIR}/${config} ]; then
/bin/rm -f ${G4WORKDIR}/${config}
fi

mom=$((${momentum}/1000))

/usr/bin/printf "#verbose \n-1 \n#rad \n" >> ${config}
/usr/bin/printf "#events \n" >> ${config}
/usr/bin/printf "1000000 \n" >> ${config}

/usr/bin/printf "//--------Proton_processes \n" >> ${config}
/usr/bin/printf "#particle \nproton \n#isHARP \n#position(mm) \n0. 0. 0. \n#direction \n0. 0. 1. \n" >> ${config}

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
/usr/bin/printf "inclxx \n" >> ${config}
/usr/bin/printf "#run \n" >> ${config}

/usr/bin/printf "#exit\n" >> ${config}


./test19 ${config}

if [ -e ${config} ]; then
/bin/rm -f ${config}
fi

$ROOTSYS/bin/root -b -p -q HARPModels.C\(\"proton\",\"${target}\",${mom}\)

if [ "x" == "x$G4REGRESSION" ]; then
   echo "Regression testing is NOT requested"
else
   echo "G4REGRESSION says: $G4REGRESSION"
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"${target}\",\"${mom}\"\)
fi

if [ "x" == "x$G4RELEASE" ] ; then
    echo "Variable G4RELEASE is not set"
else
#
. batch/copy_results.sh ${G4RELEASE} proton
#
if [ "x" == "x$G4REGRESSION" ]; then
   echo "NO regression plots"
else
# move regression plots to results area
. batch/copy_regression.sh ${G4RELEASE} proton harp-histo
#
fi
#
fi
