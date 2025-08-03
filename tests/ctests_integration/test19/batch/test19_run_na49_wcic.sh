#!/usr/bin/env bash
#

source /cvmfs/geant4-ib.opensciencegrid.org/products/setup
setup root v6_24_06b -q e20:p399:prof

# --> export G4INSTALL=/work1/g4v/yarba_j/geant4-local-builds/gcc-8.2.0/geant4-${G4RELEASE}

# setup G4 datasets
#
source /wclustre/g4p/g4p/download/g4data/g4-datasets-${G4RELEASE}.sh

echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"

# setup workdir
#
if [ "x" == "x$G4WORKDIR" ] ; then
G4WORKDIR=${SLURM_SUBMIT_DIR}/..
else
    echo "Variable says: $G4WORKDIR"
    echo "Variable SLURM_SUBMIT_DIR says: $SLURM_SUBMIT_DIR"
fi

cd ${G4WORKDIR}

./test19 test19.na49



