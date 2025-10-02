#!/usr/bin/env bash
#

G4RELEASE=${1}

if [ "x" == "x$G4RELEASE" ]; then
echo " you must specify the Geang4 release/version; exit"
exit 3
fi

if [ "x" == "x$ROOTSYS" ]; then
#source /products/setup
#setup root v5_23_01 -q "e2:prof"
source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh
fi

echo " *** !!! Before you run this, you MUST make sure that ALL your simulation jobs have finished"
echo " *** !!! You MUST also ensure that you have revisited test23/shared-root-macros/REGRESSION_TEST.h"
echo " *** !!! and have sspecified Geant4 versions you wish to compare"

$ROOTSYS/bin/root -b -p -q GammaNRegre.C\(\)

. batch/copy_regression.sh ${G4RELEASE}

