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

$ROOTSYS/bin/root -b -p -q PiMinusRegre.C
$ROOTSYS/bin/root -b -p -q KMinusRegre.C
$ROOTSYS/bin/root -b -p -q AntiProtonRegre.C

$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Al\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Si\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Ca\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Fe\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Ag\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"I\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Au\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"Pb\"\)
$ROOTSYS/bin/root -b -p -q MuMinusRegre.C\(\"S\"\)

. batch/copy_regression.sh ${G4RELEASE}
