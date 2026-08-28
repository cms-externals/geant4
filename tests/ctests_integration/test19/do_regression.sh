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


$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",\"3.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",\"3.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Ta\",\"3.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",\"5.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",\"5.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Ta\",\"5.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",\"8.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",\"8.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Ta\",\"8.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",\"12.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",\"12.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Ta\",\"12.0\"\)

. batch/copy_regression.sh ${G4RELEASE} piminus harp-histo

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",\"3.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",\"3.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Ta\",\"3.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",\"5.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",\"5.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Ta\",\"5.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",\"8.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",\"8.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Ta\",\"8.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",\"12.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",\"12.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Ta\",\"12.0\"\)

. batch/copy_regression.sh ${G4RELEASE} piplus harp-histo

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",\"3.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",\"3.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Ta\",\"3.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",\"5.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",\"5.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Ta\",\"5.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",\"8.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",\"8.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Ta\",\"8.0\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",\"8.9\"\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",\"12.0\"\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",\"12.0\"\)
# $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Ta\",\"12.0\"\)

. batch/copy_regression.sh ${G4RELEASE} proton harp-histo

$ROOTSYS/bin/root -b -p -q NA61Regre.C\(\)

. batch/copy_regression.sh ${G4RELEASE} proton na61-histo


$ROOTSYS/bin/root -b -p -q NA49Regre.C\(\)

. batch/copy_regression.sh ${G4RELEASE} proton na49-histo

#if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19 ]; then
#mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19
#fi
#if [ ! -d /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19/regre ]; then
#mkdir /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19/regre
#fi
#cp piminus*-regre*.gif /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19/regre/.
#cp piplus*-regre*.gif /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19/regre/.
#cp proton*-regre*.gif /home/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test19/regre/.


