#!/usr/bin/env bash
#

# in principle, this is not needed because the macro is typically
# executed interactively which means that Root is already setup
#
if [ "x" == "x$ROOTSYS" ] ; then
source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh
fi

G4RELEASE=${1}
G4REGRESSION=${2}

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",8000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"C\",12000\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",8000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"C\",12000\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",8000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"C\",12000\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",8000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piplus\",\"Be\",12000\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",8000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"piminus\",\"Be\",12000\)

$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",3000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",5000\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",8000\)
# ---> $ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",8900\)
$ROOTSYS/bin/root -b -p -q HARPRegre.C\(\"proton\",\"Be\",12000\)

$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"proton\",\"C\"\)
$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"piplus\",\"C\"\)
$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"piminus\",\"C\"\)

$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"proton\",\"Be\"\)
$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"piplus\",\"Be\"\)
$ROOTSYS/bin/root -b -p -q HARPModelsAllEnergies.C\(\"piminus\",\"Be\"\)

$ROOTSYS/bin/root -b -p -q NA61Models.C
$ROOTSYS/bin/root -b -p -q NA61Regre.C\(\"piminus\"\)
$ROOTSYS/bin/root -b -p -q NA61Regre.C\(\"piplus\"\)
$ROOTSYS/bin/root -b -p -q NA61Regre.C\(\"proton\"\)

$ROOTSYS/bin/root -b -p -q NA49Models.C
$ROOTSYS/bin/root -b -p -q NA49Regre.C
 
if [ "x" == "x$G4RELEASE" ] ; then
    echo "Variable G4RELEASE is not set"
else
#
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}
fi
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23 ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23
fi
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/models ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/models
fi
cp *-models*.gif /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/models/.
#
if [ "x" == "x$G4REGRESSION" ]; then
   echo "NO regression plots"
else
# move regression plots to results area
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/regre ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/regre
fi
cp *-regre*.gif /g4/g4p/pbs/g4-had-validation/results/${G4RELEASE}/test23/regre/.
fi
#
fi
