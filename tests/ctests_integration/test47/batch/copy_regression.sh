#!/usr/bin/env bash
#

release=${1}
beam=${2}

#
# if [ "x" == "x$G4REGRESSION" ]; then
#   echo "NO regression plots"
# else
# move regression plots to results area
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release}/test47/regre ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test47/regre
fi
cp ${beam}*-regre*.gif /g4/g4p/pbs/g4-had-validation/results/${release}/test47/regre/.
#
if [ ! -d /g4/g4p/pbs/g4-had-validation/regression-test-files/test47/${release} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/regression-test-files/test47/${release}
fi
cp ${beam}*.root /g4/g4p/pbs/g4-had-validation/regression-test-files/test47/${release}/.
#
# fi
