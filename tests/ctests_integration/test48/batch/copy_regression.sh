#!/usr/bin/env bash
#

release=${1}

#
# if [ "x" == "x$G4REGRESSION" ]; then
#   echo "NO regression plots"
# else
# move regression plots to results area
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release}/test48/regre ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test48/regre
fi
cp *regre*.gif /g4/g4p/pbs/g4-had-validation/results/${release}/test48/regre/.
#
if [ ! -d /g4/g4p/pbs/g4-had-validation/regression-test-files/test48/${release} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/regression-test-files/test48/${release}
fi
cp *.root /g4/g4p/pbs/g4-had-validation/regression-test-files/test48/${release}/.
#
# fi
