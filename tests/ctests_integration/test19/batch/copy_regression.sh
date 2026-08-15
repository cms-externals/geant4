#!/usr/bin/env bash
#

release=${1}
beam=${2}
dir=${3}

# move regression plots to results area
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release}/test19/regre ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test19/regre
fi
cp ${beam}*-regre*.gif /g4/g4p/pbs/g4-had-validation/results/${release}/test19/regre/.
#
if [ ! -d /g4/g4p/pbs/g4-had-validation/regression-test-files/test19/${release} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/regression-test-files/test19/${release}
fi
if [ ! -d /g4/g4p/pbs/g4-had-validation/regression-test-files/test19/${release}/${dir} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/regression-test-files/test19/${release}/${dir}
fi
cp ./${dir}/${beam}*.root /g4/g4p/pbs/g4-had-validation/regression-test-files/test19/${release}/${dir}/.

