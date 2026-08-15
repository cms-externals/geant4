#!/usr/bin/env bash
#

release=${1}

if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}
fi
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${release}/test48 ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test48
fi
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release}/test48/models ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test48/models
fi
cp *models.gif /g4/g4p/pbs/g4-had-validation/results/${release}/test48/models/.
#
