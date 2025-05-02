#!/usr/bin/env bash
#

release=${1}
beam=${2}

if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release} ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}
fi
if [ ! -d /home/g4p/pbs/g4-had-validation/results/${release}/test47 ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test47
fi
if [ ! -d /g4/g4p/pbs/g4-had-validation/results/${release}/test47/models ]; then
mkdir /g4/g4p/pbs/g4-had-validation/results/${release}/test47/models
fi
cp ${beam}*-models.gif /g4/g4p/pbs/g4-had-validation/results/${release}/test47/models/.
#
