#!/usr/bin/env bash
#

source /products/setup

if [ `echo ${PATH} | grep pbs` ]; then
echo "PBS is already set"
else
export PATH=${PATH}:/usr/local/pbs/bin
fi

#G4RELEASE=${1}
#echo "G4RELEASE = ${G4RELEASE}"
#test_list=`ls -l | grep run | awk '{print $9}'`
#for test in ${test_list} ; do
#  echo "... Processing ... ${test} ... for ${1}"
#   qsub -q grunt -v G4RELEASE=${1} ${test}
#done


harp_list=`ls -l | grep HARP | awk '{print $9}'`

for test in ${harp_list} ; do

qsub -q perf -v  G4RELEASE=${1},target=C,momentum=3000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=C,momentum=5000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=C,momentum=8000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=C,momentum=12000,physlist=NuBeam,nevents=1000000  ${test}

qsub -q perf -v  G4RELEASE=${1},target=Be,momentum=3000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=Be,momentum=5000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=Be,momentum=8000,physlist=NuBeam,nevents=1000000 ${test}
qsub -q perf -v  G4RELEASE=${1},target=Be,momentum=12000,physlist=NuBeam,nevents=1000000 ${test}

done

qsub -q perf -v  G4RELEASE=${1},target=Be,momentum=8900,physlist=NuBeam,nevents=1000000 proton_HARP.sh

qsub -q perf -v G4RELEASE=${1},target=C,momentum=31000,physlist=NuBeam,nevents=1000000 proton_NA61.sh

qsub -q perf -v G4RELEASE=${1},target=C,momentum=158000,physlist=NuBeam,nevents=1000000 proton_NA49.sh



