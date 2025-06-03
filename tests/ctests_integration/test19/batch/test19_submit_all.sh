#!/usr/bin/env bash
#

# 1st input arg (version) is mandatory;
# 2nd one (regression testing flag) is optional

g4version=${1}

if [ "x" == "x$g4version" ]; then
echo " you must specify Geant4 release/version; exit"
exit 3
fi

#source /products/setup
source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh


if [ `echo ${PATH} | grep pbs` ]; then
echo "PBS is already set"
else
export PATH=${PATH}:/usr/local/pbs/bin
fi

#
# cleanup logs/err from earlier rounds
#
/bin/rm -f *.sh.o*
/bin/rm -f *.sh.e*

test_list=`ls -l | grep run | awk '{print $9}'`

for test in ${test_list} ; do
   echo "... Processing ... ${test} ... for ${1}"
#   qsub -q grunt -v G4RELEASE=${g4version},G4REGRESSION=${2} ${test}
#   qsub -q amd32_Geant4 -l nodes=1:amd32,walltime=48:00:00 -v G4RELEASE=$G4RELEASE -A g4p test19_run_harp_proton.sh 
   qsub -q amd32_Geant4 -l nodes=1:amd32,walltime=48:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2} -A g4p ${test}
done
