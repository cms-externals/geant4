#!/usr/bin/env bash
#

# 1st input arg (version) is mandatory;
# 2nd one (regression testing flag) is optional

g4version=${1}

if [ "x" == "x$g4version" ]; then
echo " you must specify Geant4 release/version; exit"
exit 3
fi

source /products/setup

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
   qsub -q grunt -v G4RELEASE=${g4version},G4REGRESSION=${2} ${test}
done
