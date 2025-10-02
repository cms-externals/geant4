#!/usr/bin/env bash
#

# 1st input arg (version) is mandatory;
# 2nd one (regression testing flag) is optional

g4version=${1}

if [ "x" == "x$g4version" ]; then
echo " you must specify Geant4 release/version; exit"
exit 3
fi

source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh

# ---> PBS/obsolete --->
#if [ `echo ${PATH} | grep pbs` ]; then
#echo "PBS is already set"
#else
#export PATH=${PATH}:/usr/local/pbs/bin
#fi

#
# cleanup logs/err from earlier rounds
#
/bin/rm -f *.sh.o*
/bin/rm -f *.sh.e*
/bin/rm -f slurm*.out

test_list=`ls -l | grep Mu2e | awk '{print $9}'`

for test in ${test_list} ; do
  echo "... Processing ... ${test} ... for ${1}"
# ---> PBS/obsolete --->  qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2} -A g4p ${test}
#
# NOTE (JVY): if one wants more than 24hrs, one can just drop the -t arg, and the CPU will default to 48hrs
#             however, trying to explicitly request anything larger than 24hrs will make the job pend forever 
#
sbatch -p amd32_g4val_slow --export="G4RELEASE=${g4version},G4REGRESSION=${2}" -A g4p ${test}
done
