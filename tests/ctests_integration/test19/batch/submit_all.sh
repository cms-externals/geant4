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

# ---> test_list=`ls -l | grep run | awk '{print $9}'`
test_list=`ls -l | grep harp | awk '{print $9}'`

target=(C Be Ta)
momentum=(3000 5000 8000 12000)

for test in ${test_list} ; do
for (( i=0; i<${#target[@]}; ++i )) do
for (( j=0; j<${#momentum[@]}; ++j )) do
   echo "... Processing ... ${test} ... for ${1} for target ${target[$i]} and momentum ${momentum[$j]}"
## ---> PBS/obsolete --->   qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2},target=${target[$i]},momentum=${momentum[$j]} -A g4p ${test}
#
# NOTE (JVY): most of the jobs require less than 24 hours of CPU but there're a couple (Ta target) than need more
#             maybe we should just drop the -t arg and make them all default to 48 ???
#
sbatch -p amd32_g4val_slow -t 24:00:00 --export="G4RELEASE=${g4version},G4REGRESSION=${2},target=${target[$i]},momentum=${momentum[$j]}" -A g4p ${test}
done
done
done

# special case
#
# ---> PBS/obsolete ---> qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2},target=Be,momentum=8900 -A g4p test19_run_harp_proton.sh
sbatch -p amd32_g4val_slow -t 24:00:00 --export="G4RELEASE=${g4version},G4REGRESSION=${2},target=Be,momentum=8900" -A g4p test19_run_harp_proton.sh


# NA61 and NA49
#
# ---> PBS/obsolete ---> qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2},target=C -A g4p test19_run_na61.sh
sbatch -p amd32_g4val_slow -t 24:00:00 --export="G4RELEASE=${g4version},G4REGRESSION=${2},target=C" -A g4p test19_run_na61.sh
#
# ---> PBS/obsolete ---> qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2},target=C -A g4p test19_run_na49.sh
#
# If -t argument is not specified, the jobs are supposed to default to 48 hours
#
sbatch -p amd32_g4val_slow --export="G4RELEASE=${g4version},G4REGRESSION=${2},target=C" -A g4p test19_run_na49.sh

test_sasm6e=`ls -l | grep SASM6E | awk '{print $9}'`
target_sasm6e=(C Cu Pb)
#
for test in ${test_sasm6e} ; do
for (( i=0; i<${#target_sasm6e[@]}; ++i )) do
   echo "... Processing ... ${test} ... for ${1} for target ${target_sasme6[$i]}"
# ---> PBS/obsolete --->   qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v App=${test},G4RELEASE=${g4version},target=${target_sasm6e[$i]},momentum=100000 -A g4p pbs_multiCore_master.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=${target_sasm6e[$i]},momentum=100000,nevents=31250" -A g4p ${test}
done
done
