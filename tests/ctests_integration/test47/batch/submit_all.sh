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

test_list=`ls -l | grep run | awk '{print $9}'`

for test in ${test_list} ; do
  echo "... Processing ... ${test} ... for ${1}"
# ---> PBS/obsolete --->  qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2} -A g4p ${test}
#
# NOTE-1 (JVY): it looks like the jobs need a bit longer than 24 hours but there's no guarantee that they'll fit into 24 sharp
# NOTE-2 (JVY): if more thatn 24 hrs of CPU is required, just drop the -t option, as amd32_g4val_slow will default to 48 hours;
#               however, with explicitly 48 hours or just anything larger than 24 hours, the job(s) will pend forever;
#               all in all, we're just dropping the -t here and go with defaults
#
sbatch -p amd32_g4val_slow --export="G4RELEASE=${g4version},G4REGRESSION=${2}" -A g4p ${test}
done

target=(Be Al Cu Au Ta)

for (( k=0; k<${#target[@]}; ++k )) do
  echo "... Processing ... proton+${target[$k]} --> pbar"
#  qsub -q amd32_g4val -l nodes=1:slow -v App=test47_pbarprod_perCore.sh,G4RELEASE=${g4version},target=${target[$k]},momentum=9956,model=ftfp  -A g4p pbs_multiCore_master.sh 
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=${target[$k]},momentum=9956,model=ftfp,nevents=312500" -A g4p test47_pbarprod_perCore.sh
done
 
