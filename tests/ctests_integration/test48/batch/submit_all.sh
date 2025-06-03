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

# ---> PBS/obsolete ---> qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v G4RELEASE=${g4version},G4REGRESSION=${2} -A g4p test48_run_piminus.sh
#
# NOTE-1 (JVY): 24 hours seem to be sufficient, no need for 48 hours
# NOTE-2 (JVY): if more CPU is required, just drop the -t option, as amd32_g4val_slow will default to 48 hours;
#               however, explicitly specifying 48 hours or anything larger than 24 hours will make the job(s) pend forever
#
sbatch -p amd32_g4val_slow -t 24:00:00 --export="G4RELEASE=${g4version},G4REGRESSION=${2}" -A g4p test48_run_piminus.sh

target=(Al Si Ca Fe Ag I Au Pb S)
for (( k=0; k<${#target[@]}; ++k )) do
# ---> PBS/obsolete --->   qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v G4RELEASE=${g4version},target=${target[$k]},G4REGRESSION=${2} -A g4p test48_run_muminus.sh
sbatch -p amd32_g4val_slow -t 24:00:00 --export="G4RELEASE=${g4version},target=${target[$k]},G4REGRESSION=${2}" -A g4p test48_run_muminus.sh
done

