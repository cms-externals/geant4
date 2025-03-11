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

# ---> OBSOLETE
#if [ `echo ${PATH} | grep pbs` ]; then
#echo "PBS is already set"
#else
#export PATH=${PATH}:/usr/local/pbs/bin
#fi

harp_list=`ls -l | grep HARP_perCore | awk '{print $9}'`

echo " harp_list = ${harp_list} ... "

for test in ${harp_list} ; do

# ---> PBS/obsolete --->
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=3000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=5000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=8000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=12000,physlist=NuBeam -A g4p pbs_multiCore_master.sh
#
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=3000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=5000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=8000,physlist=NuBeam  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=12000,physlist=NuBeam -A g4p pbs_multiCore_master.sh
#
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=3000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=5000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=8000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=C,momentum=12000,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh
#
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=3000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=5000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=8000,physlist=ftfp_bert  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=${test},G4RELEASE=${g4version},target=Be,momentum=12000,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh

sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=3000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=5000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=8000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=12000,physlist=NuBeam,nevents=31250" -A g4p ${test}

sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=3000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=5000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=8000,physlist=NuBeam,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=12000,physlist=NuBeam,nevents=31250" -A g4p ${test}

sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=3000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=5000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=8000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=12000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}

sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=3000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=5000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=8000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=12000,physlist=ftfp_bert,nevents=31250" -A g4p ${test}

done

#
# Hold off for now - apparently, there's some issue here as most of the jobs get killed (memory use limit ?),
# and those remaining run out of 48 hours of the requested wall time 
#
# Nov.9, 2015 - try to resume these tests since with 4.10.1.ref10 they went fine
#
# Apr.13, 2016 - after some more hiccups in 10.2 cycle (in radioactive decays...), Shielding(M) seem to be back on track; resume tests as of 4.10.2.ref03
#
# ---> PBS/obsolete --->
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=Shielding  -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ShieldingM -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=qgsp_bert -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh

#
# NOTE (JVY): drop the -t arg if more than 24hrs is needed (default=48)
#
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=Shielding,nevents=31250" -A g4p proton_HARP_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ShieldingM,nevents=31250" -A g4p proton_HARP_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ftfp_bert,nevents=31250" -A g4p proton_HARP_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=qgsp_bert,nevents=31250" -A g4p proton_HARP_perCore.sh

# ---> PBS/obsolete --->
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Be,momentum=8900,physlist=NuBeam    -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v  App=proton_HARP_perCore.sh,G4RELEASE=${g4version},target=Be,momentum=8900,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh

# ---> PBS/obsolete --->
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v App=proton_NA61_perCore.sh,G4RELEASE=${g4version},target=C,momentum=31000,physlist=NuBeam    -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v App=proton_NA61_perCore.sh,G4RELEASE=${g4version},target=C,momentum=31000,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=24:00:00 -v App=proton_NA61_perCore.sh,G4RELEASE=${g4version},target=C,momentum=31000,physlist=qgsp_bert -A g4p pbs_multiCore_master.sh

sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=NuBeam,nevents=31250" -A g4p proton_NA61_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=ftfp_bert,nevents=31250" -A g4p proton_NA61_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=qgsp_bert,nevents=31250" -A g4p proton_NA61_perCore.sh

# ---> PBS/obsolete --->
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_NA49_perCore.sh,G4RELEASE=${g4version},target=C,momentum=158000,physlist=NuBeam    -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_NA49_perCore.sh,G4RELEASE=${g4version},target=C,momentum=158000,physlist=ftfp_bert -A g4p pbs_multiCore_master.sh
#qsub -q amd32_g4val -l nodes=1:slow,walltime=48:00:00 -v App=proton_NA49_perCore.sh,G4RELEASE=${g4version},target=C,momentum=158000,physlist=qgsp_bert -A g4p pbs_multiCore_master.sh

#
# NOTE (JVY): drop the -t arg if more than 24hrs is needed (default=48)
#
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=NuBeam,nevents=31250" -A g4p proton_NA49_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=ftfp_bert,nevents=31250" -A g4p proton_NA49_perCore.sh
sbatch -p amd32_g4val_slow -N1 --array=0-31 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=qgsp_bert,nevents=31250" -A g4p proton_NA49_perCore.sh


