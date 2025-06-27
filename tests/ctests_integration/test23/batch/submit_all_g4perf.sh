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

harp_list=`ls -l | grep HARP_SLURM | awk '{print $9}'`

echo " harp_list = ${harp_list} ... "

for test in ${harp_list} ; do

sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=3000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=5000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=8000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=12000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}

sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=3000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=5000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=8000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=12000,physlist=NuBeam,nevtotal=1000000" -A g4p ${test}

sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=3000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=5000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=8000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=12000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}

sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=3000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=5000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=8000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=Be,momentum=12000,physlist=ftfp_bert,nevtotal=1000000" -A g4p ${test}

done

#
# NOTE (JVY): drop the -t arg if more than 24hrs is needed (default=48)
#
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=Shielding,nevtotal=1000000" -A g4p proton_HARP_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ShieldingM,nevtotal=1000000" -A g4p proton_HARP_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=ftfp_bert,nevtotal=1000000" -A g4p proton_HARP_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=Ta,momentum=8000,physlist=qgsp_bert,nevtotal=1000000" -A g4p proton_HARP_SLURM.sh

sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=NuBeam,nevtotal=31250" -A g4p proton_NA61_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=ftfp_bert,nevtotal=31250" -A g4p proton_NA61_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 -t 24:00:00 --export="G4RELEASE=${g4version},target=C,momentum=31000,physlist=qgsp_bert,nevtotal=31250" -A g4p proton_NA61_SLURM.sh

#
# NOTE (JVY): drop the -t arg if more than 24hrs is needed (default=48)
#
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=NuBeam,nevtotal=31250" -A g4p proton_NA49_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=ftfp_bert,nevtotal=31250" -A g4p proton_NA49_SLURM.sh
sbatch --qos=g4perf -p amd32_g4perf -N1 --export="G4RELEASE=${g4version},target=C,momentum=158000,physlist=qgsp_bert,nevtotal=31250" -A g4p proton_NA49_SLURM.sh


