#!/usr/bin/env bash
#
# cleanup logs/err from earlier rounds
#
# ---> not now ---> /bin/rm -f slurm*.out

pl0=(ShieldingM Shielding ftfp_bert qgsp_bert)

for (( p=0; p<${#pl0[@]}; ++p )) do

sbatch -N 1 -n 16 -c 1 -p cpu_gce -A g4v --exclusive --qos=regular --time=23:00:00 \
 	--export="G4RELEASE=${1},beam=proton,target=Ta,momentum=8000,physlist=${pl0[$p]},EXP=HARP,nevents=62500" \
	slurm_master.sh sim_perCore_wcic.sh

done

tgt1=(C Be)
mom1=(3000 5000 8000 12000)
bm1=(proton piplus piminus)
pl1=(NuBeam ftfp_bert)

for (( i=0; i<${#tgt1[@]}; ++i )) do
for (( j=0; j<${#mom1[@]}; ++j )) do
for (( k=0; k<${#bm1[@]}; ++k )) do
for (( p=0; p<${#pl1[@]}; ++p )) do

sbatch -N 1 -n 16 -c 1 -p cpu_gce -A g4v --exclusive --qos=regular --time=23:00:00 \
 	--export="G4RELEASE=${1},beam=${bm1[$k]},target=${tgt1[$i]},momentum=${mom1[$j]},physlist=${pl1[$p]},EXP=HARP,nevents=62500" \
	slurm_master.sh sim_perCore_wcic.sh

done
done
done
done

pl2=(NuBeam ftfp_bert qgsp_bert)

for (( p=0; p<${#pl2[@]}; ++p )) do

sbatch -N 1 -n 16 -c 1 -p cpu_gce -A g4v --exclusive --qos=regular --time=23:00:00 \
 	--export="G4RELEASE=${1},beam=proton,target=C,momentum=31000,physlist=${pl2[$p]},EXP=NA61,nevents=62500" \
	slurm_master.sh sim_perCore_wcic.sh

sbatch -N 1 -n 16 -c 1 -p cpu_gce -A g4v --exclusive --qos=regular --time=23:00:00 \
 	--export="G4RELEASE=${1},beam=proton,target=C,momentum=158000,physlist=${pl2[$p]},EXP=NA49,nevents=62500" \
	slurm_master.sh sim_perCore_wcic.sh

done

