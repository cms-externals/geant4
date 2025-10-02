#!/usr/bin/env bash
#

# EL8

sbatch -N 1 -n 12 -c 1 -p wc_cpu --exclusive --qos=regular --time=16:00:00 -A g4v \
	slurm_master.sh test19_run_harp_wcic.sh proton 

sbatch -N 1 -n 12 -c 1 -p wc_cpu --exclusive --qos=regular --time=16:00:00 -A g4v \
	slurm_master.sh test19_run_harp_wcic.sh piplus

sbatch -N 1 -n 12 -c 1 -p wc_cpu --exclusive --qos=regular --time=16:00:00 -A g4v \
	slurm_master.sh test19_run_harp_wcic.sh piminus

sbatch -N 1 -n 1 -c 1 -p wc_cpu --exclusive --qos=regular --time=16:00:00 -A g4v \
	slurm_master.sh test19_run_na61_wcic.sh proton

sbatch -N 1 -n 1 -c 1 -p wc_cpu --exclusive --qos=regular --time=16:00:00 -A g4v \
	slurm_master.sh test19_run_na61_wcic.sh piplus


# SL7

#sbatch -N 1 -n 12 -c 1 -p cpu_gce --exclusive --qos=regular --time=16:00:00 -A g4v \
# 	--export="G4RELEASE=${1},beam=proton" \
#	slurm_master.sh test19_run_harp_wcic.sh 
#
#sbatch -N 1 -n 12 -c 1 -p cpu_gce --exclusive --qos=regular --time=16:00:00 -A g4v \
# 	--export="G4RELEASE=${1},beam=piplus" \
#	slurm_master.sh test19_run_harp_wcic.sh 

#sbatch -N 1 -n 12 -c 1 -p cpu_gce --exclusive --qos=regular --time=16:00:00 -A g4v \
# 	--export="G4RELEASE=${1},beam=piminus" \
#	slurm_master.sh test19_run_harp_wcic.sh 

#sbatch -N 1 -n 1 -c 1 -p cpu_gce --exclusive --qos=regular --time=16:00:00 -A g4v \
# 	--export="G4RELEASE=${1}" \
#	slurm_master.sh test19_run_na61_wcic.sh 

#sbatch -N 1 -n 1 -c 1 -p cpu_gce --exclusive --qos=regular --time=16:00:00 -A g4v \
# 	--export="G4RELEASE=${1}" \
#	slurm_master.sh test19_run_na49_wcic.sh 

#target_sasm6e=(C Cu Pb)
#beam_sasm6e=(proton piplus)

#for (( i=0; i<${#target_sasm6e[@]}; ++i )) do
#for (( j=0; j<${#beam_sasm6e[@]}; ++j )) do

#sbatch -N 1 -n 16 -c 1 -p cpu_gce --exclusive --qos=regular --time=23:00:00 -A g4v \
# 	--export="G4RELEASE=${1},beam=${beam_sasm6e[$j]},target=${target_sasm6e[$i]},nevents=62500" \
#	slurm_master.sh test19_run_sasm6e_wcic.sh 

#done
#done
