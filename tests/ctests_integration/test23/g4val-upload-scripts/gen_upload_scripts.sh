#!/usr/bin/env bash
#
# usage:
# source generate_upload_scripts.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

# NA61
./gen_upload_script_na61_models.sh proton C ${g4version} ${status}
./gen_upload_script_na61_regre.sh proton C ${g4version} ${status}

#NA49
./gen_upload_script_na49_models.sh proton C ${g4version} ${status}
./gen_upload_script_na49_regre.sh proton C ${g4version} ${status}

#HARP
./gen_upload_script_harp_models.sh proton  C  ${g4version} ${status}
./gen_upload_script_harp_regre.sh  proton  C  ${g4version} ${status}
./gen_upload_script_harp_models.sh piplus  C  ${g4version} ${status}
./gen_upload_script_harp_regre.sh  piplus  C  ${g4version} ${status}
./gen_upload_script_harp_models.sh piminus C  ${g4version} ${status}
./gen_upload_script_harp_regre.sh  piminus C  ${g4version} ${status}
#
./gen_upload_script_harp_models.sh proton  Be ${g4version} ${status}
./gen_upload_script_harp_regre.sh  proton  Be ${g4version} ${status}
./gen_upload_script_harp_models.sh piplus  Be ${g4version} ${status}
./gen_upload_script_harp_regre.sh  piplus  Be ${g4version} ${status}
./gen_upload_script_harp_models.sh piminus Be ${g4version} ${status}
./gen_upload_script_harp_regre.sh  piminus Be ${g4version} ${status}
#
./gen_upload_script_harp_models.sh proton  Ta ${g4version} ${status}
./gen_upload_script_harp_regre.sh  proton  Ta ${g4version} ${status}
