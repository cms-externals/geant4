#!/usr/bin/env bash
#

# NOTE: python v2_7_10 & numpy v1_9_2 must be setup
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

# ---> BeamDetails=( '10.0GeV=Itep proton 10GeV (1)' '10.0GeV=Itep proton 10GeV (2)' )
Targets=( 'Al' 'Au' 'Be' 'Cu' 'Ta' )
# ---> Angles=( '0' '3.5' '10.5' '10.8' '59' '61' '90' '97' '119' )

ModelDetails=( 'ftfp=FTFP' )

gdir=/g4/g4p/pbs/g4-had-validation/regression-test-files
g4version=${1}
g4vtag=${2}
access=${3}

if [ "x" == "x${access}" ]; then
   echo "ACCESS is not specified; setting it to \"public\""
   access=public
fi

if [ "x" == "x${g4version}" ]; then
   echo "please, provide geant4 version as the 1st input argument - it is mandatory"
   exit
fi

upload_dir=${gdir}/test47/${g4version}/g4val-upload-json
if [ ! -d "${upload_dir}" ]; then
/bin/mkdir ${upload_dir}
fi

for (( i=0; i<${#ModelDetails[@]}; ++i )) do
model=`echo ${ModelDetails[$i]} | awk -F '=' '{print $1}'`
mid=`echo ${ModelDetails[$i]} | awk -F '=' '{print $2}'`

for (( k=0; k<${#Targets[@]}; ++k )) do
target=${Targets[$k]}

if [ -e PbarProd-proton${target}10.0GeV-metadata-${model}.json ]; then
/bin/rm PBarProd-proton${target}10.0GeV-metadata-${model}.json
fi

sed "s/VTAG/$g4vtag/; s/MID/$mid/; s/TGT/$target/; s/ACCESS/$access/" PbarProd-metadata-template.json > PbarProd-proton${target}10.0GeV-metadata-${model}.json

python ../../uploader/DoSSiERconverter.py -c convert -o PbarProd-proton${target}10.0GeV-${model}.json \
	--metadatafile PbarProd-proton${target}10.0GeV-metadata-${model}.json \
	${gdir}/test47/${g4version}/proton${target}10.0GeV${model}.root:XSecVsMomPbarAt10.5deg \
	${gdir}/test47/${g4version}/proton${target}10.0GeV${model}.root:XSecVsMomPbarAt10.8deg 

/bin/mv PbarProd-proton${target}10.0GeV-${model}.json ${upload_dir}/.

/bin/rm PbarProd-proton${target}10.0GeV-metadata-${model}.json

done

done

