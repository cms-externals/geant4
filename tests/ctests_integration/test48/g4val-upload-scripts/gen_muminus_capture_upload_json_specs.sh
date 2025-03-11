#!/usr/bin/env bash
#

# NOTE: python v2_7_10 & numpy v1_9_2 must be setup
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

MCDetails=( 'bertini=148' )
TargetSingerDetails=( 'Al=13' 'Si=14' 'Ca=20' 'Fe=26' 'Ag=47' 'I=53' 'Au=79' 'Pb=82')
TargetSundelinDetails=( 'Si=14' 'S=16' 'Ca=20' )

gdir=/g4/g4p/pbs/g4-had-validation/regression-test-files
g4version=${1}
g4tag=${2}
access=${3}

if [ "x" == "x${access}" ]; then
   echo "ACCESS is not specified; setting it to \"public\""
   access=public
fi

if [ "x" == "x${g4version}" ]; then
   echo "please, provide geant4 version as the 1st input argument - it is mandatory"
   exit
fi

upload_dir=${gdir}/test48/${g4version}/g4val-upload-json
if [ ! -d "${upload_dir}" ]; then
/bin/mkdir ${upload_dir}
fi

for (( i=0; i<${#MCDetails[@]}; ++i )) do
model=`echo ${MCDetails[$i]} | awk -F '=' '{print $1}'`
mid=`echo ${MCDetails[$i]} | awk -F '=' '{print $2}'`

for (( k=0; k<${#TargetSingerDetails[@]}; ++k )) do
target=`echo ${TargetSingerDetails[$k]} | awk -F '=' '{print $1}'`
tgtid=`echo ${TargetSingerDetails[$k]} | awk -F '=' '{print $2}'`

if [ -e muminus-capture-Singer-${target}-metadata-${model}.json ]; then
/bin/rm muminus-capture-Singer-${target}-metadata-${model}.json
fi

sed "s/VTAG/$g4tag/; s/TGT/$target/; s/ACCESS/$access/" muminus-capture-metadata-template.json > muminus-capture-Singer-${target}-metadata-${model}.json

python ../../uploader/DoSSiERconverter.py -c convert -o muminus-capture-Singer-${target}-${model}.json \
	--metadatafile muminus-capture-Singer-${target}-metadata-${model}.json \
	${gdir}/test48/${g4version}/muminus${target}captureUpdate.root:NNeutrons 

/bin/mv muminus-capture-Singer-${target}-${model}.json ${upload_dir}/.

/bin/rm muminus-capture-Singer-${target}-metadata-${model}.json

done

for (( k=0; k<${#TargetSundelinDetails[@]}; ++k )) do
target=`echo ${TargetSundelinDetails[$k]} | awk -F '=' '{print $1}'`
tgtid=`echo ${TargetSundelinDetails[$k]} | awk -F '=' '{print $2}'`

if [ -e muminus-capture-Sundelin-${target}-metadata-${model}.json ]; then
/bin/rm muminus-capture-Sundelin-${target}-metadata-${model}.json
fi

sed "s/VTAG/$g4tag/; s/TGT/$target/; s/ACCESS/$access/" muminus-capture-metadata-template.json > muminus-capture-Sundelin-${target}-metadata-${model}.json

python ../../uploader/DoSSiERconverter.py -c convert -o muminus-capture-Sundelin-${target}-${model}.json \
	--metadatafile muminus-capture-Sundelin-${target}-metadata-${model}.json \
	${gdir}/test48/${g4version}/muminus${target}captureUpdate.root:NeutronKineticEnergy

/bin/mv muminus-capture-Sundelin-${target}-${model}.json ${upload_dir}/.

/bin/rm muminus-capture-Sundelin-${target}-metadata-${model}.json

done

done
