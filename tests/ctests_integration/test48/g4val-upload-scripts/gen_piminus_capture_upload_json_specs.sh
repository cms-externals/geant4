#!/usr/bin/env bash
#

# NOTE: python v2_7_10 & numpy v1_9_2 must be setup
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

MCDetails=( 'bertini=148' )
TargetDetails=( 'C=6' 'N=7' 'O=8' 'Al=13' 'Cu=29' 'Ta=73' 'Pb=82' )

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

for (( k=0; k<${#TargetDetails[@]}; ++k )) do
target=`echo ${TargetDetails[$k]} | awk -F '=' '{print $1}'`
tgtid=`echo ${TargetDetails[$k]} | awk -F '=' '{print $2}'`

if [ -e piminus-capture-${target}-metadata-${model}.json ]; then
/bin/rm piminus-capture-${target}-metadata-${model}.json
fi

sed "s/VTAG/$g4tag/; s/TGT/$target/; s/ACCESS/$access/" piminus-capture-metadata-template.json > piminus-capture-${target}-metadata-${model}.json

python ../../uploader/DoSSiERconverter.py -c convert -o piminus-capture-${target}-${model}.json \
	--metadatafile piminus-capture-${target}-metadata-${model}.json \
	${gdir}/test48/${g4version}/piminus${target}BertiniPreCo.root:NvsT 

/bin/mv piminus-capture-${target}-${model}.json ${upload_dir}/.

/bin/rm piminus-capture-${target}-metadata-${model}.json

done

done
