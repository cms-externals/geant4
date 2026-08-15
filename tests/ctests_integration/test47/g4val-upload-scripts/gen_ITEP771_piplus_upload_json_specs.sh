#!/usr/bin/env bash
#

# NOTE: python v2_7_10 & numpy v1_9_2 must be setup
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

#MCDetails=( 'bertini=148' 'inclxx=153' 'binary=149' 'ftfp=150'  )
#BeamDetails=( '1.40GeV=62' '5.00GeV=72' )
BeamDetails=( '1.40GeV=ITEP-771 (pi+,1.4GeV)' '5.00GeV=ITEP-771 (pi+,5GeV)' )
Targets=( 'C' 'Cu' 'Pb' 'U' )

ModelDetails=( 'bertini=Bertini' 'ftfp=FTFP' 'inclxx=INCLXX' 'binary=BIC' )

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

#for (( i=0; i<${#MCDetails[@]}; ++i )) do
#model=`echo ${MCDetails[$i]} | awk -F '=' '{print $1}'`
#mid=`echo ${MCDetails[$i]} | awk -F '=' '{print $2}'`

for (( i=0; i<${#ModelDetails[@]}; ++i )) do
model=`echo ${ModelDetails[$i]} | awk -F '=' '{print $1}'`
mid=`echo ${ModelDetails[$i]} | awk -F '=' '{print $2}'`

for (( j=0; j<${#BeamDetails[@]}; ++j )) do
bmom=`echo ${BeamDetails[$j]} | awk -F '=' '{print $1}'`
bid=`echo ${BeamDetails[$j]} | awk -F '=' '{print $2}'`

for (( k=0; k<${#Targets[@]}; ++k )) do
target=${Targets[$k]}

if [ -e ITEP771-piplus${target}${bmom}-metadata-${model}.json ]; then
/bin/rm ITEP771-piplus${target}${bmom}-metadata-${model}.json
fi

sed "s/_MODEL/$model/; s/VTAG/$g4vtag/; s/MID/$mid/; s/_BPART/piplus/; s/_BMOM/$bmom/; s/BEAM/$bid/; s/_TGT/$target/; s/ACCESS/$access/" ITEP771-metadata-template.json > ITEP771-piplus${target}${bmom}-metadata-${model}.json

python ../../uploader/DoSSiERconverter.py -c convert -o ITEP771-piplus${target}${bmom}-${model}.json \
	--metadatafile ITEP771-piplus${target}${bmom}-metadata-${model}.json \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEproton0piplus${target}${model}${bmom}\ 59.1 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEproton0piplus${target}${model}${bmom}\ 89.0 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEproton0piplus${target}${model}${bmom}119.0 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEproton0piplus${target}${model}${bmom}159.6 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEneutron0piplus${target}${model}${bmom}\ 59.1 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEneutron0piplus${target}${model}${bmom}\ 89.0 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEneutron0piplus${target}${model}${bmom}119.0 \
	${gdir}/test47/${g4version}/piplus${target}${model}${bmom}.root:KEneutron0piplus${target}${model}${bmom}159.6 

/bin/mv ITEP771-piplus${target}${bmom}-${model}.json ${upload_dir}/.

/bin/rm ITEP771-piplus${target}${bmom}-metadata-${model}.json

done

done

done


