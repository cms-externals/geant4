#!/usr/bin/env bash
#

# NOTE: python v2_7_10 & numpy v1_9_2 must be setup
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

# Phys.lists are FTFP_BERT, QGSP_BERT & NuBeam
#
# MCDetails=( 'ftfp_bert=166' 'qgsp_bert=176' 'NuBeam=196' )
ModelDetails=( 'ftfp_bert=FTFP_BERT' 'qgsp_bert=QGSP_BERT' 'NuBeam=NuBeam' )

gdir=/g4/g4p/pbs/g4-had-validation/regression-test-files
g4version=${1}
g4vtag=${2}

if [ "x" == "x${g4version}" ]; then
   echo "please, provide geant4 version - it is mandatory"
   exit
fi

upload_dir=${gdir}/test23/${g4version}/g4val-upload-json
if [ ! -d "${upload_dir}" ]; then
/bin/mkdir ${upload_dir}
fi

#for ((i=0; i<${#MCDetails[@]}; ++i )) do
#model=`echo ${MCDetails[$i]} | awk -F '=' '{print $1}'`
#mid=`echo ${MCDetails[$i]} | awk -F '=' '{print $2}'`

for (( i=0; i<${#ModelDetails[@]}; ++i )) do
model=`echo ${ModelDetails[$i]} | awk -F '=' '{print $1}'`
mid=`echo ${ModelDetails[$i]} | awk -F '=' '{print $2}'`

if [ -e NA61-proton-C-metadata-pions-${model}.json ]; then
/bin/rm NA61-proton-C-metadata-pions-${model}.json
fi

# sed "s/XXX/$mid/" NA61-proton-C-metadata-pions.json > NA61-proton-C-metadata-pions-${model}.json
sed "s/MODEL/$mid/; s/VTAG/$g4vtag/" NA61-proton-C-metadata-pions.json > NA61-proton-C-metadata-pions-${model}.json

# python ../../uploader/plot_histofiles.py -c convert -o NA61-proton-C-pions-${model}.json \
python ../../uploader/DoSSiERconverter.py -c convert -o NA61-proton-C-pions-${model}.json \
	--metadatafile NA61-proton-C-metadata-pions-${model}.json \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_0_10 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_10_20 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_20_40 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_40_60 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_60_100 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_100_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_140_180 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_180_240 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_240_300 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_300_360 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piplus_360_420 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_0_10 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_10_20 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_20_40 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_40_60 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_60_100 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_100_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_140_180 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_180_240 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_240_300 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_300_360 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:piminus_360_420

/bin/mv NA61-proton-C-pions-${model}.json ${upload_dir}/.

/bin/rm NA61-proton-C-metadata-pions-${model}.json

# sed "s/XXX/$mid/" NA61-proton-C-metadata-kaons.json > NA61-proton-C-metadata-kaons-${model}.json
sed "s/MODEL/$mid/; s/VTAG/$g4vtag/" NA61-proton-C-metadata-kaons.json > NA61-proton-C-metadata-kaons-${model}.json

#python ../../uploader/plot_histofiles.py -c convert -o NA61-proton-C-kaons-${model}.json \
python ../../uploader/DoSSiERconverter.py -c convert -o NA61-proton-C-kaons-${model}.json \
	--metadatafile NA61-proton-C-metadata-kaons-${model}.json \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_0_20 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_20_40 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_40_60 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_60_100 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_100_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_140_180 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_180_240 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus_240_300 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_0_20 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_20_40 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_40_60 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_60_100 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_100_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_140_180 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_180_240 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kminus_240_300 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_0_40 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_40_60 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_60_100 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_100_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_140_180 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:k0s_180_240 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus2piplus_20_140 \
	${gdir}/test23/${g4version}/na61-histo/protonC31.0GeV${model}.root:kplus2piplus_140_240 

/bin/mv NA61-proton-C-kaons-${model}.json ${upload_dir}/.

/bin/rm NA61-proton-C-metadata-kaons-${model}.json

done
