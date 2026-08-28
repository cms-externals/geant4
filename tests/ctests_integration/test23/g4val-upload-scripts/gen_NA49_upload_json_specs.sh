#!/usr/bin/env bash
#

# NOTE: python, numpy and ujson must be setup, all compatible with root
#       path to $ROOTSYS/lib must be added to PYTHONPATH
#       python uploader (by A.Dotti) must be installed:
#      git clone https://yarba@gitlab.cern.ch/PhysicsValidationDB/uploader.git

# Phys.lists are FTFP_BERT, QGSP_BERT & NuBeam
#
# MCDetails=( 'ftfp_bert=166' 'qgsp_bert=176' 'NuBeam=196' )
ModelDetails=( 'ftfp_bert=FTFP_BERT' 'qgsp_bert=QGSP_BERT' 'NuBeam=NuBeam'  )

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

if [ -e NA49-proton-C-metadata-integrated-spectra-${model}.json ]; then
/bin/rm NA49-proton-C-metadata-integrated-spectra-${model}.json
fi

# sed "s/XXX/$mid/" NA49-proton-C-metadata-integrated-spectra.json > NA49-proton-C-metadata-integrated-spectra-${model}.json
sed "s/MODEL/$mid/; s/VTAG/$g4vtag/" NA49-proton-C-metadata-integrated-spectra.json > NA49-proton-C-metadata-integrated-spectra-${model}.json

# python ../../uploader/plot_histofiles.py -c convert -o NA49-proton-C-integrated-spectra-${model}.json \
python ../../uploader/DoSSiERconverter.py -c convert -o NA49-proton-C-integrated-spectra-${model}.json \
	--metadatafile NA49-proton-C-metadata-integrated-spectra-${model}.json \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:piplus_dNdxF \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:piplus_pT \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:piminus_dNdxF \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:piminus_pT \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:proton_dNdxF  \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:proton_pT  \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:antiproton_dNdxF  \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:antiproton_pT  \
	${gdir}/test23/${g4version}/na49-histo/protonC158.0GeV${model}.root:neutron_dNdxF  

/bin/mv NA49-proton-C-integrated-spectra-${model}.json ${upload_dir}/.

/bin/rm NA49-proton-C-metadata-integrated-spectra-${model}.json

done
