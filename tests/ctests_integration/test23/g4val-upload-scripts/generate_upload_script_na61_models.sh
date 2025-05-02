#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_na61_models.sh proton C geant4-10.00-ref08 internal
#

beam=${1}
target=${2}
g4version=${3}
status=${4}

# g4version=geant4-10.00-ref07
# status=internal

testname=test23

momentum=(31)
secondary=(piplus piminus proton)
region=(p1 p2)
# models=(ftfp_bert NuBeam NuBeam-3-3.5 qgsp_bert)

script=${testname}-na61-${beam}-${target}-models.xml

tg=${target}
if [ ${target} = "C" ]; then
tg=Carbon
fi
if [ ${target} = "Be" ]; then
tg=Beryllium
fi

bm=${beam}
if [ ${beam} = "piplus" ]; then
bm=pi+
fi
if [ ${beam} = "piminus" ]; then
bm=pi-
fi

printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

for (( i=0; i<${#momentum[@]}; ++i )) do

for (( j=0; j<${#secondary[@]}; ++j )) do

for (( k=0; k<${#region[@]}; ++k )) do

### echo "energy = ${energy[$i]} and secondary = ${secondary[$j]}"

sec=${secondary[$j]}
if [ ${secondary[$j]} = "piplus" ]; then
sec=pi+
elif [ ${secondary[$j]} = "piminus" ]; then
sec=pi-
fi

printf "<tests>\n" >> ${script}

# printf "<testdes>test23</testdes>\n" >> ${script}
printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>production cross section vs momentum of secondary, in different theta bins</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target}</reaction>\n" >> ${script}
printf "<beam>${bm}</beam>\n" >> ${script}
printf "<beammomentum>${momentum[$i]} GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

# printf "<target>${target}</target>\n" >> ${script}
printf "<target>${tg}</target>\n" >> ${script}


# printf "<secondary>${secondary[$j]}</secondary>\n" >> ${script}
printf "<secondary>${sec}</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>${beam}-${target}-${momentum[$i]}GeV-${secondary[$j]}-models-${region[$k]}.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>FTFP_BERT, NuBeam, NuBeam-3-3.5, QGSP_BERT</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Comment</name>\n" >> ${script}
printf "<value>NuBeam has been officially added since geant4-10-01-b01</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Comment</name>\n" >> ${script}
printf "<value>NuBeam and its variants are experimental in ${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>K+/pi+ vs momentum, in different theta bins</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target}</reaction>\n" >> ${script}
printf "<beam>${bm}</beam>\n" >> ${script}
printf "<beammomentum>${momentum[$i]} GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

# printf "<target>${target}</target>\n" >> ${script}
printf "<target>${tg}</target>\n" >> ${script}


# printf "<secondary>${secondary[$j]}</secondary>\n" >> ${script}
printf "<secondary>K+,pi+</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>${beam}-${target}-${momentum[$i]}GeV-K-pi-ratio-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>FTFP_BERT, NuBeam, NuBeam-3-3.5, QGSP_BERT</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Comment</name>\n" >> ${script}
printf "<value>NuBeam has been officially added since geant4-10-01-b01</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Comment</name>\n" >> ${script}
printf "<value>NuBeam and its variants are experimental in ${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}
 
done


printf "</testslist>\n" >> ${script}
