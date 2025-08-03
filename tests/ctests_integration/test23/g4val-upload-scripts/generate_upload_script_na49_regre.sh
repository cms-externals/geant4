#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_na49_regre.sh proton C geant4-10.00-ref08 internal
#

beam=${1}
target=${2}
g4version=${3}
status=${4}

#g4version=geant4-10.00-ref07
#status=internal

testname=test23

momentum=(158)
secondary=(piplus piminus proton antiproton neutron)

models=(ftfp_bert NuBeam NuBeam-3-3.5 qgsp_bert)

#this part is not automated yet !
regression=(g4-9.6-p03 g4-10.00-p02)

script=${testname}-na49-${beam}-${target}-regre.xml

tg=${target}
if [ ${target} = "C" ]; then
tg=Carbon
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

for  (( m=0; m<${#models[@]}; ++m )) do

### echo "energy = ${energy[$i]} and secondary = ${secondary[$j]}"

sec=${secondary[$j]}
if [ ${secondary[$j]} = "piplus" ]; then
sec=pi+
elif [ ${secondary[$j]} = "piminus" ]; then
sec=pi-
fi

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>average multiplicity or average pT vs xF</observable>\n" >> ${script}
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
printf "<imageblob>${beam}-${target}-${momentum[$i]}GeV-${secondary[$j]}-${models[$m]}-regre.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${models[$m]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing:</name>\n" >> ${script}
printf "<value>geant4 versions: g4-9.6-p03,g4-10.00-p02,g4-10.00-ref07,${g4version}</value>\n" >> ${script}
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
 
done

printf "</testslist>\n" >> ${script}

# mv ${script} /home/g4p/pbs/g4-had-validation/results/${g4version}/test23/regre
