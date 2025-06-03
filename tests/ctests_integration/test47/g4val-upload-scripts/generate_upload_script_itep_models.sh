#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_itep_models.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

# g4version=geant4-10.00-ref07
# status=internal

testname=test47

beam=(proton piplus piminus)
target=(C U)

momentum=(1.40 5.00)

script=${testname}-itep-upload-models.xml

printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

for (( k=0; k<${#target[@]}; ++k )) do

tg=${target[$k]}
if [ ${target[$k]} = "C" ]; then
tg=Carbon
fi
if [ ${target[$k]} = "U" ]; then
tg=Uranium
fi

for (( i=0; i<${#beam[@]}; ++i )) do

bm=${beam[$i]}
if [ ${beam[$i]} = "piplus" ]; then
bm=pi+
fi
if [ ${beam[$i]} = "piminus" ]; then
bm=pi-
fi

for (( j=0; j<${#momentum[@]}; ++j )) do

mom=${momentum[$j]}
if [ ${beam[$i]} = "proton" ]; then
if [ ${momentum[$j]} = "5.00" ]; then
mom=7.50
fi
fi

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target[$k]}</reaction>\n" >> ${script}
printf "<beam>${bm}</beam>\n" >> ${script}
printf "<beammomentum>${mom} GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

printf "<target>${tg}</target>\n" >> ${script}

printf "<secondary>proton, neutron</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>${beam[$i]}-${target[$k]}-${mom}GeV-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
if [ ${mom} = "1.40" ]; then
	printf "<value>Bertini, Binary</value>\n" >> ${script}
else
	printf "<value>Bertini, FTF</value>\n" >> ${script}
fi
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done
 
done

done

printf "</testslist>\n" >> ${script}
