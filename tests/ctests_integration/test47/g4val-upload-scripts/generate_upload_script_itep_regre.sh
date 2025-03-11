#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_itep_regre.sh <version> <status>
# for example:
# source generate_upload_script_itep_regre.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

testname=test47

beam=(proton piplus piminus)
target=(C U)

model=(bertini binary ftfp)

#this part is still NOT autiomated !!!
regression=(g4-9.6-p03 g4-10.00-p02 g4-10.00-ref07 )

script=${testname}-itep-upload-regre.xml


printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

for (( i=0; i<${#beam[@]}; ++i )) do

bm=${beam[$i]}
if [ ${beam[$i]} = "piplus" ]; then
bm=pi+
fi
if [ ${beam[$i]} = "piminus" ]; then
bm=pi-
fi

for (( j=0; j<${#target[@]}; ++j )) do

tg=${target[$j]}
if [ ${target[$j]} = "C" ]; then
tg=Carbon
fi
if [ ${target[$j]} = "U" ]; then
tg=Uranium
fi

for (( k=0; k<${#model[@]}; ++k )) do

mom=1.40
if [ ${model[$k]} = "ftfp" ]; then
if [ ${beam[$i]} = "proton" ]; then
mom=7.50
else
mom=5.00
fi
fi

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target}</reaction>\n" >> ${script}
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
if [ ${model[$k]} = "bertini" ]; then
printf "<imageblob>${beam[$i]}-${target[$j]}-${model[$k]}-regre-p1.gif</imageblob>\n" >> ${script}
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${model[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>g4-9.6-p03,g4-10.00-p02,g4-10.00-ref07,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}
printf "</tests>\n" >> ${script}
printf "\n" >> ${script}
printf "<tests>\n" >> ${script}
printf "<testdes>${testname}</testdes>\n" >> ${script}
printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target}</reaction>\n" >> ${script}
printf "<beam>${bm}</beam>\n" >> ${script}
mom1=${mom}
if [ ${beam[$i]} = "proton" ]; then
mom1=7.50
else
mom1=5.00
fi
printf "<beammomentum>${mom1} GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}
printf "<target>${tg}</target>\n" >> ${script}
printf "<secondary>proton, neutron</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>${beam[$i]}-${target[$j]}-${model[$k]}-regre-p2.gif</imageblob>\n" >> ${script}
else
printf "<imageblob>${beam[$i]}-${target[$j]}-${model[$k]}-regre.gif</imageblob>\n" >> ${script}
fi
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${model[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>g4-9.6-p03,g4-10.00-p02,g4-10.00-ref07,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

done

printf "</testslist>\n" >> ${script}
