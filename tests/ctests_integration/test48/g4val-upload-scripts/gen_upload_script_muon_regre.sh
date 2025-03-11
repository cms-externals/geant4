#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_pion_regre.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

# g4version=geant4-10.00-ref07
# status=internal

testname=test48


# NOTE: there're also Sundelin's data on S - will add later
target=(Al Si Ca Fe Ag I Au Pb S)
tgname=(Aluminum Silicon Calcium Iron Silver Iodine Gold Lead)

script=${testname}-upload-muminus-regre.xml

printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

for (( i=0; i<${#target[@]}; ++i )) do

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}

printf "<observable>production rate or kinetic energy of sseconadry neutron</observable>\n" >> ${script}
printf "<secondary>neutron</secondary>\n" >> ${script}

printf "<reaction>capture mu- on ${target[$i]}</reaction>\n" >> ${script}

printf "<beam>mu-</beam>\n" >> ${script}
printf "<beammomentum>0.0 MeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}
printf "<target>${tgname[$i]}</target>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}

printf "<imageblob>mumin-${target[$i]}-regre.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>captureUpdate (Bertini-based)</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>g4-9.6-p04,g4-10.0-p04,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}
printf "</tests>\n" >> ${script}
printf "\n" >> ${script}

done

printf "</testslist>\n" >> ${script}
