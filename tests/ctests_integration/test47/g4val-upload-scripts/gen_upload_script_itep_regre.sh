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

model=(bertini inclxx binary ftfp)

regression=(g4-9.6 g4-10.1-p03 g4-10.2.p02)

script=${testname}-itep-upload-regre.xml
if [ -e ${script} ]; then
/bin/rm ${script}
fi

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
printf "<reaction>${bm} on ${target[$j]}</reaction>\n" >> ${script}
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
if [ ${model[$k]} = "bertini" ] || [ ${model[$k]} = "inclxx" ]; then
printf "<imageblob>${beam[$i]}-${target[$j]}-${model[$k]}-regre-p1.gif</imageblob>\n" >> ${script}
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${model[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
### --->    printf "<value>g4-9.6-p04,g4-10.0-p04,g4-10.1,${g4version}</value>\n" >> ${script}
printf "<value>geant4 versions: " >> ${script}
for (( r=0; r<${#regression[@]}; ++r )) do
printf "${regression[$r]}," >> ${script}
done
printf "${g4version}</value>\n" >> ${script}
#
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
#
### ---> printf "<value>g4-9.6-p04,g4-10.00-p04,g4-10.1,${g4version}</value>\n" >> ${script}
printf "<value>geant4 versions: " >> ${script}
for (( r=0; r<${#regression[@]}; ++r )) do
printf "${regression[$r]}," >> ${script}
done
printf "${g4version}</value>\n" >> ${script}
#
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

done

# special case for Mu2e (public releases only)

targetmu2e=(Cu Pb)
modelmu2e=(bertini inclxx ftfp)

if [ ${status} = "public" ]; then

for (( i=0; i<${#targetmu2e[@]}; ++i )) do

for (( k=0; k<${#modelmu2e[@]}; ++k )) do

tgmu2e=${targetmu2e[$i]}
if [ ${targetmu2e[$i]} = "Cu" ]; then
tgmu2e=Copper
fi
if [ ${targetmu2e[$i]} = "Pb" ]; then
tgmu2e=Lead
fi

# 7.5GeV/c proton beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>proton on ${targetmu2e[$i]}</reaction>\n" >> ${script}
printf "<beam>proton</beam>\n" >> ${script}
printf "<beammomentum>7.50 GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

printf "<target>${tgmu2e}</target>\n" >> ${script}

printf "<secondary>proton, neutron</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
#
printf "<imageblob>proton-${targetmu2e[$i]}-7.50GeV-${modelmu2e[$k]}-regre.gif</imageblob>\n" >> ${script}
#
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${modelmu2e[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>geant4 versions: " >> ${script}
for (( r=0; r<${#regression[@]}; ++r )) do
printf "${regression[$r]}," >> ${script}
done
printf "${g4version}</value>\n" >> ${script}
#
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}


# 5GeV/c pi+ beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>pi+ on ${targetmu2e[$i]}</reaction>\n" >> ${script}
printf "<beam>pi+</beam>\n" >> ${script}
printf "<beammomentum>5.00 GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

printf "<target>${tgmu2e}</target>\n" >> ${script}

printf "<secondary>proton, neutron</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
#
printf "<imageblob>piplus-${targetmu2e[$i]}-5.00GeV-${modelmu2e[$k]}-regre.gif</imageblob>\n" >> ${script}
#
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${modelmu2e[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>geant4 versions: " >> ${script}
for (( r=0; r<${#regression[@]}; ++r )) do
printf "${regression[$r]}," >> ${script}
done
printf "${g4version}</value>\n" >> ${script}
#
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

# 5GeV/c pi- beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>pi- on ${targetmu2e[$i]}</reaction>\n" >> ${script}
printf "<beam>pi-</beam>\n" >> ${script}
printf "<beammomentum>5.00 GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

printf "<target>${tgmu2e}</target>\n" >> ${script}

printf "<secondary>proton, neutron</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
#
printf "<imageblob>piminus-${targetmu2e[$i]}-5.00GeV-${modelmu2e[$k]}-regre.gif</imageblob>\n" >> ${script}
#
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${modelmu2e[$k]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing</name>\n" >> ${script}
printf "<value>geant4 versions: " >> ${script}
for (( r=0; r<${#regression[@]}; ++r )) do
printf "${regression[$r]}," >> ${script}
done
printf "${g4version}</value>\n" >> ${script}
#
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

fi

printf "</testslist>\n" >> ${script}
