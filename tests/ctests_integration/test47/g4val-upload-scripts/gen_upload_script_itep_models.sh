#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_itep_models.sh geant4-10.00-ref08 internal
#

g4version=${1}
status=${2}

testname=test47

beam=(proton piplus piminus)
target=(C U)

momentum=(1.40 5.00)

script=${testname}-itep-upload-models.xml
if [ -e ${script} ]; then
/bin/rm ${script}
fi

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
	printf "<value>Bertini, Binary, INCL++</value>\n" >> ${script}
else
	printf "<value>Bertini, FTF, INCL++</value>\n" >> ${script}
fi
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done
 
done

done

# special cases for Mu2e (public releases only)
#
targetmu2e=(Cu Pb)

if [ ${status} = "public" ]; then

for (( k=0; k<${#targetmu2e[@]}; ++k )) do

tgmu2e=${targetmu2e[$k]}
if [ ${targetmu2e[$k]} = "Cu" ]; then
tgmu2e=Copper
fi
if [ ${targetmu2e[$k]} = "Pb" ]; then
tgmu2e=Lead
fi

# 7.5GeV/c proton beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>proton on ${targetmu2e[$k]}</reaction>\n" >> ${script}
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
printf "<imageblob>proton-${targetmu2e[$k]}-7.50GeV-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>Bertini, FTF, INCL++</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

# 5GeV/c pi+ beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>pi+ on ${targetmu2e[$k]}</reaction>\n" >> ${script}
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
printf "<imageblob>piplus-${targetmu2e[$k]}-5.00GeV-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>Bertini, FTF, INCL++</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

# 5GeV/c pi- beam

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>kinetic energy of secondary proton or neutron</observable>\n" >> ${script}
printf "<reaction>pi- on ${targetmu2e[$k]}</reaction>\n" >> ${script}
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
printf "<imageblob>piminus-${targetmu2e[$k]}-5.00GeV-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>Bertini, FTF, INCL++</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

fi

printf "</testslist>\n" >> ${script}
