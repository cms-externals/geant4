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

#beam=(proton piplus piminus)
beam=(proton)
target=(Al Be Cu Au Ta)
tgname=(Aluminum Beryllium Copper Gold Tantalum)

momentum=(10.0)

angle=(10.5 10.8)

model=(ftfp)

script=${testname}-pbarprod-upload-models.xml
if [ -e ${script} ]; then
/bin/rm ${script}
fi

if [ -e ${script} ]; then
/bin/rm -f ${script}
fi

printf "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" >> ${script}
printf "\n" >> ${script}
printf "<testslist>\n" >> ${script}
printf "\n" >> ${script}

for (( k=0; k<${#target[@]}; ++k )) do

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

for (( n=0; n<${#angle[@]}; ++n )) do

if [ ${target[$k]} = "Au" ]; then
if [ ${angle[$n]} = "10.5" ]; then
continue
fi
fi

if [ ${target[$k]} = "Ta" ]; then
if [ ${angle[$n]} = "10.8" ]; then
continue
fi
fi

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>momentum of secondary antiproton, at ${angle[$n]}deg</observable>\n" >> ${script}
printf "<reaction>${bm} on ${target[$k]}</reaction>\n" >> ${script}
printf "<beam>${bm}</beam>\n" >> ${script}
printf "<beammomentum>${mom} GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}

printf "<target>${tgname[$k]}</target>\n" >> ${script}

printf "<secondary>antiproton</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>${beam[$i]}-${target[$k]}-${mom}GeV-pbar-${angle[$n]}deg-models.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>" >> ${script}
for (( m=0; m<${#model[@]}; ++m )) do
printf "${model[$m]} " >> ${script}
done
printf "</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done
 
done

done


printf "</testslist>\n" >> ${script}
