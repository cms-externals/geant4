#!/usr/bin/env bash
#
# usage:
# source generate_upload_script_harp_repre.sh <beam> <target> geant4-10.00-ref08 internal
#

beam=${1}
target=${2}
g4version=${3}
status=${4}

#g4version=geant4-10.00-ref07
#status=internal

testname=test23

momentum=(3.0 5.0 8.0 12.0)
secondary=(piplus piminus)
region=(FW LA)
# models=(ftfp_bert NuBeam NuBeam-3-3.5)
models=(ftfp_bert NuBeam)
modelsmu2e=(ftfp_bert qgsp_bert Shielding ShieldingM)

# this part is not automated yet !
#
#regression=(g4-9.6-p03 g4-10.00-p03 g4-10.01)

script=${testname}-harp-${beam}-${target}-regre.xml

tg=${target}
if [ ${target} = "C" ]; then
tg=Carbon
fi
if [ ${target} = "Be" ]; then
tg=Beryllium
fi
if [ ${target} = "Ta" ]; then
tg=Tantalum
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

if [ ${target} != "Ta" ]; then

for (( i=0; i<${#momentum[@]}; ++i )) do

for (( j=0; j<${#secondary[@]}; ++j )) do

for (( k=0; k<${#region[@]}; ++k )) do

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
printf "<observable>production cross section vs momentum of secondary</observable>\n" >> ${script}
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
printf "<imageblob>${beam}-${target}-${momentum[$i]}GeV-${secondary[$j]}-${region[$k]}-${models[$m]}-regre.gif</imageblob>\n" >> ${script}

printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${models[$m]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing:</name>\n" >> ${script}
printf "<value>geant4 versions: g4-9.6-p04,g4-10.1-p03,g4-10.2-p02,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Comment</name>\n" >> ${script}
printf "<value>NuBeam has been officially added since geant4-10-01-b01</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
# printf "<tag>\n" >> ${script}
# printf "<name>Comment</name>\n" >> ${script}
# printf "<value>NuBeam and its variants are experimental in ${g4version}</value>\n" >> ${script}
# printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

done
 
done

elif [ ${target} = "Ta" ]; then

for (( j=0; j<${#secondary[@]}; ++j )) do

sec=${secondary[$j]}
if [ ${secondary[$j]} = "piplus" ]; then
sec=pi+
elif [ ${secondary[$j]} = "piminus" ]; then
sec=pi-
fi

for  (( m=0; m<${#modelsmu2e[@]}; ++m )) do

for (( k=0; k<${#region[@]}; ++k )) do

### echo "energy = ${energy[$i]} and secondary = ${secondary[$j]}"

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>production cross section vs momentum of secondary</observable>\n" >> ${script}
printf "<reaction>proton on Ta</reaction>\n" >> ${script}
printf "<beam>proton</beam>\n" >> ${script}
printf "<beammomentum>8.0 GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}
printf "<target>${tg}</target>\n" >> ${script}
printf "<secondary>${sec}</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>proton-Ta-8.0GeV-${secondary[$j]}-${region[$k]}-${modelsmu2e[$m]}-regre.gif</imageblob>\n" >> ${script}
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${modelsmu2e[$m]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing:</name>\n" >> ${script}
printf "<value>geant4 versions: g4-9.6-p04,g4-10.1-p03,g4-10.2-p02,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

# theta-spectra

printf "<tests>\n" >> ${script}

printf "<testdes>${testname}</testdes>\n" >> ${script}

printf "<g4version>${g4version}</g4version>\n" >> ${script}
printf "<observable>production cross section vs polar angle of secondary</observable>\n" >> ${script}
printf "<reaction>proton on Ta</reaction>\n" >> ${script}
printf "<beam>proton</beam>\n" >> ${script}
printf "<beammomentum>8.0 GeV/c</beammomentum>\n" >> ${script}
printf "<beamenergy></beamenergy>\n" >> ${script}
printf "<target>${tg}</target>\n" >> ${script}
printf "<secondary>${sec}</secondary>\n" >> ${script}
printf "<status>${status}</status>\n" >> ${script}
printf "<score>\n" >> ${script}
printf "<type>expert</type>\n" >> ${script}
printf "<value>passed</value>\n" >> ${script}
printf "</score>\n" >> ${script}
printf "<imageblob>proton-Ta-8.0GeV-${secondary[$j]}-theta-spectra-${modelsmu2e[$m]}-regre.gif</imageblob>\n" >> ${script}
printf "<taglist>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Model</name>\n" >> ${script}
printf "<value>${modelsmu2e[$m]}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "<tag>\n" >> ${script}
printf "<name>Regression testing:</name>\n" >> ${script}
printf "<value>geant4 versions: g4-9.6-p04,g4-10.1-p03,g4-10.2-p02,${g4version}</value>\n" >> ${script}
printf "</tag>\n" >> ${script}
printf "</taglist>\n" >> ${script}

printf "</tests>\n" >> ${script}

printf "\n" >> ${script}

done

done

fi

printf "</testslist>\n" >> ${script}

# mv ${script} /home/g4p/pbs/g4-had-validation/results/${g4version}/test23/regre
