#/bin/csh

echo 'Run test37 at ' $HOST >! $1.out

${2} ${3}Aluminum.in >>& ${1}.out
mv Sandia.out Al_${1}.log

${2} ${3}Molybdenum.in >>& ${1}.out
mv Sandia.out Mo_${1}.log

${2} ${3}Tantalum.in >>& ${1}.out
mv Sandia.out Ta_${1}.log

${2} ${3}TaAl.in >>& ${1}.out
mv Sandia.out TaAl_${1}.log

${2} ${3}AlAuAl.in >>& ${1}.out
mv Sandia.out AlAuAl_${1}.log

${2} ${3}Beryllium.in >>& ${1}.out
mv Sandia.out Be_${1}.log

${2} ${3}Uranium.in >>& ${1}.out
mv Sandia.out U_${1}.log

${2} ${3}Silicon_15keV.in >>& ${1}.out
mv Sandia.out Si_${1}_15keV.log

${2} ${3}Silicon_20keV.in >>& ${1}.out
mv Sandia.out Si_${1}_20keV.log

${2} ${3}Silicon_30keV.in >>& ${1}.out
mv Sandia.out Si_${1}_30keV.log

${2} ${3}Silicon_40keV.in >>& ${1}.out
mv Sandia.out Si_${1}_40keV.log

${2} ${3}Silicon_50keV.in >>& ${1}.out
mv Sandia.out Si_${1}_50keV.log
