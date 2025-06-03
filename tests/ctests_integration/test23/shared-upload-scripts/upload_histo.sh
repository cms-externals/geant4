#!/usr/bin/env bash
#

#  Example usage:
#
# ./upload_histo.sh 19 geant4-10.1-p02 public <NameOfDB> <Port> <Password>

if [ "x" == "x$G4INSTALL" ]; then
echo " env. varianble G4NSTALL is NOT defined; exit"
exit 3
fi

if [ "x" == "x$ROOTSYS" ] ; then
source /g4/g4p/pbs/g4-had-validation/build-scripts/g4_set_prod.sh
fi

testid=${1}
g4version=${2}
status=${3}
# ${ROOTSYS}/bin/root -b -p -q test19_gen_histo_sql.C\(\"${g4version}\",\"${status}\"\)
rootscriptloc=${G4INSTALL}/tests/ctests_integration/test${testid}/g4val-upload-scripts
${ROOTSYS}/bin/root -b -p -q ${rootscriptloc}/test${testid}_gen_histo_sql.C\(\"${g4version}\",\"${status}\"\)

dbhost=${4}.fnal.gov
dbport=${5}
dbpasswd=${6}

# echo " dbhost = ${dbhost} "
# echo " dbport = ${dbport} "

export PGPASSWORD=${6}
# echo " PGPASSWORD = ${PGPASSWORD} "

sqlscript=test${testid}-histo.sql

# NOTE-1: -W will force asking for DB password (for the user 'g4valwriter', that is)
#         psql will ask for the password interactively
#
# --> psql -d g4validation -h fnalpgsdev.fnal.gov -p 5453 -U g4valwriter -W < test19-na49-histo.sql
#
# NOTE-2: this way it'll work WITHOUT asking for password if PGPASSWORD env.var. is set
#         ... although there's a word that PGPASSWORD option is depricated 
#             and might disappear at some point in time 
#
/usr/bin/psql -d g4validation -h ${dbhost} -p ${dbport} -U g4valwriter < ${sqlscript}
#
# NOTE-3: can NOT find 'expect' (and 'send') on tev.fnal.gov
#         in general, to run 'expect' one must include #!/usr/bin/expect header (maybe with the -f key)
#
#/usr/bin/psql -d g4validation -h ${dbhost}.fnal.gov -p ${dbport} -U g4valwriter -W < ${sqlscript}
#expect "*assword*:*" 
#{
#   send "${dbpasswd}\r"
#}

unset PGPASSWORD
# echo " PGPASSWORD = ${PGPASSWORD} "

/bin/rm ${sqlscript}

echo "Tutto bene !"


