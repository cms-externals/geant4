#!/usr/bin/env bash

THIS=$(dirname $0)
touch controlfile

j_log=jenkins.log
j_err=jenkins.err
exec 3>&1 4>&2  1> >(/usr/bin/tee -a ${j_log} >&3) \
                2> >(/usr/bin/tee -a ${j_log} > >(/usr/bin/tee ${j_err} >&4) )

echo "host: `hostname`"
date

#-- link to Jenkins job
echo ""
echo "Jenkins Console Output ${BUILD_URL}/console"
echo "" 

echo "Dumping the full environment ---------------------------------------------------------"
                  # use special char 036 (<rs>) to avoid wrong replacements
env | sort | sed 's/:/:     /g' | tr '' '\n'
echo "--------------------------------------------------------------------------------------"

#capture environment from Jenkins
cmd="${THIS}/build-release.sh \
  -ver       \"${Release}\" \
  -tag       \"${Tag}\" \
  -buildType \"${BUILDTYPE}\" \
  -thread     \"${THREAD}\" \
  -lcgver    \"${EXTERNALS}\" \
  -compiler  \"${COMPILER}\" \
  -buildOpt  \"${BuildOptions}\" \
  -extraCMakeOpt \"${ExtraCMakeOptions}\" "
echo $cmd > setup.sh
echo $cmd
$cmd   
