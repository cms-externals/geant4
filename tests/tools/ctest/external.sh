#!/bin/bash

# script driving a build of an external.
# used by job G4p-exterals-build

touch $WORKSPACE/controlfile
echo "build script running for externals...."
date
echo "running on $(hostname)"
#set -vx
/bin/pwd
env | sort
package_ver=${package}_version
package_url=${package}_URL
export VERSION=${!package_ver}
export URL=${!package_url}
export eos_install_tmp=/eos/project/g/geant4/cvmfs_install_tmp/

echo ${LABEL}
echo ${HOSTNAME}

#check for EOS being available...
if echo ${HOSTNAME} | grep -qi docker ; then
echo "will copy install files to EOS"
else
export eos_install_tmp=${WORKSPACE}/install
echo "will copy install files to ${WORKSPACE}/install"
fi   

export PATH=${PATH}:${WORKSPACE}/scripts/tests/tools/ctest
THIS=$(dirname $0)
   ${THIS}/external_download.sh
   ${THIS}/external_build.sh
   ${THIS}/external_install.sh

