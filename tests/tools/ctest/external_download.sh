#!/bin/bash

# Download package, and expand source into workspace.

unset LC_CTYPE

#set -vx

errexit()
{
  echo $1
  exit 1
}

#env | sort 

if test -z "$OSname" ; then
  OSname="unix"
  if which lsb_release > /dev/null ; then 
    lsb_release -d | grep -qi ubuntu && OSname="ubuntu"
  fi
  uname -s | grep -q Darwin && OSname="mac"
fi

HWname="x86_64"

case $OSname in 
  unix) 
    export PATH=/usr/sue/bin:/usr/local/bin:/usr/bin:/bin:/usr/local/sbin:/usr/sbin
    INSTALL_TYPE=cvmfs
    ;;
  mac)
    export PATH=/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/opt/X11/bin
    INSTALL_TYPE=local
    HWname=$(uname -m)
    test "$HWname" = "arm64" && HWname="aarch64"
    ;;
  ubuntu) 
    # check if curl is available, else try to install
    curl=$(which curl 2>/dev/null)
    if test -z "$curl" ; then
      sudo apt-get -y install curl
    fi   
    INSTALL_TYPE=cvmfs
    ;;
  *)
    errexit "OS type=$OSname not configured"
    ;;
esac

#case $package in
#  XercesC)
#       URL=${XercesC_URL}
#       VERSION=${XercesC_version}
#    ;;
#  CLHEP)
#       URL=${CLHEP_URL}
#       VERSION=${CLHEP_version}
#    ;;
#    *)  errexit "Not configured for $package"
#    ;;
#esac

URL=$(echo $URL | sed -e "s/#VERSION#/${VERSION}/" )

# if there is a keytab, authenticate
[ -r /ec/conf/sftnight.keytab ] && kinit sftnight@CERN.CH -5 -V -k -t /ec/conf/sftnight.keytab

# then Copy/expand sources to workspace
tarfile=$(basename $URL)

if test "${INSTALL_TYPE}" = "cvmfs" ; then
  # prepare for cvmfs installation, copying to eos for later transfer to cvmfs
  INSTALL_basedir=/cvmfs/geant4.cern.ch/externals/${package}/${VERSION}
  G4_dir=${WORKSPACE}
  INSTALL_OPT="DESTDIR=${eos_install_dir}"

elif test "${INSTALL_TYPE}" = "local" ; then
  # install locally 
  # identify build dir, 
  jk_job=$(echo ${JOB_NAME} |  awk -F/ '{ print $1}')
  build_dir=$(dirname $(echo ${WORKSPACE} | sed -e "s#/${jk_job}.*##" ))

  test -z "${build_dir}" && errexit "Cannot extract name of build  directory from WORKSPACE"
  test "${build_dir}" = "/" && errexit "build  directory extracted from WORKSPACE must not be /"

   G4_dir=''
  if [ -d ${build_dir}/externals ] ; then
    G4_dir=${build_dir}/externals/${package}/$VERSION
  elif [ -d /ec/externals ] ; then
    G4_dir=/ec/externals/${package}/$VERSION
  elif [ -d ${build_dir} ] ; then
    mkdir ${build_dir}/externals && G4_dir=${build_dir}/externals/${package}/$VERSION
  fi
  
  if [ -z "$G4_dir" ] ; then
    errexit "No directory found to install $package"
  fi
  
  # this is not completely correct for native compiler, but helps to cut build short. It is re-checked later when final INSTALL-dir is known
  [ ! -d $G4_dir ] && mkdir -p $G4_dir
  INSTALL_basedir=${G4_dir}
else
  errexit "Invaliad type of install : ${INSTALL_TYPE}."
fi

INSTALL_dir=${INSTALL_basedir}${HWname}-${LABEL}-${COMPILER}-opt

#Check if install already exits -> FAIL
if [ -d ${INSTALL_dir} ] ; then
  echo "Warning: install directory $G4_dir already exists"
else
  if [ ! -r ${G4_dir}/src/${tarfile} ] ; then
    lock=$(which lockfile 2>/dev/null)
    if test -n "$lock" ; then 
      $lock -r 100 -s 15 ${G4_dir}/${VERSION}.LOCK && true ||  errexit "fail to create lock file"
    fi
    test ! -d ${G4_dir}/src && mkdir -p ${G4_dir}/src || true 
    if [ ! -r ${G4_dir}/src/${tarfile} ] ; then 
      curl -o ${G4_dir}/src/${tarfile} -L  $URL || errexit "fail to download ${package}"
    fi
    [ -r ${G4_dir}/${VERSION}.LOCK ] && rm -f  ${G4_dir}/${VERSION}.LOCK  || true
  else
    echo Using tarfile ${G4_dir}/src/${tarfile}:
    ls -l ${G4_dir}/src/${tarfile}
  fi
fi

pwd
ls -al 
[ -d ${VERSION} ] && rm -rf ${VERSION} || true

ls -l ${G4_dir}/src/
file ${G4_dir}/src/${tarfile}
tar zxf ${G4_dir}/src/${tarfile}
# keep dir of unpacked, i.e. pwd/top level dir in tar file
package_dir=$(tar tzf ${G4_dir}/src/${tarfile} | head -1) 

[ -r G4_env.sh ] && rm -f G4env.sh

cat << EoI  > G4env.sh
G4_dir=$G4_dir
INSTALL_basedir=$INSTALL_basedir
INSTALL_OPT=${INSTALL_OPT}
INSTALL_TYPE=${INSTALL_TYPE}
OSname=$OSname
HWname=$HWname
package=${package}
package_dir=${package_dir}
VERSION=${VERSION}
eos_install_tmp=${eos_install_tmp}
export PATH=$PATH
unset LC_CTYPE
EoI

echo "done part 1"
