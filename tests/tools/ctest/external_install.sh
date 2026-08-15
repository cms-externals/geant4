#!/bin/bash
pwd
set -vx 

source G4env.sh

if test -r is_docker ; then
  mkdir -p ${package_dir}/install
  install_tmp=${WORKSPACE}/${package_dir}/install
  host_install=${package_dir}/install
  tardir=${package}/${VERSION}
  tarname=$(basename ${INSTALL_dir}).tgz
  cat > install.sh << EoI
cd ${host_install}
if test -d ${eos_install_tmp} ; then
  # tar cf - . | (cd ${eos_install_tmp} ; tar xf -)
  (cd ${eos_install_tmp};test -d ${tardir} && true || mkdir -p ${tardir})
  # for tar file drop initial part of install path
  cd cvmfs/geant4.cern.ch/externals
  tar zcf  ${eos_install_tmp}/${tardir}/${tarname} . 
else
  cd cvmfs/geant4.cern.ch/externals
  tar zcf ${WORKSPACE_HOST}/${package}-${VERSION}.tgz .
  echo install tar file is ${WORKSPACE_HOST}/${package}-${VERSION}.tgz
  echo "on host \$(hostname)"
fi        
EoI
else
  install_tmp=${eos_install_tmp}
fi

cd ${WORKSPACE}/build

if test "${INSTALL_TYPE}" = "cvmfs" ; then
  # prepare for cvmfs installation, copying to eos for later transfer to cvmfs
  INSTALL_OPT="DESTDIR=${install_tmp}"
elif test "${INSTALL_TYPE}" = "local" ; then
  INSTALL_OPT=""
else
  errexit "Invalid type of install : ${INSTALL_TYPE}."
fi

make install ${INSTALL_OPT} 

if test ${OSname} = "mac" ; then
  # fix name in library
  libs=$(find ${INSTALL_basedir} -type f -name \*dylib -print)
  for lib in $libs ; do
    install_name_tool -id ${lib} ${lib}
  done
fi   
     
exit $?
