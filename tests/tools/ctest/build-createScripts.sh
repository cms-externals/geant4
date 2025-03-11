#!/bin/sh

echo ${gccver}
echo ${lcgplatform}
echo ${QTDIR}
DestDir=${g4installDir}/${g4install}-${type}


# - Function to create Ba/sh setup script under CVMFS
create_sh()
{
  variant=$1
  variant_script=$2
  echo "creating setup script $DestDir/${variant}-setup.sh"

  cat > $DestDir/${variant}-setup.sh << EoI
#!/usr/bin/env bash
#
#  Script to setup environment to use CERN /cvmfs binary releases
#
# Check configuration parameters, for system and compiler
#

gccversion="${gccver}"
os="${lcgplatform}"

g++ --version | grep \$gccversion > /dev/null
if [ \$? != 0 ]
then
  echo "It looks like your compiler settings are not suitable"
  echo "The Operating system is expected to be \${os}"
  echo    "The compiler version should be g++ (GCC) \$gccversion"
  echo -n "The system reports that it is  "; g++ --version
  echo "Please set your PATH and LD_LIBRARY_PATH environment variables for this compiler"
  echo "You may use the setup script "
  echo "source /cvmfs/sft.cern.ch/lcg/contrib/gcc/\$gccversion/${platform}/setup.sh"
  echo "  to set your environment for this compiler"

else
# Geant4 Configuration parameters

   echo "Setting up the environment for Geant4 $Version for use with ${variant}"

   . ${g4install}-${type}/${variant_script}

#  extra geant4 Configuration params, not kept with CMake generated scripts

## QT
   if  test -z "\$LD_LIBRARY_PATH"  ; then
       export LD_LIBRARY_PATH=${QTDIR}/lib/
   else
       export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${QTDIR}/lib/
   fi

   export QT_QPA_PLATFORM_PLUGIN_PATH=${QTDIR}/plugins
   export QT_XKB_CONFIG_ROOT=/usr/share/X11/xkb


## HDF5

   if test -n "${HDF5DIR}" -a -x ${HDF5DIR}/bin/h5cc ; then
       export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${HDF5DIR}/lib
   fi
   
fi
EoI
}

# - Function to create T/CSH setup script under CVMFS
create_csh()
{
  variant=$1
  variant_script=$2
  variant_script_dir=$(dirname ${g4install}-${type}/${variant_script})
  echo "creating setup script $DestDir/${variant}-setup.csh"

  cat > $DestDir/${variant}-setup.csh << EoI
#!/bin/csh
#
#  Script to setup environment to use CERN /cvmfs binary releases
# 
#
# Check configuration parameters, for system and compiler
#

set gccversion="${gccver}"
set os="${lcgplatform}"

g++ --version | grep \$gccversion > /dev/null
if (\$? != 0) then
  echo "It looks like your compiler settings are not suitable"
  echo "The Operating system is expected to be \$os"
  echo    "The compiler version should be g++ (GCC) \$gccversion"
  echo -n "The system reports that it is  "; g++ --version
  echo "Please set your PATH and LD_LIBRARY_PATH environment variables for this compiler"
  echo "You may use the setup script "
  echo "source /cvmfs/sft.cern.ch/lcg/contrib/gcc/\$gccversion/${platform}/setup.csh"
  echo "  to set your environment for this compiler"

else
# Geant4 Configuration parameters

   echo "Setting up the environment for Geant4 $Version for use with ${variant}"

   source ${g4install}-${type}/${variant_script}  ${variant_script_dir}

#  extra geant4 Configuration params, not kept with CMake generated scripts

## QT
   if ( ! \${?LD_LIBRARY_PATH} ) then
       setenv LD_LIBRARY_PATH ${QTDIR}/lib/
   else
       setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:${QTDIR}/lib/
   endif

   setenv QT_QPA_PLATFORM_PLUGIN_PATH ${QTDIR}/plugins
   setenv QT_XKB_CONFIG_ROOT /usr/share/X11/xkb

## HDF5

   if ({ test -n "${HDF5DIR}" -a -x ${HDF5DIR}/bin/h5cc }) then
       setenv LD_LIBRARY_PATH \${LD_LIBRARY_PATH}:${HDF5DIR}/lib
   endif


endif
EoI
}

# - Function to create versioning file under CVMFS
create_g4versionfile()
{
  versionFile=${g4installDir}/cvmfs/geant4.cern.ch/geant4/$Version/Geant4_Version
  lockFile=${versionFile}.LOCK
  test -x /usr/bin/lockfile \
    && lockfile -s 5 -r 10 $lockFile \
    || touch $lockFile
  if test -r $lockFile -a ! -e $versionFile ; then 
    cat > $versionFile << EoI
#
# Geant4 version 
#    used by setup scripts, be cautious
#
$Version
EoI
    echo "created Geant4_Version file in ${versionFile}"
  else
    echo "Geant4_Version already exists in ${versionFile} : "
    ls -l ${versionFile}
  fi
  
  rm -f $lockFile
}

# - Create Scripts
# -- CMake
create_sh CMake bin/geant4.sh
create_csh CMake bin/geant4.csh

# -- GMake
create_sh GNUMake '/share/Geant4/geant4make/geant4make.sh'
create_csh GNUMake '/share/Geant4/geant4make/geant4make.csh'

# -- Version
create_g4versionfile
