#!/usr/bin/env bash 
echo "jk-setup.sh starting....."

# functions
errexit(){
  echo "Error: $(basename $0): $1"
  exit 1
}

#---Get the location this script (thisdir)
THIS=$(cd $(dirname ${BASH_SOURCE[0]});pwd)

#import script to get versions
source ${THIS}/versions.sh

# defaults, normally overridden from argumenents
EXTERNALS=$(get_EXTERNALS)
THREAD=Seq

# parse options
while [[ $# -gt 0 ]] ; do
  arg=$(echo $2 | sed -e's/^\"//' | sed -e's/\"$//')
  test -n "$arg" \
  && case $1 in
       -version)
         VERSION=$arg
         ;;
       -mode)
         MODE=$arg
         ;;
       -type)
         BUILDTYPE=$arg
         ;;
       -options)
         BUILDOPTIONS=$arg
         ;;
       -xoptions)  
         ExtraCMakeOptions=$arg
         ;;
       -label)
         LABEL=$arg
         ;;
       -compiler)
         COMPILER=$arg
         ;;
       -externals)
         test -n "$arg" && EXTERNALS=$arg
         ;;
       -thread)
         THREAD=$arg
         ;;
       *)    # unknown option
         errexit "Unknown option: $1"
         ;;
     esac
  shift;shift
done
touch controlfile
ARCH=$(uname -m)
#fix label: like g4-ubuntu16 -> ubuntu16
LABEL=$(echo ${LABEL} | sed -e 's/^g4-//') 

# get os name & version
OSname=''
uname -s | grep -q Darwin && OSname="mac"

lsbRelease=$(which lsb_release)
test -z "$lsbRelease" -a -r /etc/os-release && lsbRelease=$THIS/lsb-release.sh
test -z "$OSname" -a -z "$lsbRelease" && errexit "No lsb_release nor /etc/os-release found"
test -z "$OSname" && OSname=$(${lsbRelease} -d | awk '{ print $2 }')

FPE_OPT=""

#---set core file size limit
ulimit -Sc 100000000 && true || echo "failed to set core size limit"

case $OSname in
  Scientific)
    OSname="slc"
    OSvers=$(${lsbRelease}  -r -s | awk -F\. '{print $1}')
    echo "running on SLC version $(${lsbRelease}  -r -s)"
    FPE_OPT="FPE"
    ;;
  CentOS)
    OSname="centos"
    OSvers=$(${lsbRelease}  -r -s | awk -F\. '{print $1}')
    echo "running on CentOS version $(${lsbRelease}  -r -s)"
    FPE_OPT="FPE"
    ;;
  Red)
    OSname="el"
    OSvers=$(${lsbRelease}  -r -s | awk -F\. '{print $1}')
    echo "running on Red Hat Enterprise Linux version $(${lsbRelease}  -r -s)"
    FPE_OPT="FPE"
    ;; 
  AlmaLinux)
    OSname="el"
    OSvers=$(${lsbRelease}  -r -s | awk -F\. '{print $1}')
    echo "running on AlmaLinux version $(${lsbRelease}  -r -s)"
    FPE_OPT="FPE"
    ;;
  Ubuntu)
    OSname="ubuntu" 
    OSvers=$(${lsbRelease}  -r -s | tr -d '.')
    echo "running on Ubuntu version $(${lsbRelease}  -r -s)"
    FPE_OPT="FPE"
    ;;
  mac)
    OSvers=`sw_vers -productVersion | awk -F. '{ if ( $1 >= 11 ) print $1; else print $1$2; }'`
    echo "running on $ARCH, MacOS version $(sw_vers -productVersion), using label mac${OSvers}"
    ulimit -Sc 0 && true || echo "failed to set core size limit"
    ;;          
  *)
    errexit "Name of OS not known: $OSname" ;;
esac  

ulimit -a

# Add FPE to Buildoption for...
if test -n "${FPE_OPT}" ; then 
  if echo ${BUILDOPTIONS} | grep -q ${FPE_OPT} ; then
    true 
  else
    test -z "${BUILDOPTIONS}" && BUILDOPTIONS="FPE" || BUILDOPTIONS="${BUILDOPTIONS},FPE"
  fi    
fi

OS=${OSname}${OSvers}

if [ "${COMPILER}" = "native" ] ; then
  if [ "$OSname" = mac ]; then
    export CC=clang CXX=clang++
    lcgCompilerName="${CC}$(${CC} -dumpversion | awk -F. '{print $1$2}')"
  else
    gccver=`g++ --version | head -1 | awk '{ print $3 }'`
    echo "native compiler is version $gccver"
    # use major version only
    lcgCompilerName=`echo ${gccver} | awk -F. '{ printf "gcc%d",$1}'`
  fi   
else
  lcgCompilerName=$COMPILER
fi

# Set up environment...
# - check if there is a view, XercesC and clhep must be present, else use local setup script.

dir_view_qt5=/cvmfs/sft.cern.ch/lcg/views/$EXTERNALS/$ARCH-$OS-${lcgCompilerName}-opt
dir_view_qt6=/cvmfs/sft.cern.ch/lcg/views/$EXTERNALS-qt6/$ARCH-$OS-${lcgCompilerName}-opt

# Qt6 is now default in CMake configuration, so we check for qt6 first
if test -e ${dir_view_qt6} ; then
   dir_view=${dir_view_qt6}
   # We always add UseQt6 to the build options when view is used because these are
   # used to define the build name in the CDash dashboard
   echo "$BUILDOPTIONS" | grep -qe "UseQt6" && true ||   BUILDOPTIONS="${BUILDOPTIONS},UseQt6"
else
  if test -e ${dir_view_qt5} ; then
     echo "$BUILDOPTIONS" | grep -qe "UseQt5" && true ||   BUILDOPTIONS="${BUILDOPTIONS},UseQt5"
  fi
fi
echo "$BUILDOPTIONS" | grep -qe "UseQt5" && dir_view=${dir_view_qt5}
setup=${dir_view}/setup.sh

if test -r $setup  && ! ( echo $BUILDOPTIONS | grep -qe "noView" ) \
                   && test -d ${dir_view}/include/CLHEP && test -d ${dir_view}/include/xercesc ; then
  echo "Using view from $setup"
  source $setup
  
  # use clhep, expat, and zlib from view, if we find these
  test -d ${dir_view}/include/CLHEP    && Use_CLHEP="-DGEANT4_USE_SYSTEM_CLHEP=ON"
  test -d ${dir_view}/include/expat.h  && Use_EXPAT="-DGEANT4_USE_SYSTEM_EXPAT=ON"
  test -d ${dir_view}/include/zconf.h  && Use_ZLIB="-DGEANT4_USE_SYSTEM_ZLIB=ON"
  test -r ${dir_view}/bin/h5cc         && UseHDF5="-DGEANT4_USE_HDF5=ON"
  test -r ${dir_view}/bin/qmake        && Use_QT="-DGEANT4_USE_QT=ON"
   
  test -n "${CMAKE_PREFIX_PATH}" \
    && export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}:${dir_view} \
    || export CMAKE_PREFIX_PATH=${dir_view}
else
  test -d "${dir_view}" \
    && echo "No setup.sh found in view ${dir_view}, or view misses CLHEP or Xercesc" \
    || echo "No such view found: ${dir_view:-None}"
  echo "Using Geant4 setup ..."
  source ${THIS}/g4-setup.sh ${ARCH}-${OS} ${COMPILER}
fi

if test -n "$LABEL" -a -z "$CTEST_SITE"; then
  g4label=g4-$(echo $LABEL | sed -e's/lcg_docker_//')
  export CTEST_SITE=$g4label
fi
CTEST_SITE=$(echo ${CTEST_SITE} | sed -e 's/-cvmfs$//')

export CTEST_SCRIPT_DIRECTORY=${THIS}

echo "BuildOptions=$BUILDOPTIONS."

#______________________________________________________________________________
# local cleanup....
#for dir in scripts build $VERSION; do [ -d $dir ] && rm -rf $dir || true ; done

#______________________________________________________________________________
# process BuildOptions
source ${THIS}/unix-common.sh
