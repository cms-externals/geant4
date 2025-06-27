#!/bin/bash

# work around for gcc setup script....
export SHELL=/bin/bash

ScriptDir=$(cd $(dirname $0) && /bin/pwd)
echo $ScriptDir

source $ScriptDir/common_functions.sh

lcgarea=/cvmfs/sft.cern.ch/lcg/

# building of Simplified calorimeter
build_SC=1
build_static=0
skipCMake=0
gccver=""
eos_install_tmp=/eos/project-g/geant4/cvmfs_install_tmp

#______________________________________________________________

case ${CMAKE_BUILD_TYPE} in 
  Release)
    type=opt
    ;;
  RelWithDebInfo)
    type=optdeb
    ;;    
  TestRelease)
    type=dbg
    ;;
  *)
    errexit "Error:  build type $CMAKE_BUILD_TYPE not configured"
    ;;
esac        

#build multithreaded G4
echo "$THREAD" | grep -q "MT" && build_MT="ON" || build_MT="OFF"
# for MT type will be opt-MT or dbg-MT:
[ ${build_MT} == "ON" ] && type=${type}-MT

test -z "${Version}" && errexit "Error: version to build not given"

echo "Using tag name $Branch"

#---------------------------------------------------------------------------------------
        
hw=`uname -i`
case $hw in 
  x86_64)
    bits=64
    ;;
  i*86)
    bits=32
    ;;
  *)
    errexit "unknown hardware :$hw."
    ;;
esac

lsbRelease=$(which lsb_release)
test -z "$lsbRelease" -a -r /etc/os-release && lsbRelease=$ScriptDir/lsb-release.sh

${lsbRelease} -r -s

lsb_major_version=$(${lsbRelease} -r -s | awk -F. '{ print $1 }')

case "${lsb_major_version}" in
  7)
    platform=${hw}-centos7
    ;;
  8)
    platform=${hw}-centos8
    ;;
  9)
    platform=${hw}-el9
    ;;
  *)
    errexit " not configured for: `cat /etc/issue`"
    ;;
esac

echo "=======Compiler setup for ${COMPILER}============================="

case ${COMPILER} in 
  gcc*)
    gccver=$(echo ${COMPILER} | sed -e's/gcc//')
    if echo $gccver | grep -q -e ^[7-9] ; then
      #gcc 7 and up use single digit
      gccshort=`echo ${gccver} | awk -F. '{ printf"gcc%d",$1 }'`
    elif echo "10 11 12 13 14" | grep -q -e "$gccver" ; then 
      gccshort=`echo ${gccver} | awk -F. '{ printf"gcc%d",$1 }'`
    else
      errexit "Compiler ${COMPILER} not configured"
    fi
    lcgplatform=${platform}-${gccshort}-opt
    ;;
  *)
    errexit " not configured for compiler ${COMPILER}"  
    ;;    
esac

if test -r /cvmfs/sft.cern.ch/lcg/views/${EXTERNALS}/${lcgplatform}/setup.sh ; then 
   source /cvmfs/sft.cern.ch/lcg/views/${EXTERNALS}/${lcgplatform}/setup.sh
   #unset G4 variables set by LCG... 
   unset $(env | grep ^G4 | awk -F= '{ print $1 }')
else   
   errexit "Need a LCG view to build" 
fi

echo " CC=$CC"
echo " CXX=$CXX"
echo $PATH
echo "-------------------"
date="`date  +%F-%R`"
echo "build for platform=$platform-${gccshort} in cvmfs"
echo " using local dir : $WORKSPACE" 
echo " using gcc : "
which gcc
gcc -v
echo "-------------------"

lcgrel=/cvmfs/sft.cern.ch/lcg/releases
lcgview=/cvmfs/sft.cern.ch/lcg/views
lcgcontrib=/cvmfs/sft.cern.ch/lcg/contrib
g4ext=/cvmfs/geant4.cern.ch/externals
g4install=/cvmfs/geant4.cern.ch/geant4/${Version}/${platform}-${gccshort}
g4optinstall=/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/bin

test -z "${WORKSPACE_HOST}" \
  && g4share_ref=${WORKSPACE}/geant4 \
  || g4share_ref=/geant4 

g4share_cvmfs=/cvmfs/geant4.cern.ch/geant4/${Version}/share/

echo $Version | grep -q -e ref -e cand \
  && g4share=$g4share_ref \
  || g4share=$g4share_cvmfs

if test -n "${SourceDirectory}" ; then
  echo Using source code from ${WORKSPACE}/share
  g4share=${WORKSPACE}/share
fi              

# Always use tests and verificatrion from checkout
tests_src=${g4share_ref}
verif_src=${g4share_ref}/verification
       
g4data=/cvmfs/geant4.cern.ch/share/data

g4installDir=${WORKSPACE}/install
install_opt="DESTDIR=${g4installDir}"

#-------------------------------------------------------------------------------------
# return directory of package from lcg release
rootdependency()
{
  dir=""
  rootlib=$1
  shift
  packag=$1
  shift
  lcgver=$1
  shift
  platform=$1
  if ldd ${rootlib} | grep -q ${packag} ; then
     lib=$(ldd ${rootlib} | grep ${packag}  | cut -d ' ' -f 3 - )
     dir=$(dirname $(readlink ${lib} ))
  fi
  echo "${dir}"
}

echo "=================================================================="
echo "Using Source code from ${g4share}"
ls -l ${g4share}/                
echo "=================================================================="

#---Xerces-C-------------------------------------------------------
XERCESC_ROOT_DIR=$lcgview/${EXTERNALS}/${lcgplatform}
testdir XERCESC_ROOT_DIR
export XERCESC_ROOT_DIR
PrependPath CMAKE_PREFIX_PATH ${XERCESC_ROOT_DIR}
PrependPath LD_LIBRARY_PATH ${XERCESC_ROOT_DIR}/lib

#---CLHEP----------------------------------------------------------
CLHEP_ROOT_DIR=$lcgview/${EXTERNALS}/${lcgplatform}
testdir CLHEP_ROOT_DIR
echo "Using CLHEP from ${CLHEP_ROOT_DIR}"
PrependPath LD_LIBRARY_PATH ${CLHEP_ROOT_DIR}/lib
PrependPath CMAKE_PREFIX_PATH ${CLHEP_ROOT_DIR}

#---hdf5----------------------------------------------------------------
# Search for hdf5, used by test03
HDF5DIR=$lcgview/${EXTERNALS}/${lcgplatform}
if test -x ${HDF5DIR}/bin/h5cc ; then
  echo "Using HDF5 from $HDF5DIR"
  PrependPath  CMAKE_PREFIX_PATH ${HDF5DIR} 
  export UseHDF5="-DGEANT4_USE_HDF5=ON"
else
  echo "Not using HDF5"
fi                    
echo "HDF5 UseHDF5=$UseHDF5."

echo "Using cmake from `which cmake`"
echo "qmake is found in:"
which qmake
`which qmake` -query "QT_INSTALL_PREFIX"

ls -l $g4share || errexit "Missing directory or symlink to g4share $g4share"

#env
buildDir=$WORKSPACE/$type

[ -n "$gccver" ] && buildDir=${buildDir}-$gccver
if [ -d ${buildDir} ] ; then
  echo "Warning: existing workdirectory ${buildDir} will be deleted"
  rm -rf ${buildDir}
fi
 
mkdir -p ${buildDir}

test -d ${buildDir}  && true || errexit "ERROR: Build directory ${buildDir} not found!"

cd ${buildDir}

echo "Starting cmake configure on `date +%F-%T`"

echo "build files are in ${buildDir}"

if [ $build_SC == 1 ] ; then
  # ROOT - pre-check, actual call of thisroot.sh is done after build of Geant4, just prior to build of grid executables.
  echo "Setting up ROOT for platform: ${lcgplatform}"
  if [ -z "${ROOTSYS}" ] ; then
    errexit "environment variable ROOTSYS not defined or empty" 
  else
    thisroot=${ROOTSYS}/bin/thisroot.sh
    rootdir=${ROOTSYS}
  fi   

  if [ ! -e ${thisroot} ];then
    echo "ERROR: No ROOT setup file found ($rootdir) Cannot continue"
    echo "ERROR: No ROOT setup file found ($rootdir) Cannot continue" 
    exit 1
  fi
  
  echo "Using root setup $thisroot" 
  ##### Using root -- check C++ level required by root
  std=$(grep -e "ROOT_CXX_FLAGS" ${rootdir}/cmake/ROOTConfig.cmake | \
        awk -F' |"' '{split($0,fields);for (f in fields) if ($f ~ "-std=") print substr($f,6) }')
  test -z "$std" &&   \
    std=$(grep -e "ROOT_CXX_STANDARD" ${ROOTSYS}/cmake/ROOTConfig.cmake | head -1 |\
     tr -d ')' | awk '{print $NF}' )
  echo "C++ standard used by root: ${std:=none}"

  case ${std:="none"} in          # := to protect for empty std   
    c++1z|c++17|17) # is default
      ;;
    c++20|20)
      CXXSTD="CXX20"
      ;;
    c++23|23) 
      # only with CMake >= 3.20
      CXXSTD="CXX23"
      ;;
    none)  
      # nothing to add
      ;; 
    *) 
      errexit "Invalid CXX standard given by ROOT: $std"
      ;;
  esac

  #fix up CXXSTD to be number only
  CXXSTD=$(echo ${CXXSTD} | sed -s '/CXX/s/CXX//')
fi      

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++End to create setup scripts++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if [ $skipCMake == 0 ] ; then
  exportCXXFLAGS=-DG4FPE_DEBUG

  echo "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}"

  cmakeCMD="cmake ${CXXSTD:+-DCMAKE_CXX_STANDARD=${CXXSTD}} \
    -DCLHEP_ROOT_DIR=${CLHEP_ROOT_DIR} \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_G3TOG4=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_USE_QT_QT5=ON \
    -DGEANT4_USE_XM=ON \
    -DGEANT4_USE_RAYTRACER_X11=ON \
    -DXERCESC_ROOT_DIR=${XERCESC_ROOT_DIR} \
    -DCMAKE_INSTALL_PREFIX=${g4install}-${type} \
    -DGEANT4_INSTALL_DATA=OFF \
    -DGEANT4_INSTALL_DATASETS_TENDL=ON      \
    -DGEANT4_INSTALL_DATASETS_NUDEXLIB=ON   \
    -DGEANT4_INSTALL_DATASETS_URRPT=ON      \
    -DGEANT4_BUILD_MULTITHREADED=${build_MT} \
    -DGEANT4_INSTALL_DATADIR=${g4data} \
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
    ${UseHDF5} \
    ${extraCMakeOptions} \
    ${g4share}"

  echo "CMake command used: $cmakeCMD" 
  $cmakeCMD 
fi  # skipCmake  

maxload=`cat /proc/cpuinfo | grep process | wc -l`
maxjobs=`echo "$maxload * 1.5" | bc | awk  -F \. '{ print $1 }'`
echo "Using make -j $maxjobs -l $maxload  ...."
echo "Starting make on `date +%F-%T`"

make -j $maxjobs -l $maxload   

echo "Starting make install on `date +%F-%T`"

make install $install_opt 

echo "install on finished on  `date +%F-%T`"

#===============================================================================
# + Start to create setup scripts
source ${ScriptDir}/build-createScripts.sh

# + Start to build tests
if [ $build_SC == 1 ] ; then
  [ -d ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type} ] \
    && true \
    || mkdir -p  ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}
   
  echo "softlink for thisroot: ln -s ${thisroot}  ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/thisroot.sh"
   
  test -L ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/thisroot.sh \
    && rm ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/thisroot.sh
  
  ln -s ${thisroot} ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/thisroot.sh
  # set up root  - this will reset thisroot variable
  test -z "${rootenv}" && source ${thisroot} || ${rootenv}
  echo "Using ROOT from: $(which root)" 
  
#set -vx
  #check for additional lib root needs ... not given by
  if test -r ${ROOTSYS}/lib/libImt.so ; then 
      TBB_LIB_DIR=$(rootdependency ${ROOTSYS}/lib/libImt.so libtbb $EXTERNALS ${lcgplatform})
      AppendPath LD_LIBRARY_PATH $TBB_LIB_DIR
  fi    
  
  if test -r ${ROOTSYS}/lib/libROOTDataFrame.so ; then 
      VDT_LIB_DIR=$(rootdependency ${ROOTSYS}/lib/libROOTDataFrame.so libvdt $EXTERNALS ${lcgplatform})
      AppendPath LD_LIBRARY_PATH $VDT_LIB_DIR
  
      DAVIX_LIB_DIR=$(rootdependency ${ROOTSYS}/lib/libROOTDataFrame.so libdavix $EXTERNALS ${lcgplatform})
      AppendPath LD_LIBRARY_PATH $DAVIX_LIB_DIR
   
      SQLITE_LIB_DIR=$(rootdependency ${ROOTSYS}/lib/libROOTDataFrame.so libsqlite $EXTERNALS ${lcgplatform})
      AppendPath LD_LIBRARY_PATH $SQLITE_LIB_DIR
  fi
fi      

#-------------------------------------------------------=====================================================================================================
Build_test(){
  testName=$1
  execName=$2
  srcPath=$3
 
  echo "Using sources from $srcPath and $testName"
  ls ${srcPath}/${testName}
  echo "======="
  # do we need to install source code?
  installSRC=1
  if ( test $# -ge 4 ) ; then
    test "x-$4" != "x-no" && installSRC=0 || true
  fi    

  # build the test to be installed into opt
  testbuildDir=${buildDir}/tests/ctests_integration/${testName}
  test -d ${testbuildDir} || mkdir -p ${testbuildDir}
  test -d ${testbuildDir} && cd ${testbuildDir} || errexit "No directory to build tests: ${testbuildDir}"

  echo "build files are in ${testbuildDir} "
  test -d build && rm -rf build
  mkdir build
  cd build
  echo "Compiling ${testName} binary in directory:"
  pwd 
  cmake -DCMAKE_BUILD_TYPE=TestRelease -DGeant4_DIR=${buildDir} ${srcPath}/${testName} 

  #---- fix up RPATH in link.txt....
  rpath=$(cat CMakeFiles/*.dir/link.txt | tr ' ' '\n' | grep -e 'Wl,-rpath' | awk -F, '{ print $3 }')
  #echo "Initial rpath: $rpath" 
  IFS=':' read -r -a a_rpath <<< "${rpath}"

  for index in "${!a_rpath[@]}" ; do
    if echo "${a_rpath[index]}" | grep -q ${buildDir} ; then
      # replace by installed 
      a_rpath[index]=${g4install}-${type}/lib64
    fi  
  done
  
  rpath=$(echo ${a_rpath[*]} | tr ' ' ':')
  # echo "Intermediate rpoath: $rpath "
  if echo $rpath | grep -q ROOT ; then  
    AppendPath rpath $TBB_LIB_DIR
    AppendPath rpath $VDT_LIB_DIR
    AppendPath rpath $DAVIX_LIB_DIR
    AppendPath rpath $SQLITE_LIB_DIR
  fi   
  echo "Final rpath $rpath" 
# create edit script for link.txt, to change rpath, and to add libtbb.so; without libtbb.so
#   being explicitely linked, it will not be resolved by rpath.
  cat >sedfile <<EoI
/Wl,-rpath/c \
-Wl,-rpath,$rpath
\=${ROOTSYS}/lib/libImt.so=c \
${ROOTSYS}/lib/libImt.so ${TBB_LIB_DIR}/libtbb.so
EoI
  cat CMakeFiles/*.dir/link.txt | tr ' ' '\n' | sed -f sedfile | tr '\n' ' ' > link.txt
  mv link.txt CMakeFiles/*.dir/link.txt
  rm sedfile

  #--- end RPATH

  make -j $maxjobs 
  
  if ( test $installSRC -eq 0 ); then
    # add to install area
    echo "Copy source files to install area: ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version" 
    if ( test ! -d ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/share/${testName} ); then 
      mkdir -p ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/share/${testName}
      cd ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/share/${testName}
      (cd ${srcPath} ; tar cf - ${testName}) | tar xf -
    fi
  fi

  mkdir -p  ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/bin
  cp -p ${testbuildDir}/build/${execName} \
    ${g4installDir}/cvmfs/geant4.cern.ch/opt/$Version/${platform}-${gccshort}-${type}/bin/
  rc=$?        
  echo "=============== Done with build for $testName = rc=${rc}  ============================" 
  if test ${rc} -eq 0 ; then
    echo "=============== Successful build for $testName ============================="
  else
    echo "=============== Failed build for $testName ============================="
    FailedBuilds="$FailedBuilds $testName"
  fi
}

#-------------------------------------------
if [ $build_SC == 1 ] ; then
  # Name  src  binary  source-from  install-source
  Build_test SimplifiedCalorimeter StatAccepTest ${verif_src}/SimplifiedCalorimeter   Yes
  
  if test -d ${tests_src}/tests/ctests_integration ; then
     testSRC=tests/ctests_integration
  else
     testSRC=tests
  fi   
  Build_test test30 test30 ${tests_src}/${testSRC}/ Yes
  Build_test test46 test46 ${tests_src}/${testSRC}/ Yes

  Build_test TestEm2 TestEm2 ${tests_src}/examples/extended/electromagnetic
  Build_test TestEm3 TestEm3 ${tests_src}/examples/extended/electromagnetic
  Build_test TestEm9 TestEm9 ${tests_src}/examples/extended/electromagnetic

  Build_test Hadr00 Hadr00 ${tests_src}/examples/extended/hadronic

  Build_test MSCL3          MSCL3          ${verif_src}/electromagnetic Yes
  Build_test FluctTest      FluctTest      ${verif_src}/electromagnetic Yes    
  Build_test ElecBackScat   ElecBackScat   ${verif_src}/electromagnetic Yes
  Build_test EmSamplingCalo EmSamplingCalo ${verif_src}/electromagnetic Yes
  Build_test MscHanson      MscHanson      ${verif_src}/electromagnetic Yes
  Build_test MscTest        MscTest        ${verif_src}/electromagnetic Yes
  Build_test MscTest1       MscTest1       ${verif_src}/electromagnetic Yes
  Build_test SiTest         SiTest         ${verif_src}/electromagnetic Yes
  Build_test TestDEDX2      TestDEDX2      ${verif_src}/electromagnetic Yes
  Build_test tileatlas      tileatlas      ${verif_src}/electromagnetic Yes
  Build_test Tracker        Tracker        ${verif_src}/electromagnetic Yes
  Build_test zmumu          zmumu          ${verif_src}/electromagnetic Yes
  Build_test FragTest       FragTest       ${verif_src}/hadronic        Yes
fi

if test -n "$FailedBuilds" ; then 
  echo '========================================================'
  echo "The following executables failed to build: $FailedBuilds"
  echo '========================================================'
else
  echo "All tests successfully built"
fi   

#---provide install.sh to copy install area to eos in host-----------------------------------------
set -vx 
cd $WORKSPACE
ls -l
echo "prepare script to copy to EOS"
cat > install.sh << EoI
    if test -d ${eos_install_tmp} ; then
      cd ${eos_install_tmp}; test -d ${Version} && true || mkdir ${Version}
      cd ${g4installDir}/cvmfs/geant4.cern.ch
      tar zcf - . | (cd ${eos_install_tmp}/${Version} ; dd of=${platform}-${gccshort}-${type}.tgz)  \
        && true \
        || echo "Fail to copy install file to Eos dir ${eos_install_tmp}/${Version}"
    else
      echo "Fail to copy install file to Eos dir ${eos_install_tmp}"
    fi        
EoI
chmod +x install.sh
