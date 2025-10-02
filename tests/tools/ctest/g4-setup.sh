#!/usr/bin/env bash
# common setup of environment.

#directory of script
ScriptDir=$(cd $(dirname ${BASH_SOURCE[0]});pwd)

#import common functions
source $ScriptDir/common_functions.sh

#import functions to get version numbers, like get_CLHEP_version
source $ScriptDir/versions.sh

if [ $# -ne 2 ] ; then
  echo "Error, wrong number of arguments"
  echo " arguments lcgplatform, e.g x86_64-slc5 and compiler version, e.g 4.3"
  exit 1
else
  lcgbase=$1
  shift
  gcc=$1
  shift        # do not leave args, these will be taken if we source compiler setup
fi

configGcc() {
  if IsNotInPath "gcc$gccver" ; then
    if echo ${lcgbase} | grep -e centos8 ; then
      # pick corrected compiler for centos8, in centos 8.2 cannot use compiler compiled with 8.1
      gccfullver=$(basename $(cd ${lcgcontrib}/gcc/${gccver}; /bin/pwd ))
      gcc_centos82=${lcgcontrib}/gcc/${gccfullver}.c82/${lcgbase}
      [ -d ${gcc_centos82} ] && gcc_dir=${gcc_centos82} || true
    fi
    
    [ -z "$gcc_dir" ] && gcc_dir=${lcgcontrib}/gcc/${gccver}/${lcgbase} || true
    [ -d $gcc_dir ] && true || gcc_dir=${lcgcontrib}/gcc/${gccver}/${lcgopt}
    [ -d $gcc_dir ] && true || gcc_dir=${lcgcontrib}/gcc/${gccver}/${lcgoptcentos}
    
    # as of gcc5, links may miss, so try adding .0
    if test $gcc_major_ver -ge 5 ; then
      [ -d $gcc_dir ] && true || gcc_dir=${lcgcontrib}/gcc/${gccver}.0/${lcgbase}
      [ -d $gcc_dir ] && true || gcc_dir=${lcgcontrib}/gcc/${gccver}.0/${lcgopt}
      [ -d $gcc_dir ] && true || gcc_dir=${lcgcontrib}/gcc/${gccver}.0/${lcgoptcentos}
    fi  
    
    echo Try gcc from $gcc_dir
    if [ -d $gcc_dir ] ; then
      . $gcc_dir/setup.sh
      echo "Using g++ from $gcc_dir"
    else
      echo "Error: gcc $gccver not found"  
    fi 
    unset gcc_dir    
  fi
}

configClang() {
  clangversion=$(echo $1 | awk '{ len=length($1)-6;maj=substr($1,6); printf"%d",maj }')

  # on centos (and slc6?) with clang (>=80) comnpiler, root creates shadowing warnings from enums --> disable root
  export ROOTSYS="-"

  Clang_dir=/cvmfs/sft.cern.ch/lcg/contrib/clang
  Clang_dir=$(ls -Ldtr1 ${Clang_dir}/${clangversion}/${lcgbase} | tail -1)

  # now set up clang
  .  ${Clang_dir}/setup.sh

  # .. and also tell clang about g++ used: 
  gccbase=`which g++ | sed -e "s=/bin/g++=="`
      
  # and tell CMake about the compiler
  export CXX="${Clang_dir}/bin/clang++"
  export  CC="${Clang_dir}/bin/clang"

  echo "              and g++ from: `which g++`" 
  echo "              and gcc toolchain from $gccbase"

  export gccver=$(which g++ | awk -F/ '{print substr($7,1,3)}')   # only keep maj.minor version
  echo "clang uses gcc version $gccver"
  echo THIS=${THIS}
  CTF=${THIS}/clang-toolchain.cmake
  cat > ${CTF} << EoI     
set(CMAKE_CXX_COMPILER_EXTERNAL_TOOLCHAIN "$gccbase")
set(CMAKE_C_COMPILER_EXTERNAL_TOOLCHAIN   "$gccbase")
EoI
  
  test -z "${ExtraCMakeOptions}" \
    && ExtraCMakeOptions="-DCMAKE_TOOLCHAIN_FILE=${CTF}" \
    || ExtraCMakeOptions="${ExtraCMakeOptions};-DCMAKE_TOOLCHAIN_FILE=${CTF}"
  export CMAKE_TOOLCHAIN_FILE=${CTF}      
  CMake_version=3.23.2
}

configMacClang() {
  # using clang/clang++
  cid=`clang --version | head -1`
  if echo $cid | grep -q sed ; then
    gccver=`echo $cid | awk '{ print $NF }' | sed -e"s/svn)//" | awk -F\. '{print $1"."$2}'`
  else
    gccver=`echo $cid  | awk '{ print $(NF-1) }'  | awk -F\. '{print $1"."$2}'`
  fi      

  export compiler_name=clang
  export CC=`which clang`
  export CXX=`which clang++`
}

configIntel() {
  OS=$(echo $lcgbase | awk -F- '{ print $2 }')
  case ${OS} in
    slc6)
      gccver=4.9
      ;;
    centos7)
      gccver=9
      ;;
    centos8) # to be checked...
      gccver=10   
      ;;  
    el9)     # on lxplus, native g++ fail in icpx -E -x c++ -v - < /dev/null
      gccver=11
      ;;  
    *)
      errexit "Error: OS:  ${OS} not configured for Intel compiler"
      ;;      
  esac

  #---Compiler------------------------------------------------------
  # Must first set up gcc, then icc 
  . /cvmfs/sft.cern.ch/lcg/contrib/gcc/${gccver}/${lcgbase}/setup.sh
  which g++

  case $1 in 
    icpx*)
      local oneAPIvers=2024
      ;;
    icc*)   #icc not in 2024....
      local oneAPIvers=2023
      ;;
    *)
      errexit "Intel compiler $1 not configured"
  esac

  source /cvmfs/projects.cern.ch/intelsw/oneAPI/linux/x86_64/${oneAPIvers}/setvars.sh

  case $1 in 
    icpx*)
      export CC=`which icx`
      export CXX=`which icpx`
      icpx -E -x c++ -v - < /dev/null
      ;;
    icc*)
      export CC=`which icc`
      export CXX=`which icc`
      icc -E -x c++ -v - < /dev/null
      ;;
    *)
      errexit "Intel compiler $1 not configured"
  esac      

  # CMake is looking for libs no longer part of icc19, skip this check in CMAKE
  #        test -z "${ExtraCMakeOptions}" \
  #            && ExtraCMakeOptions="-DCMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS=ON" \
  #            || ExtraCMakeOptions="${ExtraCMakeOptions};-DCMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS=ON"
  #        export ExtraCMakeOptions

  # on centos (and slc6?) with icc19 comnpiler, root triggers compilation errors --> disable root
  export ROOTSYS="-"

  #----- Force use of internal CLHEP----------------
  echo ${BUILDOPTIONS} | grep -q UseInternalCLHEP \
    && true                                     \
    || export BUILDOPTIONS="${BUILDOPTIONS}${BUILDOPTIONS:+,}UseInternalCLHEP"

  # CMake < 3.19.3 has bug when used with INTEL oneAPI, icpx needs 3.20
  CMake_version=3.20.0
}


lcgpackageroot()
{
  # return root directory of package from lcg release LCG_ver for lcgplatform
  packag=$1
  shift
  lcgver=$1
  shift
  read -r -a lcgplt <<< "$*"
  for platform in ${A_lcgopt[@]} ; do
    dir=`ls -trd1 ${lcgrel}/${LCG_ver}/${packag}/*/${platform} 2> /dev/null | tail -1`
    test -n "${dir}" && break || true
  done   
  echo "${dir}"
}

#-----------------------------------------------------------------

OS=$(echo $lcgbase | awk -F- '{ print $2 }')
compiler_name=gcc
#---Compiler-setup-----------------------------------------------------
case "$gcc" in 
  native)
    if echo "$OS" | grep -q mac; then 
      configMacClang
    else
      # docker ubuntu ccache mimics gcc in usr/local/bin/gcc and fails
      DeletePath PATH /usr/local/sbin
      DeletePath PATH /usr/local/bin
           
      [ -z "${CXX}" ] && export CXX=`which g++`
      [ -z "${CC}"  ] && export CC=`which gcc`
      gccver=`g++ --version | head -1 | awk '{ print $3 }'`
      echo "native compiler is version $gccver"
      # reduce version from 3 digit to 2
      gccver=`echo ${gccver} | awk -F. '{ printf "%d.%d",$1,$2}'`
    fi
    ;;
  clang*)
    configClang ${gcc}
    ;;
  icc*|icpx*)
    configIntel ${gcc}
    ;;
  gcc4* | gcc5* | gcc6*)
    gccver=`echo ${gcc} | awk '{ maj=substr($1,4,1); printf"%d",maj }'`
    test $gccver -lt 7 \
      && gccver=`echo ${gcc} \
       | awk '{ maj=substr($1,4,1);min=substr($1,5);printf"%d.%d",maj,min }'`
    ;;
  gcc*)
    gccver=$(echo ${gcc} | sed -e 's/gcc//')
    # full setup only after setting variables for lcg...
    ;;
  *)
    errexit "Not configured for compiler $gcc"
    ;;        
esac

#--- initial lcg setup, needed for compiler ---------------------

lcgcontrib="/cvmfs/sft.cern.ch/lcg/contrib"
lcg="${lcgbase}-${compiler_name}`echo ${gccver} | tr -d '.'`"
gcc_major_ver=`echo ${gccver} | awk -F. '{ print $1 }'`
gccver1=`echo ${gccver} | awk -F. 'NF==2 {printf"gcc%d%d",$1,($2-1>0)?($2-1):0 } NF==1 { printf"gcc%d",$1 }'`
gccver2=`echo ${gccver} | awk -F. 'NF==2 {printf"gcc%d%d",$1,($2-2>0)?($2-2):0 } NF==1 { printf"gcc%d",$1 }'`
lcgopt="${lcg}-opt"
lcgopt0="${lcgbase}-gcc${gcc_major_ver}${binutils}-opt"     #major gcc version only
lcgopt1="${lcgbase}-${gccver1}-opt"                         #reduce minor gcc version by 1, ie for gcc63 try gcc62
lcgopt2="${lcgbase}-${gccver2}-opt"                         # reduce minor gcc version by 2
lcgoptcentos=$lcgopt
echo $lcgopt | grep -e '-cc7' && lcgoptcentos=`echo $lcgopt | sed 's/-cc7/-centos7/'`

declare A_lcgopt=("$lcgopt" "$lcgoptcentos" "$lcgopt0" "$lcgopt1" "$lcgopt2")

#---Compiler-explicit-version-from-LCG---------------------------------------------------
     
if [ -z "${CXX}" ] ; then
  configGcc
  export CXX=`which g++`
  export CC=`which gcc`
fi

echo "Using Compiler:  $CXX"
$CXX --version | head -1

#--- lcg ----------------------------------------------------------
#debug.......
env | sort
lcgrel="/cvmfs/sft.cern.ch/lcg/releases"
lcgapp="/cvmfs/sft.cern.ch/lcg/app/releases"
g4lcgapp="/cvmfs/geant4.cern.ch/externals"

# clhep/XercesC locally installed in ...
# identify build dir, i.e. dir where workspace is  
jk_job=$(echo ${JOB_NAME} |  awk -F/ '{ print $1}')
build_dir=$(dirname $(echo ${WORKSPACE} | sed -e "s#/${jk_job}.*##" ))
# for a bad value, set to default 
test -z "${build_dir}"    && build_dir='/build'
test "${build_dir}" = "/" && build_dir='/build'

set +vx

declare -a g4local=("${build_dir}/externals" "/build/externals" "/ec/externals")

[ -z ${rootlcg} ] && rootlcg=${lcgopt}

LCG_ver=$(get_LCG_version)
#^^^^^^^^^^^^^-------------------------------

#---CMAKE_PREFIX_PATH----
# a search path for CMake Find..., set to "" if not set
CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:-""}
#------------------------
#===========================================================================================================================
# set OS specific things
case "$OS" in
  mac*) # fix for freetype, may miss link for /usr/X11R6 to /opt/X11 (known for macitois21
    AppendPath CMAKE_PREFIX_PATH /opt/X11
    ;;
  *)
    true
    ;;
esac

#===========================================================================================================================
#---XercesC----------------------------------------------------------
whichXercesC_LCG()
{
  for platform in ${A_lcgopt[@]}; do
    if test -d ${lcgrel}/${LCG_ver}/XercesC/*/${platform} ; then
      XercesC_dir=`ls -dt1 ${lcgrel}/${LCG_ver}/XercesC/*/${platform}`
      break
    fi
  done        
}

whichXercesC()
{
  version=$1;
  shift;
  for place in "$@"; do
    for platform in ${A_lcgopt[@]}; do
      echo "looking for XercesC in ${place}/XercesC/${version}/${platform}"
      if test -d ${place}/XercesC/${version}/${platform} ; then
        XercesC_dir=${place}/XercesC/${version}/${platform}
        break 2
      fi
    done
   done        
}

#---Xerces-C-------------------------------------------------------
if [ -z "${XERCESC_ROOT_DIR}" ] ; then
  whichXercesC_LCG
  #g4 version of  of XercesC, use as backup for compilers/builds not in lcg..
  if [ -z "${XercesC_dir}" ] ; then
    echo "looking for XercesC in ${g4lcgapp}, and in ${g4local[@]}"
    whichXercesC  3.3.0 ${g4lcgapp} ${g4local[@]}
  fi
  # try older XercesC version
  if [ -z "${XercesC_dir}" ] ; then
    echo "looking for XercesC in ${g4lcgapp}, and in ${g4local[@]}"
    whichXercesC  3.2.5 ${g4lcgapp} ${g4local[@]}
  fi
  if [ -z "${XercesC_dir}" ] ; then
    whichXercesC  3.2.3 ${g4lcgapp} ${g4local[@]}
  fi
  
  if [ -n "$XercesC_dir" -a -d "$XercesC_dir" ] ; then
    export XERCESC_ROOT_DIR=${XercesC_dir}
    AppendPath CMAKE_PREFIX_PATH ${XercesC_dir}
    AppendPath LD_LIBRARY_PATH ${XercesC_dir}/lib
    echo "Using XercesC from $XercesC_dir"
  else
    echo "Not using XercesC, XercesC_dir=${XercesC_dir:-Null}"
  fi
  unset XercesC_dir
fi

#---Inventor--------------------------------------------------------
#if  IsNotInPath coin3d ; then
#  Inventor=${lcgext}/coin3d/3.1.3.p1/${lcgopt}
#  if [ -d $Inventor ] ; then 
#    export PATH=${Inventor}/bin:${PATH}
#    export LD_LIBRARY_PATH=${lcgext}/lib:${LD_LIBRARY_PATH}
#    echo "Using Inventor in $Inventor"
#  fi
#  unset Inventor  
#fi

#---CLHEP----------------------------------------------------------
if [ -z "${CLHEP_Version}" ] ; then
  CLHEP_Version=$(get_CLHEP_version)
else 
  echo "Attempt to use requested CLHEP version ${CLHEP_Version}."    
fi

findCLHEP () {
  version=$1
  shift
  for dir in $@; do 
    for platform in ${A_lcgopt[@]}; do
      if test -d ${dir}/clhep/${version}/${platform}; then
        CLHEP_ROOT_DIR=${dir}/clhep/${version}/${platform}
        break 2
      fi  
    done   
  done   
}

#---CLHEP----------------------------------------------------------
if [ -z "${CLHEP_ROOT_DIR}" ] ; then
  findCLHEP ${CLHEP_Version} ${g4lcgapp} ${g4local[@]}
fi  
if [ -n "${CLHEP_ROOT_DIR}" -a "${CLHEP_ROOT_DIR}" != "-" -a -d "${CLHEP_ROOT_DIR}" ] ; then
  export CLHEP_ROOT_DIR
  PrependPath LD_LIBRARY_PATH ${CLHEP_ROOT_DIR}/lib
  Use_CLHEP="-DGEANT4_USE_SYSTEM_CLHEP=ON;-DCLHEP_ROOT_DIR=${CLHEP_ROOT_DIR}"
  echo "Using system CLHEP from ${CLHEP_ROOT_DIR}"
else
  echo "No CLHEP found in ${CLHEP_ROOT_DIR}"
  unset CLHEP_ROOT_DIR
  Use_CLHEP="-DGEANT4_USE_SYSTEM_CLHEP=OFF"
  echo "Using internal CLHEP"
fi
export Use_CLHEP
#------------------------------------------------------------------


whichQT()
{
  # Find 5 in prefence to 6
  for QTver in qt5 qt6; do
    for platform in ${A_lcgopt[@]}; do
      if test -d ${lcgrel}/${LCG_ver}/${QTver}/*/${platform}; then
        QTDIR=$(ls -dt1 ${lcgrel}/${LCG_ver}/${QTver}/*/${platform})
        QTVER=${QTver}
        break 2
      fi  
    done   
  done   
}

#---Qt-------------------------------------------------------------
# skip system provided qt3: 
QTDIR=''

if [ -z "${QTDIR}" ] ; then
  whichQT
fi    
if [ -n "${QTDIR}" -a "${QTDIR}" != "-" -a -d "${QTDIR}" ] ; then
  Use_QT="-DGEANT4_USE_QT=ON"
  export QTDIR
  AppendPath CMAKE_PREFIX_PATH ${QTDIR}
  if [ "${QTVER}" = "qt5" ] ; then
    export BUILDOPTIONS="${BUILDOPTIONS}${BUILDOPTIONS:+,}UseQt5"
  elif [ "${QTVER}" = "qt6" ] ; then
    export BUILDOPTIONS="${BUILDOPTIONS}${BUILDOPTIONS:+,}UseQt6"
  else
    echo "Unknown QT version ${QTVER}"
  fi
  echo "Using QT from ${QTDIR}"
else
  unset QTDIR
  unset QTVER
  echo "Not using QT"
fi
export Use_QT

#---CMake----------------------------------------------------------
test -z "$CMake_version" && CMake_version=$(get_CMake_version)
case ${OS} in
  slc6|centos*|ubuntu*|el*)  
    PrependPath PATH ${lcgcontrib}/CMake/${CMake_version}/Linux-x86_64/bin
    PrependPath CMAKE_MODULE_PATH ${lcgcontrib}/CMake/${CMake_version}/Linux-x86_64/share/cmake-*/Modules
    echo "Using LCG provided cmake/ctest: $(which cmake)"
    ;;
  mac*) 
    AppendPath PATH /usr/local/bin
    echo "Adding /usr/local/bin to PATH for CMake/CTest"
    ;; 
  *) 
    echo "Using local cmake/ctest: $(which cmake)"
    ;;
esac
cmake --version       


whichRoot()
{
  for platform in  ${A_lcgopt[@]}; do
    if [ -d ${lcgrel}/${LCG_ver}/ROOT/*/${platform} ] ; then
      root_dir=`ls -dt1 ${lcgrel}/${LCG_ver}/ROOT/*/${platform}`
      echo "found ROOT in $root_dir"
      root=${root_dir}/ROOT-env.sh
      break
    fi
  done        
}

#---ROOT-----------------------------------------------------------
# when not using a view, disable root;
echo $BUILDOPTIONS | grep -e "noView" && ROOTSYS="-" || true

if [ -z "${ROOTSYS}" ] ; then
  #  find root in given LCG release
  whichRoot
  if [ -n "${root}" -a -r "${root}" ] ; then
    echo "Setting up root using ${root}"
    . ${root}

    if test -r ${ROOTSYS}/lib/libImt.so ; then 
      if ldd ${ROOTSYS}/lib/libImt.so | grep -q libtbb ; then
        TBB_LIB_DIR=$(lcgpackageroot tbb $LCG_ver ${A_lcgplatform[@]})/lib
        AppendPath LD_LIBRARY_PATH $TBB_LIB_DIR
      fi
    fi    
    if test -r ${ROOTSYS}/lib/libROOTDataFrame.so ; then 
      if ldd ${ROOTSYS}/lib/libROOTDataFrame.so | grep -q libvdt ; then
        VDT_LIB_DIR=$(lcgpackageroot vdt $LCG_ver ${A_lcgplatform[@]})/lib
        AppendPath LD_LIBRARY_PATH $VDT_LIB_DIR
      fi
      if ldd ${ROOTSYS}/lib/libROOTDataFrame.so | grep -q libdavix ; then
        DAVIX_LIB_DIR=$(lcgpackageroot Davix $LCG_ver ${A_lcgplatform[@]})/lib64
        AppendPath LD_LIBRARY_PATH $DAVIX_LIB_DIR
      fi
      if ldd ${ROOTSYS}/lib/libROOTDataFrame.so | grep -q sqlite3 ; then
        SQLITE_LIB_DIR=$(lcgpackageroot  sqlite $LCG_ver ${A_lcgplatform[@]})/lib
        AppendPath LD_LIBRARY_PATH $SQLITE_LIB_DIR
      fi
      export LD_LIBRARY_PATH
    fi    
  else
    echo "Not using root, no setup script: ${root}"
  fi
  unset root
fi

whichValgrind()
{
  for platform in  ${A_lcgopt[@]}; do
    if test -x ${lcgrel}/${LCG_ver}/valgrind/*/${platform}/bin/valgrind ; then
      valgrind=`ls -dt1 ${lcgrel}/${LCG_ver}/valgrind/*/${platform}`
      break
    fi
  done
  #backup, if not (yet?) in LCG release, try to pick directly from lcg/releases/valgrind
  test -z "$valgrind" && valgrind=$(ls -1t ${lcgrel}/valgrind/*/${lcgopt}/bin/valgrind 2>/dev/null | head -1)    
}
#---Valgrind----------------------------------------------------------
whichValgrind
if [ -n "$valgrind" ] ; then
  vg_dir=`dirname $valgrind`
  echo "Using valgrind from $vg_dir"
  AppendPath  PATH ${vg_dir}
  VALGRIND_ROOTSUPP=""
  if [ -n "$ROOTSYS" -a -d "$ROOTSYS" ] ; then
    test -r $ROOTSYS/etc/valgrind-root.supp  \
      && export VALGRIND_ROOTSUPP=$ROOTSYS/etc/valgrind-root.supp \
      || true
  fi  
  export VALGRIND_LIB=`echo ${valgrind} | sed -e 's=/bin/valgrind=/lib/valgrind='`
  echo "  and VALGRIND_LIB = $VALGRIND_LIB"
  echo "  and Root valgrind supression file : $VALGRIND_ROOTSUPP"
else
  echo "Not using valgrind"
fi 
#----------------------------------------------------------------------

whichPython()
{
  for platform in  ${A_lcgopt[@]}; do
    if test -d ${lcgrel}/${LCG_ver}/Python/*/${platform} ; then
      pythonDir=`ls -dt1 ${lcgrel}/${LCG_ver}/Python/*/${platform}`
      echo "found python in $pythonDir"
      break
    fi
  done        
}
#---Python-------------------------------------------------------------
#Needed for PhysicsChecks: python compatible with ROOT
if IsNotInPath Python ; then
  whichPython
  if test -n "${pythonDir}" -a -d "${pythonDir}"; then
    echo "Using python from ${pythonDir}"
    PrependPath PATH ${pythonDir}/bin
    PrependPath LD_LIBRARY_PATH ${pythonDir}/lib
  else
    echo "Using Python provided by system"
  fi
  unset pythonDir
else
  echo "Using Python as found in PATH"
fi
#-----------------------------------------------------------------------

whichHDF5()
{
  for platform in  ${A_lcgopt[@]}; do
    if test -d ${lcgrel}/${LCG_ver}/hdf5/*/${platform} ; then
      HDF5DIR=`ls -dt1 ${lcgrel}/${LCG_ver}/hdf5/*/${platform}`
      break
    fi
  done        
}
#---hdf5----------------------------------------------------------------
# Search for hdf5, used by test03
if IsNotInPath h5cc ; then
  whichHDF5
  if test -x ${HDF5DIR}/bin/h5cc &&  ! echo ${HDF5DIR} | grep -qe 1\.14\. ; then
    echo "Using HDF5 from $HDF5DIR"
    AppendPath CMAKE_PREFIX_PATH ${HDF5DIR}
    export UseHDF5="-DGEANT4_USE_HDF5=ON"
  else
    echo "Not using HDF5"
  fi                    
else
  echo "Using HDF5 from PATH, h5cc is: `which h5cc`"
  export UseHDF5="-DGEANT4_USE_HDF5=ON"
fi
echo "HDF5 UseHDF5=$UseHDF5."
#-----------------------------------------------------------------------

export CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH
export PATH
