#!/usr/bin/env bash

# Setup 

#directory of script
ScriptDir=$(cd $(dirname ${BASH_SOURCE[0]});pwd)

#import common functions
source $ScriptDir/common_functions.sh

test -n "$VERSION_BENCHMARKS" && echo "Using benchmarks tag: $VERSION_BENCHMARKS"
 
# jk uses workspace, and use/adapt BUILDDIR
# docker has source etc mounted as ReadOnly, so ...
if [ -n "${WORKSPACE_HOST}" -a -d "/geant4" ]; then
  export SOURCE=/geant4
  export BINARY=${WORKSPACE}
  echo "Docker build using source from ${SOURCE} and building in ${BINARY}."  
elif [ -n "${WORKSPACE}" -a -d "${WORKSPACE}" ]; then
  BUILDDIR=`/bin/pwd`
  export SOURCE=$(/bin/pwd)/geant4
  export BINARY=$(/bin/pwd)/build
  echo "Using jenkings workspace .${WORKSPACE}. and  build dir: .${BUILDDIR}."
else 
  echo "Error: No Directory to build in, using BUILDDIR=${BUILDDIR}"
fi

#
#default, to be modified for non default buildtypes
#
if [ -z "$G4_XOPTS" ]; then
  G4_XOPTS="-DGEANT4_USE_G3TOG4=ON;-DGEANT4_USE_RAYTRACER_X11=ON;-DGEANT4_INSTALL_DATA=OFF"

  case $OSname in
    mac)
      test "$OSvers" -ge 1014 && G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_OPENGL_X11=ON"
      ;;    
    slc|centos|ubuntu|el) 
      # do not enable GL, nor xm for centos 8
      if test "${OSvers}" != "8" ; then
        G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_OPENGL_X11=ON;-DGEANT4_USE_XM=ON"
      fi   
      ;;
    *)
      errexit  "not configured for OS $OSname"
      ;;
  esac    
fi

# number of cores n_cpus
test -r /proc/cpuinfo && n_cpus=$(grep processor  /proc/cpuinfo  | wc -l)   || n_cpus=""

#add optional setting from geant4 setup.
test -n "${Use_CLHEP}" && G4_XOPTS="${G4_XOPTS};${Use_CLHEP}"
test -n "${Use_EXPAT}" && G4_XOPTS="${G4_XOPTS};${Use_EXPAT}"
test -n "${Use_ZLIB}"  && G4_XOPTS="${G4_XOPTS};${Use_ZLIB}"
test -n "${Use_QT}"    && G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_QT=${Use_QT}"
test -n "${UseHDF5}"   && G4_XOPTS="${G4_XOPTS};${UseHDF5}"

#Always make TENDL data available - if we find Tendl in G4DatasetDefinitions.cmake
dsdefs="${SOURCE}/cmake/Modules/G4DatasetDefinitions.cmake"
if grep -q -e "G4TENDL" ${dsdefs} ; then  
  G4_XOPTS="${G4_XOPTS};-DGEANT4_INSTALL_DATASETS_TENDL=ON"
  needParticleHP=""
else
  needParticleHP="needed"
fi

# when using ROOT, must use -std used by root  
# -- set here, and not in g4-setup.sh,  as this also applies when using view
if test -n "$ROOTSYS" -a -d "$ROOTSYS" ; then
  std=$(grep -e "ROOT_CXX_FLAGS" ${ROOTSYS}/cmake/ROOTConfig.cmake | \
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
      # nothing to add, use G4 default
      ;; 
    *)
      errexit "Invalid CXX standard given by ROOT: $std"
      ;;
  esac            

  if echo $BUILDOPTIONS | grep -q CXX ; then 
    echo "Buildoption has CXX option already, override to use $CXXSTD, i.e. -std option used by root: -std=$std"
    BUILDOPTIONS=$(echo $BUILDOPTIONS | sed -e "s/CXX../$CXXSTD/")
  else   
    [ -z "${BUILDOPTIONS}" ] \
      && BUILDOPTIONS="${CXXSTD}" \
      || BUILDOPTIONS="$BUILDOPTIONS,${CXXSTD}"
  fi   
fi

#
# parse Buildoptions
#
echo "BuildsOptions: ${BUILDOPTIONS}"
for option in ${BUILDOPTIONS//,/ } ; do
  case "${option}" in 
    staticlibs)
      G4_XOPTS="${G4_XOPTS};-DBUILD_STATIC_LIBS=ON;-DBUILD_SHARED_LIBS=OFF"
      export CMAKE_BUILD_PARALLEL_LEVEL=$(reduce_MAX_CPUS_USE $(( $n_cpus /  4 )) )
      echo "set CMAKE_BUILD_PARALLEL_LEVEL = ${CMAKE_BUILD_PARALLEL_LEVEL}, to reduce load in building"
      export MAX_CPUS_USE=${CMAKE_BUILD_PARALLEL_LEVEL}
      ;;
    EpCheck)
      export G4Hadronic_epReportLevel=-3
      ;;
    BoundsCheck)
      export CXXFLAGS="${CXXFLAGS:- } -DG4FPE_DEBUG -D_GLIBCXX_DEBUG"
      export CTEST_TIMEOUT=5000
      G4_XOPTS="${G4_XOPTS};-DCMAKE_DISABLE_FIND_PACKAGE_ROOT=ON"
      echo ${THREAD} | grep -q MT && export MAX_CPUS_USE=$(reduce_MAX_CPUS_USE 2) || true
      echo ${THREAD} | grep -q MT && echo "setting number of tests run in parallel to ${MAX_CPUS_USE}."
      ;;
    FPE)
      export CXXFLAGS="${CXXFLAGS:- } -DG4FPE_DEBUG"
      ;;
    MemCheck)
      echo "Using Memory checking"
      export WITH_MEMCHECK=1
      ;;
    DisableVerbose)
      echo "Disabling verbose code (no definition of G4VERBOSE)" 
      G4_XOPTS="${G4_XOPTS};-DGEANT4_BUILD_VERBOSE_CODE=OFF"
      ;;
    GranularCheck)
      G4_XOPTS="${G4_XOPTS};-D__GEANT4_LIBRARY_DEFINITION_FILE=${SOURCE}/tests/code_checks/GranularCheck.cmake"
      ;;
    UseUsolids|UseVecGeom)
      G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_USOLIDS=ON"
      ;;
    VecGeom_DIR*)
      echo "VecGeomDIR = ${option}"
      VecGeomDir=$(echo $option | awk -F= '{ print $2 }')
      G4_XOPTS="${G4_XOPTS};-D${option}"
      ;;     
    UseInternalCLHEP)
      echo "Switching to internal CLHEP, as per BUILDOPTIONS"
      echo ${G4_XOPTS} | grep -q GEANT4_USE_SYSTEM_CLHEP \
        && G4_XOPTS=`echo ${G4_XOPTS} | sed -e 's/GEANT4_USE_SYSTEM_CLHEP=ON/GEANT4_USE_SYSTEM_CLHEP=OFF/'` \
        || G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_SYSTEM_CLHEP=OFF"
      ;;
    UseGranularCLHEP)
      if ( echo ${G4_XOPTS} | grep -q GEANT4_USE_SYSTEM_CLHEP ) ; then
        G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_SYSTEM_CLHEP_GRANULAR=ON"
      else
        echo "Not using system CLHEP, therefore ignoring option to use granular CLHEP"
      fi
      ;;
    UseTBB)
      echo "Using TBB switch GEANT4_USE_TBB"
      G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_TBB=ON"
      ;;            
    UseQt6)
      G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_QT_QT6=ON"   # view with Qt6 is set in jk-setup.sh.
      ;;
    UseVTK)
      G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_VTK=ON"
      ;;
    CXX17)
      echo "CXX17 is default, ignoring superfluous option"
      ;;
    CXX20)
      G4_XOPTS="${G4_XOPTS};-DCMAKE_CXX_STANDARD=20"
      echo " Compiling against c++20"
      ;;
    CXX23)
      G4_XOPTS="${G4_XOPTS};-DCMAKE_CXX_STANDARD=23"
      echo " Compiling against c++23"
      ;;
    UseLTO)
      G4_XOPTS="${G4_XOPTS};-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON"
      echo " Compiling with LTO enabled"
      ;;
    TestExampleGuidelines)
      G4_XOPTS="${G4_XOPTS};-DGEANT4_TEST_EXAMPLE_GUIDELINES=ON"
      echo " Enabling Tests of Example Coding Guidelines"
      # This is currently experimental, so we override the MODE (which is really the CTest track under Nightly)
      export MODE="experimental"
      echo " - Example Guideline tests are experimental so overriding MODE to 'experimental'"
      ;;
    DownloadData)
      echo "Turn on download of data files"
      echo ${G4_XOPTS} | grep -q GEANT4_INSTALL_DATA \
        && G4_XOPTS=`echo ${G4_XOPTS} | sed -e 's/GEANT4_INSTALL_DATA=OFF/GEANT4_INSTALL_DATA=ON/'` \
        || G4_XOPTS="${G4_XOPTS};-DGEANT4_INSTALL_DATA=ON"
      ;;
    #Options for RunManagerType 
    RM_Serial)
      echo "Use Serial run manager"
      export G4FORCE_RUN_MANAGER_TYPE="Serial"
      ;;
    RM_MT)
      echo "Use MT run manager"
      export G4FORCE_RUN_MANAGER_TYPE="MT"
      ;;
    RM_Tasking)
      echo "Use Tasking run manager"
      export G4FORCE_RUN_MANAGER_TYPE="Tasking"
      #reduce_MAX_CPUS_USE $(( $n_cpus /  4 ))
      export MAX_CPUS_USE=1
      ;;
    RM_TBB)
      echo "Use  TBB run manager"
      export G4FORCE_RUN_MANAGER_TYPE="TBB"        
      G4_XOPTS="${G4_XOPTS};-DGEANT4_USE_TBB=ON"
      ;;
    # options to ignore
    none)
      true   #nothing to do for none
      ;;
    noView)
      true   #nothing to do for noView, used in jk-setup.sh to turn off using view
      ;;
    *)
      errexit "Error: unknown option passed in BUILDOPTIONS: ${option}"
      ;;
  esac
done

# Mark if we implictely use CLHEP, when not found in lcg...
if echo ${BUILDOPTIONS} | grep -e 'UseInternalCLHEP' > /dev/null ; then
  echo "use internal CLHEP"
else
  if echo ${G4_XOPTS} | grep -e  USE_SYSTEM_CLHEP=OFF > /dev/null ; then
    BUILDOPTIONS="${BUILDOPTIONS},UseInternalCLHEP"
    echo "marking build as using internal CLHEP"
  fi    
fi

# If using existing data, check these exist, and then specify location:
if echo ${G4_XOPTS} | grep -e 'GEANT4_INSTALL_DATA=OFF' > /dev/null ; then
  g4data=/cvmfs/geant4.cern.ch/share/data
  if test -r $g4data ; then
    G4_XOPTS="${G4_XOPTS};-DGEANT4_INSTALL_DATADIR=${g4data}"
    echo "Using existing data from $g4data"
    test -n "${needParticleHP}" && export G4PARTICLEHPDATA=${g4data}/G4TENDL1.3.1
  else
    G4_XOPTS=`echo ${G4_XOPTS} | sed -e 's/GEANT4_INSTALL_DATA=OFF/GEANT4_INSTALL_DATA=ON/'`
    echo "Local data not found, requesting download"
    BUILDOPTIONS="${BUILDOPTIONS},DownloadData"
  fi
fi

# when using VecGeom, make sure VecGeom_Dir is set
if echo ${G4_XOPTS} | grep -e 'GEANT4_USE_USOLIDS=ON' > /dev/null ; then
  if echo ${G4_XOPTS} | grep -e 'VecGeom_DIR' > /dev/null ; then
    echo using VecGeom_DIR=${VecGeomDir}
  else
    VecGeomDir="${WORKSPACE}/build/VecGeom/install/lib64/cmake/VecGeom"
    echo using default VecGeom_DIR=${VecGeomDir}
    G4_XOPTS="${G4_XOPTS};-DVecGeom_DIR=${VecGeomDir}"
  fi
fi         

# When using MemCheck, must use ROOT suppressions file
if [ -n "${WITH_MEMCHECK}" ] ; then
  VALGRIND_ROOTSUPP=""
  if [ -n "$ROOTSYS" -a -d "$ROOTSYS" ] ; then
    test  -r $ROOTSYS/etc/valgrind-root.supp  \
      && export VALGRIND_ROOTSUPP=$ROOTSYS/etc/valgrind-root.supp \
      || true
  fi

  echo "Valgrind MemCheck will use suppressions file: ${VALGRIND_ROOTSUPP}"
fi

#reformat for CMake for use in build name 
for option in ${BUILDOPTIONS//,/ } ; do
  case "${option}" in
    VecGeom_DIR*)
      echo ${option} | grep -qi vecgeom && BuildName="${BuildName}${BuildName+-}VecGeom"  || true  
      ;;
    none)
      true # skip over this option, do not add to name 
      ;;  
    *)
      BuildName="${BuildName}${BuildName+-}${option}"
      ;;
  esac
done    
export BUILDOPTIONS=${BuildName}

echo " THREAD=${THREAD}."
case "${THREAD}" in 
  Seq)
    G4_XOPTS="${G4_XOPTS};-DGEANT4_BUILD_MULTITHREADED=OFF"
    #override tasking default  
    export G4FORCE_RUN_MANAGER_TYPE="Serial"
    ;;
  MT)
    G4_XOPTS="${G4_XOPTS};-DGEANT4_BUILD_MULTITHREADED=ON"
    ;;
  MTmax)
    G4_XOPTS="${G4_XOPTS};-DGEANT4_BUILD_MULTITHREADED=ON"
    if [ -n "${n_cpus}" ] ; then 
      #threads=$((15 * ${n_cpus} / 10))
      threads="max"
    else 
      threads="max"
    fi    
    export G4FORCENUMBEROFTHREADS="$threads"
    ;;
  *) 
    echo "Error: unknown option passed in THREAD: ${THREAD}"
    exit 1;
    ;;
esac

if test -n "$ExtraCMakeOptions" ; then
  echo "Adding extra CMake OPtions : $ExtraCMakeOptions"
  G4_XOPTS="${G4_XOPTS};$ExtraCMakeOptions"
fi
    
export G4_XOPTS
echo "Geant4 CMake options: ${G4_XOPTS}"
