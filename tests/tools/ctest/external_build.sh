#!/bin/bash

source G4env.sh


#ls -la
errexit()
{
   echo $1
   exit 1
}

test -d ${package_dir} || errexit "No source code tree found"

set -vx

config_opts=""
label=$LABEL
# for alma the LCG label is el
echo $label | grep '^alma' > /dev/null && label=$(echo $label | sed -e 's/alma/el/')

echo "label: ${label}"

lcg_cvmfs=/cvmfs/sft.cern.ch/lcg/contrib

if ( test "${OSname}" = "mac" ) ; then
 # MAC
   macvers=$(sw_vers -productVersion | awk -F. '{ if ( $1 >= 11 ) print $1; else print $1$2; }')
   label="mac$macvers"
   [ "${COMPILER}" != "native" ] && errexit "On Mac, only native compiler allowed" || true
   # using clang/clang++
   cid=`clang --version | head -1`
   if echo $cid | grep -q sed ; then
      cc_version=`echo $cid | awk '{ print $NF }' | sed -e"s/svn)//" | awk -F\. '{print $1$2}'`
   else
      cc_version=`echo $cid  | awk '{ print $(NF-1) }' | awk -F\. '{ printf"%d%d",$1,$2}'`
   fi      
   export COMPILER=clang${cc_version}
   export CC=`which clang`
   export CXX=`which clang++`
   export CMAKE=cmake
   ncpu=` /usr/sbin/sysctl hw.ncpu | awk '{ print $2 }'`
else
 #linux
   case  $OSname in 
      unix) [ -d ${lcg_cvmfs}/CMake/3.14.2 ] \
          && CMAKE=${lcg_cvmfs}/CMake/3.14.2/Linux-x86_64/bin/cmake
           ;;
      ubuntu) CMAKE=cmake
           os_release=`lsb_release -r | awk '{ print $2 }' | awk -F. '{ print $1$2}'`  # produces eg. 2004
	   #os_release=`lsb_release -r | awk '{ print $2 }' | awk -F. '{ print $1}'`     # produces eg. 20
           label=${OSname}${os_release}
           ;;
        *) ;;   
   esac
   if [ "$COMPILER" == "native" ] ; then
      cc_version=`g++ --version | head -1 | awk '{ print $3 }' | awk -F\. '{print $1$2}'` 
      export COMPILER=gcc${cc_version}
   else
      case ${COMPILER} in
        gcc4* | gcc5* | gcc6*) 
	      echo $COMPILER | grep -qe binutils \
                && gccversion=`echo $COMPILER | sed -e's/gcc//'` \
                || gccversion=`echo $COMPILER |\
                  awk '{ v=substr($1,4); gsub("[0-9]","&.",v) ; print substr(v,0,length(v)-1) }'`
           ;;
        gcc*)  gccversion=`echo $COMPILER | sed -e's/gcc//'`
	   ;;
        clang*) clangversion=`echo $COMPILER | awk \
             '{ maj=substr($1,6,1);min=substr($1,7);printf"%d.%d",maj,min }'`
             case ${clangversion} in
               3.7 | 3.8 | 3.9)
	              gccversion=4.9
	            ;;
              *) errexit "Error, not configured for compiler clang ${clangver}"
                ;;
              esac
           ;;
        *) errexit "Error, not configured for compiler $1"
           ;;
      esac
      
      gcc_dir=${lcg_cvmfs}/gcc
      
      setup=${gcc_dir}/${gccversion}/x86_64-${label}/setup.sh    
      [ -r ${setup} ] && source ${setup} || errexit "Cannot find setup for compiler in ${setup}" || true
      
      if [ ${clangversion:-NO} != "NO" ] ; then
         Clang_dir=/afs/cern.ch/sw/lcg/external/llvm/${clangversion}/x86_64-${label}

         # clang  setup uses libs from gcc specific versions, see case above ...
         # .. so also tell clang about this: 
         gccbase=`which g++ | sed -e "s=/bin/g++=="`
         export CXXFLAGS="${CXXFLAGS:- } --gcc-toolchain=$gccbase"
         export  CCFLAGS="${CCFLAGS:- } --gcc-toolchain=$gccbase"

         # now set up clang
         .  ${Clang_dir}/setup.sh

         # and make sure to use the correct libs and includes....
         export CXX=${Clang_dir}/bin/clang++
         export  CC=${Clang_dir}/bin/clang
      else 
         export CC=`which gcc`
         export CXX=`which g++`
      fi
   fi

   echo " Using compiler : $CXX"
   echo "        g++ from: `which g++`" 


   [ -r /proc/cpuinfo ] && ncpu=`grep processor  /proc/cpuinfo  | wc -l`  || ncpu=1

fi
 
#re-check install dir, to catch case of native compiler or MAC with local install dir
INSTALL_dir=${INSTALL_basedir}/${HWname}-${label}-${COMPILER}-opt
if [ -d ${INSTALL_dir} ] ; then
   echo "============================================================="
   echo " Error: install directory (${INSTALL_dir}) already exists"
   echo "============================================================="
   exit 1
fi
#---------

echo ${INSTALL_dir}

test -d build && rm -rf build
mkdir build
cd build
${CMAKE} -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=${INSTALL_dir} \
                          -DCMAKE_BUILD_TYPE=Release \
                          ../${package_dir}

                # # from example on web??-Dmessage-loader=icu
make -j ${ncpu}
unset LC_CTYPE
make test

cd $WORKSPACE

cat << EoI  >> G4env.sh

INSTALL_dir=$INSTALL_dir

EoI

#  next step make install
