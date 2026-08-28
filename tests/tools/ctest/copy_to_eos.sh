#!/usr/bin/env bash

#set -vx

errexit(){
  echo "Error: $(basename $0): $1"
  exit 1
}

if test -d ${WORKSPACE}/build ; then 
  echo "Will copy build area to EOS"
else
  echo "no build dir found, quit..."
  exit 0
fi
   
cd ${WORKSPACE}/build

IFS=',' read -a jbn <<< $JOB_BASE_NAME

declare -a jbo

jbo="BUILDTYPE=${BUILDTYPE}"
jbo="BuildOptions=${BuildOptions} ${jbo[@]}"
jbo="DAY=$(LC_ALL=en_US.UTF-8 date +%a) ${jbo[@]}"
jbo="${jbn[@]} ${jbo[@]}"

declare -A options

for v in ${jbo[@]} ; do 
  IFS='=' read -a kv <<< $v
  echo ${kv[0]} = ${kv[1]}
  options[${kv[0]}]=${kv[1]}
done

eosDir=/eos/project/g/geant4/nightlies
eosDir=$(cd $eosDir; /bin/pwd)
cd $eosDir || errexit "Directory for nightlies in EOS not available"

srcDir=${options["DAY"]}/geant4
test -d $srcDir || mkdir -p $srcDir

LOCK=LOCKED

# remove old lock file
if test -r ${LOCK} ; then
  if [ `stat --format=%Y ${LOCK}` -le $(( `date +%s` - 3600 )) ]; then    
    rm -f $LOCK
  fi
fi      

# if there is no lockfile use true....
test -x /usr/bin/lockfile && lockcmd="lockfile -10 -r 10 ${LOCK}" || lockcmd="true"

if $lockcmd ; then
  if test -f ${srcDir}/CREATED ; then 
    # the following does a basic check, but if we re-build within 6 hours the old source code will be kept.
    if [ `stat --format=%Y ${srcDir}/CREATED` -le $(( `date +%s` - 6 * 3600 )) ]; then 
      echo "remove old source dir $srcDir"
      rm -rf $srcDir
    fi
  fi

  test -d $srcDir || mkdir -p $srcDir
  cd $srcDir || errexit "Directory for source code, $srcDir, not available, or not created"
  if ! test -f CREATED ; then
    touch CREATED
    (cd ${WORKSPACE}/geant4 ; tar cf - . ) | tar xf -
    touch CREATED
  fi
   
  rm -f ${eosDir}/${LOCK}
fi

cd $eosDir || errexit "Directory for nightlies in EOS not available"

IFS='/' read -a dirs <<< $(echo $WORKSPACE) 
# part of WORKSPACE following $JOB_NAME
dir_work=""
for ((i= ${#dirs[@]}-1; $i;--i)) ; do 
  echo ${dirs[$i]};item=${dirs[$i]};test "${item}" == "${JOB_NAME}" \
    && break \
    || dir_work="$item/$dir_work" 
done

dir=${options["DAY"]}/${dir_work}
mkdir -p ${dir} || errexit "Fail to create directory in EOS: $dir"
cd ${dir}       || errexit "Fail to cd to directory in EOS: $dir"

if test -f ${WORKSPACE}/controlfile -o -f controlfile ; then
  # remove previous build
  for item in build controlfile geant4 scripts setup.sh; do
    test -L $item && rm -f $item
    test -d $item && rm -rf $item
    test -f $item && rm -f $item
  done   
  test -f ${WORKSPACE}/controlfile \
    && cp -p ${WORKSPACE}/controlfile . \
    || touch controlfile
fi

ls -l ${WORKSPACE}/build
if test -d ${WORKSPACE}/build ; then 
  ln -s ${eosDir}/${srcDir} geant4
  #  the following are also in build.tgz, the extra copy helps to look at these without unpacking build.tgz
  #  location of setup.sh depends on build on hardware or docker
  test -r ${WORKSPACE}/setup.sh && cp -p ${WORKSPACE}/setup.sh .
  test -r ${WORKSPACE}/build/setup.sh && cp -p ${WORKSPACE}/build/setup.sh .
  test -r ${WORKSPACE}/rundocker.sh  &&  cp -p ${WORKSPACE}/rundocker.sh .
  (cd ${WORKSPACE} ; tar cf - scripts ) | tar xf -
  rc=$?
  echo "Copy return code for scripts/: $rc"
  items="build"
  test -r ${WORKSPACE}/setup.sh && items="$items setup.sh"
  test -r ${WORKSPACE}/scripts  && items="$items scripts"
  test -r ${WORKSPACE}/rundocker.sh  && items="$items rundocker.sh"
   
  cd ${WORKSPACE}; tar zcf $eosDir/$dir/build.tgz $items
else   
  echo "No build dir found to copy to EOS"
fi 

exit
