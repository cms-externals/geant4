#!/bin/bash -e
# Script to run full Geant4 builds inside docker containers
# Beware the code here highly relies on multiple variable defined inside the
# jenkins job configuration
# Docker command
# -w directory                                     Working directory inside the container
# -e variable=value                                Environment variable to define inside the container
# -u username                                      User to run all build instructions
# --name $NAME                                     Name to identify the container
# --hostname $HOSTNAME-docker                      Name to identify the docker host (used for cdash)
# --cpus=$DOCKER_CPUS                              CPU container restriction
# -v host_path:container_path                      Folder to bind from docker host to docker container
# gitlab-registry.cern.ch/sft/docker:<image_name>  Image and tag to run in the docker container
# <command>                                        Command to execute inside the container (here: script with build instructions followed by params)

export WORKSPACE_HOST=$WORKSPACE
export WORKSPACE='/build/jenkins/workspace'
touch $WORKSPACE_HOST/controlfile

errexit(){
  echo "Error: $(basename $0): $1"
  exit 1
}

DOCKER_IMAGE=''

#If OS is set, use this over LABEL
if test -n "${OS}" ; then
  DOCKER_IMAGE="sft/docker/${OS}"
  export LABEL=${OS}

  # special cases for non default image to use
  case $OS in
    centos8)  
      DOCKER_IMAGE="sft/docker/${OS}-g4"
      ;;
  esac
fi

# Choose the correct docker image, see gitlab.cern.ch/sft/docker/
export OSname="unix"

if test -z "${DOCKER_IMAGE}" ; then
  case "$LABEL" in
    lcg_docker_cc7)
      DOCKER_IMAGE="sft/docker:lcg-cc7"
      ;;
    lcg_docker_slc6)
      DOCKER_IMAGE="sft/docker/slc6"
      ;;
    lcg_docker_cc7_noafs)
      DOCKER_IMAGE="sft/docker/centos7"
      export LABEL="centos7"
      ;;
    lcg_docker_c8_noafs)
      DOCKER_IMAGE="sft/docker/centos8-g4"
      export LABEL="centos8"
      ;;
    docker-thin)
      DOCKER_IMAGE="sft/docker:cc7-thin"
      ;;
    *)
      errexit "Docker image for $LABEL not configured"
      ;;
  esac
fi

# Extract name from build tag
export NAME=`echo $BUILD_TAG | tr "," "-" | tr "=" "-"`
# Prepare number of cpus to use in the container
TOTALCPU=`nproc --all`
CONTAINERS_LIMIT=$(($EXECUTOR_NUMBER+1))
if [ -z $DOCKER_CPUS ]; then
  DOCKER_CPUS=$(($TOTALCPU/$CONTAINERS_LIMIT))
fi

# Pull image from gitlab
docker pull gitlab-registry.cern.ch/$DOCKER_IMAGE

# Find out script to run inside the container, default to jk-all, then add script directory to each part
script_dir=/scripts/tests/tools/ctest

test -z "$build_script" && build_script=jk-all || true
echo ${build_script} | sed -e's/\&\&/\&\&#/g' -e's/||/||#/g' -e's/;/;#/g'
IFS='#' read -a scripts <<<$(echo ${build_script} | sed -e's/\&\&/\&\&#/g' -e's/||/||#/g' -e's/;/;#/g')

echo ${scripts[0]} == ${scripts[1]} == ${scripts[2]} ==

build_script=""
# don't use scripts[@], as this will (re-)split items at spaces
for ((i=0; i< ${#scripts[*]};i++)) ; do
  build_script="${build_script} ${script_dir}/$(sed 's/^[[:space:]]*//' <<< "$var"${scripts[$i]})"
done
echo "build scripts=$build_script."

# parse extra_env for variable to pass into docker...
extra_env=""
for x in $docker_env; do
  extra_env="$extra_env -e $x=${!x}"
done
echo "extra env: $extra_env"

#Preparation for test
if [ -e $WORKSPACE_HOST/docker ]; then
  rm -rf $WORKSPACE_HOST/docker
fi
mkdir -p $WORKSPACE_HOST/docker

# make sure needed cvmfs volumes are mounted 
stat /cvmfs/geant4.cern.ch/geant4
stat /cvmfs/sft.cern.ch/lcg

cat > rundocker.sh << EOF
#!/bin/bash
# rename build back to docker, undo rename below
test -d build && mv build docker
docker run -it --rm --network=host \
  -e WORKSPACE=$WORKSPACE \
  -e WORKSPACE_HOST=$WORKSPACE_HOST \
  -e JOB_NAME=$JOB_NAME \
  -e VERSION=$VERSION \
  -e MODE=$MODE \
  -e BUILDTYPE=$BUILDTYPE \
  -e EXTERNALS=$EXTERNALS \
  -e BuildOptions="$BuildOptions" \
  -e ExtraCMakeOptions="$ExtraCMakeOptions" \
  -e COMPILER=$COMPILER \
  -e LABEL=$LABEL \
  -e THREAD=$THREAD \
  -e GIT_BRANCH=$GIT_BRANCH \
  -e gitlabUserName="$gitlabUserName" \
  -e MergeRequestId=$MergeRequestId \
  -e MergeRequestLastCommit=$MergeRequestLastCommit \
  $extra_env \
  -e BUILD_URL=${BUILD_URL} \
  -e OSname=$OSname \
  -u sftnight \
  --name $NAME \
  --hostname $HOSTNAME-docker \
  --cpus=$DOCKER_CPUS \
  --security-opt seccomp:unconfined \
  -v /ccache:/ccache \
  -v /var/run/nscd:/var/run/nscd \
  -v /tmp:/tmp \
  -v /cvmfs:/cvmfs \
  -v $PWD/geant4:/geant4 \
  -v $PWD/scripts:/scripts \
  --mount type=bind,source="$WORKSPACE_HOST/docker",target=/build/jenkins/workspace \
  gitlab-registry.cern.ch/$DOCKER_IMAGE \
  bash -c "cd $WORKSPACE && source setup.sh && bash" 
EOF
chmod u+x rundocker.sh

# Run container and full build inside
docker run --rm --network=host \
  -e WORKSPACE=$WORKSPACE \
  -e WORKSPACE_HOST=$WORKSPACE_HOST \
  -e JOB_NAME=$JOB_NAME \
  -e VERSION=$VERSION \
  -e MODE=$MODE \
  -e BUILDTYPE=$BUILDTYPE \
  -e EXTERNALS=$EXTERNALS \
  -e BuildOptions="$BuildOptions" \
  -e ExtraCMakeOptions="$ExtraCMakeOptions" \
  -e COMPILER=$COMPILER \
  -e LABEL=$LABEL \
  -e THREAD=$THREAD \
  -e GIT_BRANCH=$GIT_BRANCH \
  -e gitlabUserName="$gitlabUserName" \
  -e MergeRequestId=$MergeRequestId \
  -e MergeRequestLastCommit=$MergeRequestLastCommit \
  $extra_env \
  -e BUILD_URL=${BUILD_URL} \
  -e OSname=$OSname \
  -u sftnight \
  --name $NAME \
  --hostname $HOSTNAME-docker \
  --cpus=$DOCKER_CPUS \
  --security-opt seccomp:unconfined \
  -v /ccache:/ccache \
  -v /var/run/nscd:/var/run/nscd \
  -v /tmp:/tmp \
  -v /cvmfs:/cvmfs \
  -v $PWD/geant4:/geant4 \
  -v $PWD/scripts:/scripts \
  --mount type=bind,source="$WORKSPACE_HOST/docker",target=/build/jenkins/workspace \
  gitlab-registry.cern.ch/$DOCKER_IMAGE \
  bash -c "df -h; ls -l /geant4; cd $WORKSPACE; touch $WORKSPACE/is_docker; $build_script; echo $?"
echo "Docker return code $?"

# If there is a install.sh, run it
if test -n "${eos_install_tmp}" ; then
  cd $WORKSPACE_HOST/docker
  export WORKSPACE=$WORKSPACE_HOST
  source install.sh  
else
  # prepare to copy build area to eos - rename docker to build
  mv docker build
fi
