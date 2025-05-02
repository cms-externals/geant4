#!/bin/bash

# on Mac try to mount missing cvmfs volumes.
# inspired from lcgcmake/.../wait_for_cvmfs.sh

mount() {

        vol=$1.cern.ch
        current_rev=`xattr -p user.revision /cvmfs/$vol`
        if [ $? -ne 0 ]; then
            # Probably the CVMFS volume is not mounted. Try to mount it...
            echo "Mounting CVMFS volume $vol"
            sudo mount -t cvmfs $vol /Users/Shared/cvmfs/$vol
        fi
}

# only check on Mac
test  $(uname -s) == Darwin  && true || exit 0

# Check there is /cvmfs, ~/../Shared, and it has cvmfs, if not exit with Success
test -d /cvmfs -o -L /cvmfs     && true || exit 0
test -d ${HOME}/../shared/cvmfs && true || exit 0

# we need to mount cvmfs-config before we can mount the other repositories

if mount cvmfs-config ; then   # if this fails: don't try other volumes.
   mount sft
   mount geant4
else
   echo "failed to mount cvmfs-config: giving up"
fi
exit 0
