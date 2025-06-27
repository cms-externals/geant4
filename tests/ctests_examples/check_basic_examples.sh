#!/bin/bash
#
# This script should be called from the geant4/examples directory
# Usage:
# check_basic_examples.sh
#
# Script for checking of coding guidelines for all basic examples
# (by processing check_example.sh)
#
# By I. Hrivnacova, IJCLab Orsay

#set -x
this_script="$0"
calldir=`pwd`
if [ ${this_script##/} = ${this_script##~} ] ; then
  this_script="$calldir/$this_script"
fi
SCRIPT="$(dirname $this_script)/check_example.sh"
EXAMPLES_BASEDIR="$(dirname $this_script)/../../examples"

for DIR in basic; do
  BASIC_DIR=$EXAMPLES_BASEDIR/$DIR
  for EXAMPLE in B1 B2/B2a B2/B2b B3/B3a B3/B3b B4/B4a B4/B4b B4/B4c B4/B4d B5; do
    cd $BASIC_DIR/$EXAMPLE
    RESULT=`$SCRIPT -silent; echo $?`
    if [ "$RESULT" == "1" ]; then
      echo "${PWD#*/extended}" >&2
      $SCRIPT -violations
    fi
    cd $BASIC_DIR
  done
done
