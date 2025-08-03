#!/bin/bash
#
# This script should be called from the geant4/examples directory
# Usage:
# check_extended_examples.sh
#
# Script for checking of coding guidelines for all extended examples
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

# extended
#CATEGORY errorpropagation; polarisation
for DIR in extended; do
  #  echo ... processing $DIR
  EXTENDED_DIR=$EXAMPLES_BASEDIR/$DIR
  for CATEGORY in `ls $EXTENDED_DIR`; do
    # select directories with .README.txt 
    if [ -f ${CATEGORY}/.README.txt ]; then
      #echo ... processing $CATEGORY
      cd $CATEGORY
      CATEGORY_DIR=`pwd`
      for FILE in `find . -name .README.txt`; do
        EXAMPLE_DIR=`echo $FILE | sed sY/.README.txtYYg`
        NOT_EXAMPLE=`echo ". ./dna ./gdml ./pythia ./MPI ./HepMC ./ParticleFluence ./FlukaCern" | grep $EXAMPLE_DIR`
        if [ "${NOT_EXAMPLE}" = "" ]; then
          # echo "Go to process $EXAMPLE_DIR"
          cd $EXAMPLE_DIR
          RESULT=`$SCRIPT -silent; echo $?`
          if [ "$RESULT" == "1" ]; then
            echo "${PWD#*/extended}" >&2
            $SCRIPT -violations
          fi
          cd $CATEGORY_DIR
        fi
      done
    fi
    cd $EXTENDED_DIR
  done
done
