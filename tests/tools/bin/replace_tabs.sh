#!/bin/bash
#
# Usage:
# replace_tabs.sh
#
# Script for fixing Geant4 examples coding guideline (avoid tabulations).

#set -x

CURDIR=`pwd`

for SOURCE in `ls *.cc include/*.hh src/*.cc`
do
  RESULT=`cat $SOURCE | sed 's/\t/TAB_TO_BE_REMOVED/g' | grep TAB_TO_BE_REMOVED`
  if [ ! "$RESULT" = "" ];
  then
    echo "Updating file $SOURCE ... "
    sed 's/\t/        /g' $SOURCE  > $SOURCE.1
    mv $SOURCE.1 $SOURCE
  fi
done 

cd $CURDIR
