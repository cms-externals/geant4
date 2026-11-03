#!/bin/sh
cd build || ( echo "ERROR> No build directory present "; exit 1 )
echo "Run on `hostname`, which is a `uname -a` machine"
echo `date`
for i in ../*.cc
do
  j=`basename $i .cc`
  echo Test output for $j...
  ./$j
done

