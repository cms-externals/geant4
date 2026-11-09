#!/bin/sh
mkdir -p build || ( echo "ERROR> Failed to create build directory build "; exit 1 )
cd build
echo "Run on `hostname`, which is a `uname -a` machine" | tee -a test.out
echo `date` >> test.out
cmake -DCMAKE_PREFIX_PATH=$GEANT4_BUILD_DIR ..
make -j8
for i in ../*.cc
do
  j=`basename $i .cc`
  echo Test output for $j... >> test.out
  ./$j >> test.out 2>&1;
done

