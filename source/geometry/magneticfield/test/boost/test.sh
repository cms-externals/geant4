#!/bin/sh
#
# A simple script to run all the tests in this directory and check
# their results against the expected (previous) results
#
# $Name: geant4-09-02-ref-05 $
#

# Find the Geant4 environment
if [[ "x" != "x$G4BUILD_DIR" ]]; then
   . $G4BUILD_DIR/geant4make.sh
fi

echo "Running on `hostname`, which is a `uname -a` machine" 
host=`hostname`

#  Choice of build engine
# MAKE=make
MAKE=ninja
#####

#  Choice of build target directory BUILD_DIR
BUILD_DIR=$G4BUILD_DIR
#G4BIN=$G4WORKDIR/bin/
#G4BIN=$G4BIN
#G4BIN=$G4BIN_GMAKE/bin/
G4BIN=$G4BIN_CMAKE/

#   Location of executable binaries  BINDIR
#BINDIR=$G4BIN/$G4SYSTEM/
BINDIR=$G4BIN/

target=testDoPri5
echo  "Compiling $target ... "
( ( cd $BUILD_DIR ; $MAKE $target || exit ) ) 2>&1 | grep -v 'multiple rules generate'
if [[ ! -x $BINDIR/$target ]] ; then
  echo "Could not find executable $target in directory $BINDIR" 
  exit 1
fi

echo  "Executing $target ..."

$BINDIR/$target $1 

