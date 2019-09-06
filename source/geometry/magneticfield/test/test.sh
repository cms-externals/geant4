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

sampleLines=8 ##  How many lines of same 'error' to display

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

OUTPUT_DIR=Outputs

target=testInterpolation
echo  "Compiling $target ... "
( ( cd $BUILD_DIR ; $MAKE $target )  || exit ) 2>&1 | grep -v 'multiple rules generate'
if [[ ! -x $BINDIR/$target ]] ; then
  echo "Could not find executable $target in directory $BINDIR" 
  exit 1
fi
echo  "Executing $target ..."
$BINDIR/$target $1 

exit 0

echo ""; echo "##############################################################################"
target=testProElectroMagField
echo  "Compiling $target ... "
## ( cd $BUILD_DIR ; $MAKE $target )   || exit
( ( cd $BUILD_DIR ; $MAKE $target )   || exit ) 2>&1 | grep -v 'multiple rules generate'

echo  "Executing $target ..."
for n in 1 2 3 4 8   23 45 56 78 145 745
do

  NEWOUT_FILE=$OUTPUT_DIR/$target.newout$n
  NEWERR_FILE=$OUTPUT_DIR/$target.newerr$n
  OUT_FILE=$OUTPUT_DIR/$target.out$n
  ERR_FILE=$OUTPUT_DIR/$target.err$n

  echo "Executing with stepper choice $n ..."
  $BINDIR/$target $n  > $NEWOUT_FILE 2> $NEWERR_FILE

  #---- Old code to check stdout & stderr
  # echo  ".. difference from expected output: "
  # diff -wb $OUT_FILE $NEWOUT_FILE
  # sleep 1;
  # echo  ".. difference from expected error: "
  # diff -wb $ERR_FILE $NEWERR_FILE
  # sleep 1;
  # echo  " "

  if [[ -f $OUT_FILE ]]; then
     if [[ `cmp --silent $OUT_FILE $NEWOUT_FILE` ]]; then
	 echo  ".. difference(s) from expected output: "
	 diff -wb $OUT_FILE $NEWOUT_FILE
	 sleep 2;
     fi
  else
     echo  " Expected output *not* found. Time to create " $OUT_FILE
  fi
  # if [[ -f $ERR_FILE ]]; then
  if [[ -f $NEWERR_FILE ]]; then
      if [[ `wc -l $NEWERR_FILE | awk ' { print $1; } '` != 0 ]]; then
	  if [[ -f $ERR_FILE ]]; then
              if [[ `cmp --silent $ERR_FILE $NEWERR_FILE` ]]; then	      
     	         echo "Differences seen from expected error output: "
		 diff -wb $ERR_FILE $NEWERR_FILE
	      else
		 echo "Same error output."
	      fi
              echo "Sample of (current) error output: "
	      head -$sampleLines  $NEWERR_FILE	      
	  else
   	      echo  "Unexpected error output: (max 50 lines - top)"
	      head -50 $NEWERR_FILE
	  fi
	  sleep 5;
      else
	  ## echo " Deleting empty error file ... "
	  rm $NEWERR_FILE
      fi
  else
      echo "WARNING> 'Error' output file $NEWERR_FILE not found."
  fi      
  # echo  " "
done

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo ""; echo "##############################################################################"    
  echo  "Compiling $target ... "
  ## ( cd $BUILD_DIR ; $MAKE $target )  || exit
  ( ( cd $BUILD_DIR ; $MAKE $target )   || exit ) 2>&1 | grep -v 'multiple rules generate'  
  echo  "Executing $target ..."
  for n in  4 3 2 1 0  # 23 45 56 78 145 745
  do
  
    NEWOUT_FILE=$OUTPUT_DIR/$target.newout$n
    NEWERR_FILE=$OUTPUT_DIR/$target.newerr$n
    OUT_FILE=$OUTPUT_DIR/$target.out$n
    ERR_FILE=$OUTPUT_DIR/$target.err$n
  
    echo "Executing with stepper choice $n ..."
    $BINDIR/$target $n  > $NEWOUT_FILE 2> $NEWERR_FILE
    echo  ".. difference from expected output: "
    diff -wb $OUT_FILE $NEWOUT_FILE
    sleep 1;
    echo  ".. difference from expected error: "
    diff -wb $ERR_FILE $NEWERR_FILE
    sleep 1;
    # echo  " "
  done
done

exit

for i in *Spin.cc
do
  target=`basename $i .cc`
  echo  "Compiling $target ... "
  ## gmake G4TARGET=$target || exit
  ( ( cd $BUILD_DIR ; $MAKE $target )   || exit ) 2>&1 | grep -v 'multiple rules generate'    
  echo  "Executing $target ..."
  
  NEWOUT_FILE=$OUTPUT_DIR/$target.newout$n
  NEWERR_FILE=$OUTPUT_DIR/$target.newerr$n
  OUT_FILE=$OUTPUT_DIR/$target.out$n
  ERR_FILE=$OUTPUT_DIR/$target.err$n
  
  $BINDIR/$target $n  > $NEWOUT_FILE 2> $NEWERR_FILE
  echo  ".. difference from expected output: "
  diff -wb $OUT_FILE NEWOUT_FILE
  sleep 1;
  echo  ".. difference from expected error: "
  diff -wb $ERR_FILE $NEWERR_FILE
  sleep 1;
  echo  " "
done

exit
