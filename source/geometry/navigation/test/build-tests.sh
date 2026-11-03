#!/bin/sh
#
# Generic per-directory build+run driver for standalone G4TARGET-style
# geometry unit tests (source/geometry/*/test/test*.cc). Invoke from within
# the test source directory. Pass explicit target name(s), or "all" to build
# and run every test*.cc in this directory.
#
echo "Run on `hostname`, which is a `uname -a` machine" | tee -a test.out
echo "------------------------------------------------" | tee -a test.out
echo `date` >> test.out

original_dir=`pwd`
test_dir=`pwd | perl -pe "s#$G4SOURCE_DIR/##;"`

echo "test_dir= " $test_dir

targets=${*:-"all"}

if [[ "$targets" == "all" ]]; then
  # Expand *.cc against the source directory (original_dir), not the build
  # directory we are about to cd into - a glob evaluated after cd'ing to
  # $G4BUILD_DIR/$test_dir silently fails to expand (no .cc files live
  # there), leaving the literal string "*.cc" and corrupting every
  # subsequent path built from it.
  targets=`cd "$original_dir" && ls test*.cc 2>/dev/null`
fi

if [[ ! -x $G4BUILD_DIR ]]; then
  echo "ERROR> Must use environment variable G4BUILD_DIR for build directory."
  exit 1
else
  echo "Building tests in build directory " $G4BUILD_DIR
fi

cd $G4BUILD_DIR || exit 1

mkdir -p $test_dir || ( echo " ERROR> Cannot create directory " $test_dir " in " $G4BUILD_DIR ; exit 1 )

echo "Base location to build tests: " $G4BUILD_DIR/$test_dir
cd $test_dir  || exit 1

echo "Will build and run the following targets: "
echo "  " $targets
echo

declare -a results

for target in $targets
do
   g4target=`basename $target .cc`

   echo
   echo "target " $g4target
   echo "----------------------------------------"

   # Always delete and recreate the test's build directory before
   # configuring/building it. A reused build directory can silently keep
   # object files compiled against an older version of the Geant4 headers
   # (e.g. G4LogicalVolume.hh) than the shared libraries it links against,
   # producing an ABI/layout mismatch that crashes at runtime instead of
   # failing to build - a stale-binary bug, not a real regression, but one
   # that is otherwise easy to hit and hard to diagnose.
   if [[ -d ${g4target}.dir ]]; then
      echo " 1. Removing existing test directory ${g4target}.dir to avoid stale objects"
      rm -rf ${g4target}.dir || exit 1
   fi
   echo " 1. creating test directory "
   mkdir -p ${g4target}.dir  || exit 1

   pushd ${g4target}.dir > /dev/null || exit 1

   echo " 2. Configuring with cmake"
   cmake -DG4TARGET=$g4target $G4SOURCE_DIR/$test_dir > ../${g4target}.cmake.log 2>&1
   cmake_rc=$?

   if [[ $cmake_rc -ne 0 ]]; then
      echo " CMAKE FAILED (see ${g4target}.dir/../${g4target}.cmake.log)"
      results+=("$g4target : FAILED (cmake configure error)")
      popd > /dev/null
      continue
   fi

   echo " 3. make               ---------"
   make -j 2 > ../${g4target}.build.log 2>&1
   make_rc=$?

   if [[ $make_rc -ne 0 ]]; then
      echo " BUILD FAILED (see ${g4target}.dir/../${g4target}.build.log)"
      tail -30 ../${g4target}.build.log
      results+=("$g4target : FAILED (build error)")
      popd > /dev/null
      continue
   fi

   echo " 4. Running the test   ---------"
   ./$g4target 2>&1 | tee -a ${original_dir}/test.out
   run_rc=${PIPESTATUS[0]}

   if [[ $run_rc -eq 0 ]]; then
      results+=("$g4target : OK")
   else
      results+=("$g4target : FAILED (exit=$run_rc)")
   fi

   echo " To duplicate run: "
   echo "   " `pwd`/$g4target

   popd > /dev/null
done

echo
echo "=================== SUMMARY ==================="
for r in "${results[@]}"; do
  echo "  $r"
done
echo "================================================="
