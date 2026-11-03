target=${1:-"testExitNormalNav"}
echo "Target = " $target

if [[ ! -f ${target}.cc ]]; then
  echo "Cannot find the target's C++ source file " ${target}.cc
  exit 1
fi

test_dir="source/geometry/volumes/test"

mkdir -p $G4BUILD_DIR/unit-tests/$test_dir || ( echo "ERROR> Failed to create build directory $G4BUILD_DIR/unit-tests/$test_dir "; exit 1 )

cd       $G4BUILD_DIR/unit-tests/$test_dir

echo "Building target " $target
if [[ ! -d $target ]]; then
   echo "Creating sub-directory for target"
   mkdir $target
   cd    $target
   cmake -G Ninja -DGeant4_DIR=${G4BUILD_DIR} -DG4TARGET=${target} ${G4SOURCE_DIR}/${test_dir}
fi
ninja -j 8

echo "To run it: "
echo      $G4BUILD_DIR/unit-tests/${test_dir}/${target}/${target}
echo "Run results:"
$G4BUILD_DIR/unit-tests/${test_dir}/${target}/${target}
