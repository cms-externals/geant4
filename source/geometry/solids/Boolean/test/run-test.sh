
#
# Build and run a single test on Mac (works also with VecGeom) 
# 
# -- by default the test is testG4SubtractionSolid
#
#  J. Apostolakis  2021.05.04 and 05.28
#
target=${1:-"testG4SubtractionSolid"}
g4target=`basename $target .cc`

src_dir="source/geometry/solids/Boolean/test/"

if [[ -x $G4BUILD_DIR ]] ; then
   build_dir="$G4BUILD_DIR/${src_dir}/${g4target}.dir"
else
   echo "G4BUILD_DIR is not set - expect to add this test there!"
   exit
fi

if [[ -x ./other-env.sh ]] ; then
   source ./other-env.sh
fi   

echo "build_dir= " $build_dir
echo
echo "mkdir build_dir "
mkdir -p  $build_dir || exit
echo "cd build_dir "
cd $build_dir || exit

cmake $G4SOURCE_DIR/$src_dir -DG4TARGET=$g4target && make -j 4 && ./$g4target
