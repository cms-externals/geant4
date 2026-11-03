epsilon="0.001"
beta="0.1"

n=${1:-"001"}

#unType="v2-useMid"
runType="v4-accurate"

echo "Run parameters (fixed) : epsilon= $epsilon beta= $beta"

for fr in  1.0 0.5 0.25 0.125 0.0625 0.03125 
do
   outfile="out-b${beta}-fr${fr}-${runType}.${n}.log"
   echo "run with fraction= $fr"
   $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 $epsilon 0.1 $fr >  $outfile
   chmod -w $outfile
done
