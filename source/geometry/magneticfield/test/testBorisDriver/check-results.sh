   ( cd $G4BUILD_DIR/ ; ninja -j 8 && ( cd source/geometry/magneticfield/test/field08 ; ninja -j 8 ) )  || exit 
   set fr="0.5" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr   >  out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 %  4 == 0 &&0 <= $1 && $1 < 1000 ) '  out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log

   set fr="0.25"  ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 %  8 == 0 &&0 <= $1 && $1 < 1000 ) '  out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log

   set fr="0.125" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 %  8 == 0 &&0 <= $1 && $1 < 1000 ) '  out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 % 16 == 0 &&0 <= $1 && $1 < 1000 ) '  out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log

   set fr="0.0625" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 % 32 == 0 &&0 <= $1 && $1 < 1000 ) ' out-b0.1-fr0.0625-v4c-AccurateAdv.${runId}.log

   set fr="0.03125" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 % 64 == 0 &&0 <= $1 && $1 < 1000 ) ' out-b0.1-fr0.03125-v4c-AccurateAdv.${runId}.log

   set fr="0.015625" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 % 128 == 0 &&0 <= $1 && $1 < 10000 ) ' out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log

   set fr="0.0078125" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
   awk ' /^Step#/ || ( $1 % 256 == 0 &&0 <= $1 && $1 < 10000 ) ' out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log

   set fr="0.00390625" ; $G4BUILD_DIR/source/geometry/magneticfield/test/field08/field08 0.001 0.1 $fr > out-b0.1-fr${fr}-v4c-AccurateAdv.${runId}.log
