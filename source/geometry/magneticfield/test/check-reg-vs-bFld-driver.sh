
nRun="3"
epsLab="7"  #  Label for epsilon

alias testPropagateMagField.exe="$G4BIN_CMAKE/testPropagateMagField"

( cd $G4BUILD_DIR/ ; ninja testPropagateMagField testProElectroMagField 2>&1 | egrep -v 'builds involving this target' )  \
   || exit "Failed compilation"

for Stepper in 45 56 78 
do
   label="eps$epsLab.st-DP$Stepper.n$nRun"
   outfl="out-reg-drv.$label.log"
   errfl="err-reg-drv.$label.log"
   ( testPropagateMagField.exe  $Stepper  > $outfl ) 2>&1 | tee $errfl
done
#  testPropagateMagField.exe  78 > out-reg-drv.eps8.st-DP78.n1.log

for Stepper in 101 102 103 
do
   label="eps$epsLab.st-DP$Stepper.n$nRun"
#             bField-drv.eps8.st101.n1.log 
   sufx="bField-drv.$label.log"
   outfl="out-$sufx"
   errfl="err-$sufx"
   ( testPropagateMagField.exe  $Stepper > $outfl ) 2>&1 | tee $errfl
done
