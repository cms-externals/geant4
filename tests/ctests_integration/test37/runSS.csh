#/bin/csh

cd $VFEM/test37
mkdir -p $REFERENCE
cd $REFERENCE

set    work = "$G4BWORK/test37"
set    dir  = "$G4INSTALL/tests/ctests_integration/test37/"

setenv PHYSLIST    emstandardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir}

