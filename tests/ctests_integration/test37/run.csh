#/bin/csh

cd $VFEM/test37
mkdir -p $REFERENCE
cd $REFERENCE

set    work = "$G4BWORK/test37"
set    dir  = "$G4INSTALL/tests/ctests_integration/test37/"

rm -rf *.gz *.out *.log

setenv PHYSLIST    emstandard
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandard_opt3
set    phys = "opt3"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandard_opt4
set    phys = "opt4"
source ${dir}run_single.csh ${phys} ${work} ${dir}

setenv PHYSLIST    emstandardWVI
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir}
