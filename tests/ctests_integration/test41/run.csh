#/bin/csh

cd $VFEM/test41
mkdir -p $REFERENCE
cd $REFERENCE

rm -f *.*

set    work = "$G4BWORK/test41"
set    dir  = "$G4INSTALL/tests/ctests_integration/test41/"

setenv PHYSLIST    emstandard_opt0
set    phys = "opt0"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res0.out

setenv PHYSLIST    emstandardUB
set    phys = "optUB"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resUB.out

setenv PHYSLIST    emstandard_optG
set    phys = "optG"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resWVI.out

setenv PHYSLIST    emstandardSS
set    phys = "optS"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& resSS.out

setenv PHYSLIST    emstandard_option4
set    phys = "opt4"
source ${dir}run_single.csh ${phys} ${work} ${dir} >& res4.out

