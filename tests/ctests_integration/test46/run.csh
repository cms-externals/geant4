#/bin/csh

cd $VFEM/test46
mkdir -p $REFERENCE

cd $REFERENCE

setenv PRIMARYBEAM e-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

setenv PHYSLIST FTFP_BERT_EMM

source $G4INSTALL/tests/ctests_integration/test46/run_em.csh e50gev

cd $REFERENCE
setenv PRIMARYBEAM pi-
mkdir -p $PRIMARYBEAM
cd $PRIMARYBEAM

source $G4INSTALL/tests/ctests_integration/test46/run_em.csh pi-5gev
source $G4INSTALL/tests/ctests_integration/test46/run_em.csh pi-9gev

#
