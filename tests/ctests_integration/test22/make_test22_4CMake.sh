#!/bin/sh
#
# Last update: 06-Mar-2020
# 
# Bash shell script to build all the executables of test22 using the
# CMake system.
#
# After the use of  EXCLUDE_FROM_ALL  on  test22/CMakeLists.txt
# (suggeseted by Ben on 26-Jul-2013 to solve a CMake problem in Windows
# platforms), the test22 executables need to be built explicitly,
# by doing "make test22_xxx". 
# This is done by this script for all the tests.
# Notice that if you don't run this script, then you get an error during
# "make install".
# 
# You need to copy this script in the cmake-build-directory/tests/ctests_integration/test22
# before you execute it as:
#   ./make_test22_4CMake.sh
#
echo "  "
echo " --- hA --- "
echo "Building test22_E941 ..."
make test22_E941
#
echo "Building test22_hA100 ..."
make test22_hA100
#
echo "Building test22_hA120 ..."
make test22_hA120
#
echo "Building test22_HARP ..."
make test22_HARP
#
echo "Building test22_NA35 ..."
make test22_NA35
#
echo "Building test22_NA49 ..."
make test22_NA49
#
echo "Building test22_NA61 ..."
make test22_NA61
#
echo "Building test22_pA200 ..."
make test22_pA200
#
echo "  "
echo " --- KmPchan --- "
echo "Building test22_KmPchan ..."
make test22_KmPchan
#
echo "  "
echo " --- KpPchan --- "
echo "Building test22_KpPchan ..."
make test22_KpPchan
#
echo "  "
echo " --- PbarA --- "
echo "Building test22_PbC_608 ..."
make test22_PbC_608
#
echo "Building test22_PbU_608 ..."
make test22_PbU_608
#
echo "Building test22_PbA_Rest ..."
make test22_PbA_Rest
#
echo "Building test22_PbDNe_Rest ..."
make test22_PbDNe_Rest
#
echo "Building test22_Nb_C_750 ..."
make test22_Nb_C_750
#
echo "Building test22_PbA_1p22_nE ..."
make test22_PbA_1p22_nE
#
echo "Building test22_PbTa_4 ..." 
make test22_PbTa_4
#
echo "  "
echo " --- PbarA_X --- "
echo "Building test22_PbarA_X ..."
make test22_PbarA_X
#
echo "  "
echo " --- PbarPchan --- "
echo "Building test22_PbarPchan ..."
make test22_PbarPchan
#
echo "Building test22_atRest ..."
make test22_atRest
#
echo "Building test22_RHO ..."
make test22_RHO
#
echo "  "
echo " --- PbarHyperonV --- "
echo "Building test22_dSigdTetPiK ..."
make test22_dSigdTetPiK
#
echo "Building test22_dSigdTet ..."
make test22_dSigdTet
#
echo "Building test22_dSigdTLLbar ..."
make test22_dSigdTLLbar
#
echo "Building test22_PbarYXPt2 ..."
make test22_PbarYXPt2
#
echo "Building test22_PbPHyperonSig ..."
make test22_PbPHyperonSig
#
echo "Building test22_PbarXe_200 ..."
make test22_PbarXe_200
#
echo "  "
echo " --- PbarPinclusive --- "
echo "Building test22_PbarP4p6 ..."
make test22_PbarP4p6
#
echo "Building test22_PbarP5p7 ..."
make test22_PbarP5p7
#
echo "Building test22_PbarP7p3 ..."
make test22_PbarP7p3
#
echo "Building test22_PbarP9p1 ..."
make test22_PbarP9p1
#
echo "Building test22_PbarP14p8 ..."
make test22_PbarP14p8
#
echo "Building test22_PbarP22p4 ..."
make test22_PbarP22p4
#
echo "Building test22_PbarP32 ..."
make test22_PbarP32
#
echo "Building test22_PbarP100 ..."
make test22_PbarP100
#
echo " --- PimPchan --- "
echo "Building test22_PimPchan ..."
make test22_PimPchan
#
echo "  "
echo " --- PiPinclusive --- "
echo "Building test22_PipP8 ..."
make test22_PipP8
#
echo "Building test22_PimP8 ..."
make test22_PimP8
#
echo "Building test22_PipP16 ..."
make test22_PipP16
#
echo "Building test22_PimP16 ..."
make test22_PimP16
#
echo "Building test22_PipP18 ..."
make test22_PipP18
#
echo "Building test22_PimP18 ..."
make test22_PimP18 
#
echo "Building test22_PipP32 ..."
make -f test22_PipP32
#
echo "Building test22_PimP40 ..."
make test22_PimP40
#
echo "Building test22_PimP58 ..."
make test22_PimP58
#
echo "Building test22_PipP100 ..."
make test22_PipP100
#
echo "Building test22_PimP100 ..."
make test22_PimP100
#
echo "Building test22_PipP175 ..."
make test22_PipP175
#
echo "Building test22_PipP250 ..."
make test22_PipP250
#
echo "Building test22_PimP360 ..."
make test22_PimP360
#
echo "  "
echo " --- PipPchan --- "
echo "Building test22_PipPchan ..."
make test22_PipPchan
#
echo "  "
echo " --- PPchan --- "
echo "Building test22_PPchan ..."
make test22_PPchan
#
echo "  "
echo " --- PPinclusive --- "
echo "Building test22_PP12 ..."
make test22_PP12
#
echo "Building test22_PP24 ..."
make test22_PP24
#
echo "Building test22_PP69 ..."
make test22_PP69
#
echo "Building test22_PP100 ..."
make test22_PP100
#
echo "Building test22_PP158 ..."
make test22_PP158
#
echo "Building test22_PP175 ..."
make test22_PP175
#
echo "Building test22_PP205 ..."
make test22_PP205
#
echo "Building test22_PP360 ..."
make test22_PP360
#
echo "Building test22_PP400 ..."
make test22_PP400
#
echo "  "
echo " --- hA_neutron --- "
echo "Building test22_Ishiba ..."
make test22_Ishiba
#
echo "Building test22_ITEP ..."
make test22_ITEP
#
echo "Building test22_ITEPx ..."
make test22_ITEPx
#
echo "Building test22_JINR ..."
make test22_JINR
#
echo "Building test22_Leray ..."
make test22_Leray
#
echo "Building test22_Nmult ..."
make test22_Nmult
#
echo "Building test22_PbA_1p22_nE ..."
make test22_PbA_1p22_nE
#
echo "  "
echo " --- NA61 --- "
#
echo "Building test22_pC ..."
make test22_pC
#
echo "Building test22_pimC ..."
make test22_pimC 
#
echo "Building test22_pimCres ..."
make test22_pimCres
#
echo "Building test22_pP ..."
make test22_pP
#
echo "Building test22_ppLam ..."
make test22_ppLam
#
