cmake_minimum_required(VERSION 3.16...3.21)
project(test22)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#--- Library of common code ---
file(GLOB test22_common_sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.icc)
file(GLOB test22_common_headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_library(test22_common STATIC ${test22_common_sources} ${test22_common_headers})
target_include_directories(test22_common PUBLIC include src)
target_link_libraries(test22_common PUBLIC ${Geant4_LIBRARIES})

#--- Individual applications ---
set(test22_EXE_SOURCES
  #- hA
  hA/test22_E941.cc
  hA/test22_HARP.cc
  hA/test22_NA35.cc
  hA/test22_NA49.cc
  hA/test22_NA61.cc
  hA/test22_hA100.cc
  hA/test22_hA120.cc
  hA/test22_pA200.cc
	# KmPchan
	KmPchan/test22_KmPchan.cc
  #- KpPchan
  KpPchan/test22_KpPchan.cc
  #- PbarA
  PbarA/test22_Nb_C_750.cc
  PbarA/test22_PbA_1p22_nE.cc
  PbarA/test22_PbA_Rest.cc
  PbarA/test22_PbC_608.cc
  PbarA/test22_PbDNe_Rest.cc
  PbarA/test22_PbTa_4.cc
  PbarA/test22_PbU_608.cc
  #- PbarA_X
  PbarA_X/test22_PbarA_X.cc
  #- PbarPchan
  PbarPchan/test22_PbarPchan.cc
  PbarPchan/test22_atRest.cc
  PbarPchan/test22_RHO.cc
  #- PbarHyperonV
  PbarHyperonV/test22_dSigdTetPiK.cc
  PbarHyperonV/test22_dSigdTet.cc
  PbarHyperonV/test22_dSigdTLLbar.cc
  PbarHyperonV/test22_PbarYXPt2.cc
  PbarHyperonV/test22_PbPHyperonSig.cc
  PbarHyperonV/test22_PbarXe_200.cc
  #- PbarPinclusive
  PbarPinclusive/test22_PbarP100.cc
  PbarPinclusive/test22_PbarP14p8.cc
  PbarPinclusive/test22_PbarP22p4.cc
  PbarPinclusive/test22_PbarP32.cc
  PbarPinclusive/test22_PbarP4p6.cc
  PbarPinclusive/test22_PbarP5p7.cc
  PbarPinclusive/test22_PbarP7p3.cc
  PbarPinclusive/test22_PbarP9p1.cc
  #- PimPchan
  PimPchan/test22_PimPchan.cc
  #- PiPinclusive
  PiPinclusive/test22_PimP100.cc
  PiPinclusive/test22_PimP16.cc
  PiPinclusive/test22_PimP40.cc
  PiPinclusive/test22_PipP100.cc
  PiPinclusive/test22_PipP16.cc
  PiPinclusive/test22_PipP175.cc
  PiPinclusive/test22_PipP8.cc
  PiPinclusive/test22_PimP8.cc
  PiPinclusive/test22_PimP18.cc
  PiPinclusive/test22_PipP18.cc
  PiPinclusive/test22_PipP32.cc
  PiPinclusive/test22_PimP58.cc
  PiPinclusive/test22_PipP250.cc
  PiPinclusive/test22_PimP360.cc
  #- PipPchan
  PipPchan/test22_PipPchan.cc
  #- PPchan
  PPchan/test22_PPchan.cc
  #- PPinclusive
  PPinclusive/test22_PP100.cc
  PPinclusive/test22_PP12.cc
  PPinclusive/test22_PP158.cc
  PPinclusive/test22_PP175.cc
  PPinclusive/test22_PP205.cc
  PPinclusive/test22_PP24.cc
  PPinclusive/test22_PP360.cc
  PPinclusive/test22_PP400.cc
  PPinclusive/test22_PP69.cc
  #- hA_neutron
  hA_neutron/test22_Ishiba.cc
  hA_neutron/test22_ITEP.cc
  hA_neutron/test22_ITEPx.cc
  hA_neutron/test22_JINR.cc
  hA_neutron/test22_Leray.cc
  hA_neutron/test22_Nmult.cc
  hA_neutron/test22_PbA_1p22_nE_bis.cc
  #- NA61
  NA61/test22_pC.cc
  NA61/test22_pimC.cc
  NA61/test22_pimCres.cc
  NA61/test22_pP.cc
  NA61/test22_ppLam.cc
 	)

foreach(_test22_src ${test22_EXE_SOURCES})
  get_filename_component(_test22_exe ${_test22_src} NAME_WE)
	add_executable(${_test22_exe} ${_test22_src})
	target_link_libraries(${_test22_exe} PRIVATE test22_common)
endforeach()

# - Testing applications
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test22_NA49-exerciser
    COMMAND test22_NA49 ${CMAKE_CURRENT_SOURCE_DIR}/hA/NA49pC.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
