cmake_minimum_required(VERSION 3.16...3.21)
project(test12)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test12
  test12.cc
  src/Tst12ActionInitialization.cc
  src/Tst12DetectorConstruction.cc
  src/Tst12DetectorMessenger.cc
  src/Tst12PrimaryGeneratorAction.cc
  src/Tst12RunAction.cc
  src/Tst12SteppingAction.cc)
target_include_directories(test12 PRIVATE include)
target_link_libraries(test12 PRIVATE ${Geant4_LIBRARIES})

add_executable(test12_hadronic_exerciser test12_hadronic_exerciser.cc)
target_link_libraries(test12_hadronic_exerciser PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test12_hadronic_exerciser
    COMMAND test12_hadronic_exerciser
    OUTPUT  test12.in)
  geant4_add_test(test12
    COMMAND test12 test12.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test12_hadronic_exerciser)
  geant4_add_test(test12-FTFBIC
    COMMAND test12 test12.in FTF_BIC
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test12
    TIMEOUT 2000)
  geant4_add_test(test12-QGSBIC
    COMMAND test12 test12.in QGS_BIC
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test12)
  geant4_add_test(test12-QGSP_FTFP_BERT
    COMMAND test12 test12.in QGSP_FTFP_BERT
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test12)
endif()
