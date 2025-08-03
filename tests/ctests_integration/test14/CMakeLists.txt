# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(test14)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test14
  test14.cc
  src/PhysicsList.cc
  src/PhysicsListMessenger.cc
  src/StepMax.cc
  src/StepMaxMessenger.cc
  src/Tst14ActionInitializer.cc
  src/Tst14CalorHit.cc
  src/Tst14CalorimeterSD.cc
  src/Tst14DetectorConstruction.cc
  src/Tst14DetectorMessenger.cc
  src/Tst14PrimaryGeneratorAction.cc
  src/Tst14ProcCallSA.cc
  src/Tst14RunAction.cc
  src/Tst14SteppingAction.cc
  src/Tst14TrackingAction.cc)
target_include_directories(test14 PRIVATE include)
target_link_libraries(test14 PRIVATE ${Geant4_LIBRARIES})

set(INPUTS test14_lowE.large_N.in test14_penelope.large_N.in test14_polarised.large_N.in)

foreach(_input ${INPUTS})
  configure_file(
	  ${PROJECT_SOURCE_DIR}/${_input}
	  ${PROJECT_BINARY_DIR}/${_input}
	  COPYONLY)
endforeach()

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test14-lowE
    COMMAND test14 test14_lowE.large_N.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test14-penelope
    COMMAND test14 test14_penelope.large_N.in
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test14-polarised
    COMMAND test14 test14_polarised.large_N.in
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

