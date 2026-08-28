cmake_minimum_required(VERSION 3.16...3.21)
project(test13)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test13
  test13.cc
  src/Tst13ActionInitialization.cc
  src/Tst13DetectorConstruction.cc
  src/Tst13DetectorMessenger.cc
  src/Tst13PhysicsList.cc
  src/Tst13PrimaryGeneratorAction.cc
  src/Tst13RunAction.cc
  src/Tst13SteppingAction.cc)
target_include_directories(test13 PRIVATE include)
target_link_libraries(test13 PRIVATE ${Geant4_LIBRARIES})

add_executable(test13_hadronic_exerciser test13_hadronic_exerciser.cc)
target_link_libraries(test13_hadronic_exerciser PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test13_hadronic_exerciser
    COMMAND test13_hadronic_exerciser
    OUTPUT  test13.in)
  geant4_add_test(test13
    COMMAND test13 test13.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test13_hadronic_exerciser)
endif()
