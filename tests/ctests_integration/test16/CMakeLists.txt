cmake_minimum_required(VERSION 3.16...3.21)
project(test16)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test16
  test16.cc
  src/Tst16ActionInitialization.cc
  src/Tst16DetectorConstruction.cc
  src/Tst16DetectorMessenger.cc
  src/Tst16PhysicsList.cc
  src/Tst16PrimaryGeneratorAction.cc
  src/Tst16RunAction.cc
  src/Tst16SteppingAction.cc)
target_include_directories(test16 PRIVATE include)
target_link_libraries(test16 PRIVATE ${Geant4_LIBRARIES})

add_executable(test16_hadronic_exerciser test16_hadronic_exerciser.cc)
target_link_libraries(test16_hadronic_exerciser PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test16_hadronic_exerciser
    COMMAND test16_hadronic_exerciser
    OUTPUT  test16.in)
  geant4_add_test(test16
    COMMAND test16 test16.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test16_hadronic_exerciser)
endif()
