cmake_minimum_required(VERSION 3.16...3.21)
project(test08)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test08
  test08.cc
  src/ActionInitialization.cc
  src/DetectorConstruction.cc
  src/EventAction.cc
  src/PrimaryGeneratorAction.cc
  src/ProcCounterAccumulable.cc
  src/RunAction.cc
  src/SteppingAction.cc
  src/TrackingAction.cc)
target_include_directories(test08 PRIVATE include)
target_link_libraries(test08 PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test08
    COMMAND test08 ${CMAKE_CURRENT_SOURCE_DIR}/test08.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
