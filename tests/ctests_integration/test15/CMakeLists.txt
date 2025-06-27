# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(test15)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test15
  test15.cc
  src/Test15ActionInitialization.cc
  src/Test15DetectorConstruction.cc
  src/Test15EventAction.cc
  src/Test15PrimaryGeneratorAction.cc
  src/Test15PVolumeStore.cc
  src/Test15RunAction.cc
  src/Test15Run.cc
  src/Test15ShellDetectorConstruction.cc
  src/Test15StackingAction.cc
  src/Test15SteppingAction.cc)
target_include_directories(test15 PRIVATE include)
target_link_libraries(test15 PRIVATE ${Geant4_LIBRARIES})

#--Test definitions -------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test15
    COMMAND test15 -m ${CMAKE_CURRENT_SOURCE_DIR}/test15.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
