cmake_minimum_required(VERSION 3.16...3.21)
project(test17)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test17
  test17.cc
  src/PhysicsList.cc
  src/PhysicsListMessenger.cc
  src/StepMax.cc
  src/StepMaxMessenger.cc
  src/Test17DetectorConstruction.cc
  src/Test17DetectorMessenger.cc
  src/Test17EventAction.cc
  src/Test17PrimaryGeneratorAction.cc
  src/Test17PrimaryGeneratorMessenger.cc
  src/Test17RunAction.cc
  src/Test17SteppingAction.cc)
target_include_directories(test17 PRIVATE include)
target_link_libraries(test17 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  #---Test definitions----------------------------------------------------------------------------------
  geant4_add_test(test17
    COMMAND test17 ${CMAKE_CURRENT_SOURCE_DIR}/test17.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  geant4_add_test(test17-largeN
    COMMAND test17 ${CMAKE_CURRENT_SOURCE_DIR}/test17.large_N.in
    DEPENDS test17
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
