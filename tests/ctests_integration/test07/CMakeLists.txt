cmake_minimum_required(VERSION 3.16...3.21)
project(test07)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test07
  test07.cc
  src/T07CalorHit.cc
  src/T07CalorimeterSD.cc
  src/T07DetectorConstruction.cc
  src/T07DetectorMessenger.cc
  src/T07EventAction.cc
  src/T07EventActionMessenger.cc
  src/T07PhysicsList.cc
  src/T07PhysicsListMessenger.cc
  src/T07PrimaryGeneratorAction.cc
  src/T07PrimaryGeneratorMessenger.cc
  src/T07RunAction.cc
  src/T07SteppingAction.cc)
target_include_directories(test07 PRIVATE include)
target_link_libraries(test07 PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test07
    COMMAND test07 ${CMAKE_CURRENT_SOURCE_DIR}/test07.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  geant4_add_test(test07-largeN
    COMMAND test07 ${CMAKE_CURRENT_SOURCE_DIR}/test07.large_N.in
    DEPENDS test07
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
