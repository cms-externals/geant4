cmake_minimum_required(VERSION 3.16...3.21)
project(test09)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test09
  test09.cc
  src/Tst09DetectorConstruction.cc
  src/Tst09DetectorMessenger.cc
  src/Tst09PhysicsList.cc
  src/Tst09PrimaryGeneratorAction.cc
  src/Tst09RunAction.cc
  src/Tst09SteppingAction.cc
  src/Tst09TrackingAction.cc)
target_include_directories(test09 PRIVATE include)
target_link_libraries(test09 PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test09
    COMMAND test09 ${CMAKE_CURRENT_SOURCE_DIR}/test09.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  geant4_add_test(test09-largeN
    COMMAND test09 ${CMAKE_CURRENT_SOURCE_DIR}/test09.large_N.in
    DEPENDS test09
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
