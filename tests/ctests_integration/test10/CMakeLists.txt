cmake_minimum_required(VERSION 3.16...3.21)
project(test10)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test10
  test10.cc
  src/Tst10DetectorConstruction.cc
  src/Tst10DetectorMessenger.cc
  src/Tst10EventAction.cc
  src/Tst10EventActionMessenger.cc
  src/Tst10PhysicsList.cc
  src/Tst10PrimaryGeneratorAction.cc
  src/Tst10RunAction.cc)
target_include_directories(test10 PRIVATE include)
target_link_libraries(test10 PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test10
    COMMAND test10 ${CMAKE_CURRENT_SOURCE_DIR}/test10.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  geant4_add_test(test10-largeN
    COMMAND test10 ${CMAKE_CURRENT_SOURCE_DIR}/test10.large_N.in
    DEPENDS test10
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
