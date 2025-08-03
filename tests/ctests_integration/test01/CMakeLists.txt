cmake_minimum_required(VERSION 3.16...3.21)
project(test01)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test01 test01.cc
  src/Tst01DetectorConstruction.cc
  src/Tst01PrimaryGeneratorMessenger.cc
  src/Tst01DetectorMessenger.cc
  src/Tst01RunAction.cc
  src/Tst01PhysicsList.cc
  src/Tst01SteppingAction.cc
  src/Tst01PrimaryGeneratorAction.cc
  )
target_include_directories(test01 PRIVATE include)
target_link_libraries(test01 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test01
    COMMAND test01 ${CMAKE_CURRENT_SOURCE_DIR}/test01.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test01-largeN
    COMMAND test01 ${CMAKE_CURRENT_SOURCE_DIR}/test01.large_N.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test01)
endif()
