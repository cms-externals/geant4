cmake_minimum_required(VERSION 3.16...3.21)
project(test05)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test05
  test05.cc
  src/Tst05DetectorConstruction.cc
  src/Tst05PhysicsList.cc
  src/Tst05PrimaryGeneratorAction.cc
  src/Tst05RunAction.cc)
target_include_directories(test05 PRIVATE include)
target_link_libraries(test05 PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test05
    COMMAND test05 ${CMAKE_CURRENT_SOURCE_DIR}/test05.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  geant4_add_test(test05-largeN
    COMMAND test05 ${CMAKE_CURRENT_SOURCE_DIR}/test05.largeN.in
    DEPENDS test05
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
