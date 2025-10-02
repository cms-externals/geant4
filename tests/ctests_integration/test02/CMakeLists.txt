cmake_minimum_required(VERSION 3.16...3.21)
project(test02)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test02
  test02.cc
  src/Tst02DetectorConstruction.cc
  src/Tst02RunAction.cc
  src/Tst02DetectorMessenger.cc
  src/Tst02SteppingAction.cc
  src/Tst02PrimaryGeneratorAction.cc)

target_include_directories(test02 PRIVATE include)
target_link_libraries(test02 PRIVATE ${Geant4_LIBRARIES})

add_executable(test02_had test02_hadronic_exerciser.cc)
target_link_libraries(test02_had PRIVATE ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test02
    COMMAND test02 ${CMAKE_CURRENT_SOURCE_DIR}/test02.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test02-largeN
    COMMAND test02 ${CMAKE_CURRENT_SOURCE_DIR}/test02.large_N.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test02-gen
    COMMAND test02_had
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    OUTPUT test02-had.in)

  geant4_add_test(test02-had
    COMMAND test02 test02-had.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test02-gen)
endif()
