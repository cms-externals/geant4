cmake_minimum_required(VERSION 3.16...3.21)
project(test65)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test65 test65.cc ${sources} ${headers})
target_include_directories(test65 PRIVATE include)
target_link_libraries(test65 PRIVATE ${Geant4_LIBRARIES})

add_executable(test65_hadronic_exerciser test65_hadronic_exerciser.cc)
target_link_libraries(test65_hadronic_exerciser PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING AND NOT WIN32)
  #-----------------------------------------------------------------------
  # Test definitions
  #
  list(APPEND GEANT4_TEST_ENVIRONMENT "G4LENDDATA=${PROJECT_BINARY_DIR}/Data/g4lend")

  geant4_add_test(test65-download-data
    COMMAND ${CMAKE_COMMAND} -DALTERNATE_URL="${GEANT4_TEST_REFERENCES_URL}/tests/test65" -P ${PROJECT_SOURCE_DIR}/DownloadLENDData.cmake
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    TIMEOUT 3600)

  geant4_add_test(test65_hadronic_exerciser
    COMMAND test65_hadronic_exerciser
    OUTPUT  test65.in
    DEPENDS test65-download-data)

  geant4_add_test(test65
    COMMAND test65 test65.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test65_hadronic_exerciser test65-download-data)
endif()
