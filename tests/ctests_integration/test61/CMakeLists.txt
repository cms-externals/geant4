cmake_minimum_required(VERSION 3.16...3.21)
project(test61)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test61 test61.cc ${sources} ${headers})
target_include_directories(test61 PRIVATE include)
target_link_libraries(test61 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  # QMD tests
  geant4_add_test(test61
    COMMAND test61 ${CMAKE_CURRENT_SOURCE_DIR}/test61.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  # LightIonQMD tests
  geant4_add_test(test61-liqmd
    COMMAND test61 ${CMAKE_CURRENT_SOURCE_DIR}/test61.in -liqmd
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    LABELS Nightly Continuous)
endif()
