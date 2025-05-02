cmake_minimum_required(VERSION 3.16...3.21)
project(test40)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test40 test40.cc ${sources} ${headers})
target_include_directories(test40 PRIVATE include)
target_link_libraries(test40 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test40
    COMMAND test40 ${CMAKE_CURRENT_SOURCE_DIR}/test40.in 4
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test40-largeN
    COMMAND test40 ${CMAKE_CURRENT_SOURCE_DIR}/test40.large_N.in 4
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
