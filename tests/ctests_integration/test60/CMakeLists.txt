cmake_minimum_required(VERSION 3.16...3.21)
project(test60)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test60 test60.cc ${sources} ${headers})
target_include_directories(test60 PRIVATE include)
target_link_libraries(test60 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  foreach(_i 2 4 6 8)
    geant4_add_test(test60-${_i}
    COMMAND test60 -m ${CMAKE_CURRENT_SOURCE_DIR}/test60.in -p ${_i}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()

  geant4_add_test(test60-largeN
    COMMAND test60 -m ${CMAKE_CURRENT_SOURCE_DIR}/test60.large_N.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
