cmake_minimum_required(VERSION 3.16...3.21)
project(test37)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test37 test37.cc ${sources} ${headers})
target_include_directories(test37 PRIVATE include)
target_link_libraries(test37 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test37
    COMMAND test37 ${CMAKE_CURRENT_SOURCE_DIR}/test37.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
