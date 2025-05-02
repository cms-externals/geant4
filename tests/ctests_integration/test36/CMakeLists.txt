cmake_minimum_required(VERSION 3.16...3.21)
project(test36)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test36 test36.cc ${sources} ${headers})
target_include_directories(test36 PRIVATE include)
target_link_libraries(test36 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test36
    COMMAND test36 ${CMAKE_CURRENT_SOURCE_DIR}/test36.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
