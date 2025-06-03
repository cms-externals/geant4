cmake_minimum_required(VERSION 3.16...3.21)
project(test27)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test27 test27.cc ${sources} ${headers})
target_include_directories(test27 PRIVATE include)
target_link_libraries(test27 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test27
    COMMAND test27 ${CMAKE_CURRENT_SOURCE_DIR}/test27.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
