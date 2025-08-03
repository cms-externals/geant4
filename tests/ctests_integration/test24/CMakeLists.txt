cmake_minimum_required(VERSION 3.16...3.21)
project(test24)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test24 test24.cc ${sources} ${headers})
target_include_directories(test24 PRIVATE include)
target_link_libraries(test24 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test24
    COMMAND test24 ${CMAKE_CURRENT_SOURCE_DIR}/test24.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
