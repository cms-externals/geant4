cmake_minimum_required(VERSION 3.16...3.21)
project(test30)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test30 test30.cc ${sources} ${headers})
target_include_directories(test30 PRIVATE include)
target_link_libraries(test30 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test30
    COMMAND test30 ${CMAKE_CURRENT_SOURCE_DIR}/test30.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
