cmake_minimum_required(VERSION 3.16...3.21)
project(test52)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test52 test52.cc ${sources} ${headers})
target_include_directories(test52 PRIVATE include)
target_link_libraries(test52 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test52
    COMMAND test52 ${CMAKE_CURRENT_SOURCE_DIR}/test52.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

