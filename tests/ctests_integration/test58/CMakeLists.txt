cmake_minimum_required(VERSION 3.16...3.21)
project(test58)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test58 test58.cc ${sources} ${headers})
target_include_directories(test58 PRIVATE include)
target_link_libraries(test58 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test58
    COMMAND test58 ${CMAKE_CURRENT_SOURCE_DIR}/test58.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
