cmake_minimum_required(VERSION 3.16...3.21)
project(test74)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test74 test74.cc ${sources} ${headers})
target_include_directories(test74 PRIVATE include)
target_link_libraries(test74 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test74
    COMMAND test74 ${CMAKE_CURRENT_SOURCE_DIR}/test74.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test74-largeN
    COMMAND test74 ${CMAKE_CURRENT_SOURCE_DIR}/test74.large_N.in
    DEPENDS test74
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
