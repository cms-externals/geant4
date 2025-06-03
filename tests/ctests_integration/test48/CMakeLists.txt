cmake_minimum_required(VERSION 3.16...3.21)
project(test48)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. --> test48 disabled")
  return()
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test48 test48.cc ${sources} ${headers})
target_include_directories(test48
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${INCLUDE_DIRECTORIES}
    ${EXTRA_INCLUDE_DIRECTORIES}
    ${ROOT_INCLUDE_DIRS})
target_link_libraries(test48
  PRIVATE
    ${EXTRA_LIBRARIES}
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES}
    ${XercesC_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts, input files to the build directory, i.e. the directory in which we
# build test48.
#
#  Note that it is not recommend to use GLOB to collect a list of
#  source files etc... from your source tree. If no CMakeLists.txt file
#  changes when a source is added or removed then the generated build
#  system cannot know when to ask CMake to regenerate.
file(GLOB TEST48_FILES RELATIVE ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/test48.*minus*
  ${PROJECT_SOURCE_DIR}/ctest48-picapture-Pb.config
  ${PROJECT_SOURCE_DIR}/*.C
  ${PROJECT_SOURCE_DIR}/antiproton/*.dat
  ${PROJECT_SOURCE_DIR}/kminus/*.dat
  ${PROJECT_SOURCE_DIR}/muminus/*.dat
  ${PROJECT_SOURCE_DIR}/piminus/*.dat
  test48.antiproton
  ${PROJECT_SOURCE_DIR}/*.sh
  ${PROJECT_SOURCE_DIR}/batch/*.sh
  ${PROJECT_SOURCE_DIR}/g4val-upload-scripts/gen_*.sh
  README)

foreach(_file ${TEST48_FILES})
  configure_file(${PROJECT_SOURCE_DIR}/${_file} ${PROJECT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test48-piminus
    COMMAND test48 ${CMAKE_CURRENT_SOURCE_DIR}/ctest48-picapture-Pb.config
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test48-muminus
    COMMAND test48 ${CMAKE_CURRENT_SOURCE_DIR}/ctest48-mucapture-Si.config
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
