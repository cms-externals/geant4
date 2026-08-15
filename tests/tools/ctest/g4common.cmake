#---Utility Macros----------------------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4macros.cmake)

#---Make sure that VERBOSE is OFF to avoid screwing up the build performance
unset(ENV{VERBOSE})

#---Primary Configuration---------------------------------------------------
#---Dashboard
set(CTEST_PROJECT_NAME "Geant4")
set(CTEST_NIGHTLY_START_TIME "00:00:00 CET")
set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "cdash.cern.ch")
set(CTEST_DROP_LOCATION "/submit.php?project=Geant4")
set(CTEST_DROP_SITE_CDASH TRUE)

#---Run/Source/Binary Directories
get_pwd(pwd)
# Can probably use one of:
message(STATUS "cpwd = ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "cpwd = ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "Running CTest from '${pwd}'")

set(CTEST_SOURCE_DIRECTORY "$ENV{SOURCE}")
if(NOT CTEST_SOURCE_DIRECTORY)
  set(CTEST_SOURCE_DIRECTORY "${pwd}/geant4")
endif()

set(CTEST_BINARY_DIRECTORY "$ENV{BINARY}")
if(NOT CTEST_BINARY_DIRECTORY)
  set(CTEST_BINARY_DIRECTORY "${pwd}/build")
endif()

#--- Site
set(CTEST_SITE "$ENV{CTEST_SITE}")
if(NOT CTEST_SITE)
  cmake_host_system_information(RESULT host QUERY HOSTNAME)
  set(CTEST_SITE "${host}")
endif()
message(STATUS "Running build and test on ${host}")

#--- Notes
set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME})

# Jenkins.log and jenkins.err keep Jenkins log and error, see jk-all; do not exist for Windows. 
if(EXISTS ${pwd}/jenkins.err)
  message(STATUS "jenkins.err file exists")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${pwd}/jenkins.err)
endif()
if(EXISTS ${pwd}/jenkins.log)
  message(STATUS "jenkins.log file exists")
  set(CTEST_NOTES_FILES ${CTEST_NOTES_FILES} ${pwd}/jenkins.log)
endif()
 
#---General Step Configuration---------------------------------------------------
#---Update: Discover commit we're building only in "update" step
find_package(Git)
set(CTEST_UPDATE_VERSION_ONLY ON)
set(CTEST_UPDATE_COMMAND "${GIT_EXECUTABLE}")

#--- Configure
set(CTEST_START_WITH_EMPTY_BINARY_DIRECTORY_ONCE 1)
set(CTEST_CONFIG_OPTIONS
  -DGEANT4_ENABLE_TESTING=ON
  -DGEANT4_INSTALL_DATA=ON
  -DGEANT4_USE_GDML=ON
  -DGEANT4_BUILD_BUILTIN_BACKTRACE=ON
  $ENV{G4_XOPTS})


#--- Build
cmake_host_system_information(RESULT ncpu QUERY NUMBER_OF_PHYSICAL_CORES)

if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64" AND CMAKE_SYSTEM_NAME MATCHES "Linux")
   cmake_host_system_information(RESULT ncpu QUERY NUMBER_OF_LOGICAL_CORES)
   message(STATUS "Using number of LOGICAL cores on aarch64/Linux")
endif()

#---Allow to reduce parallelism, useful for MTmax builds
# limit for build using make -j ...
set(opt_make_j "$ENV{CMAKE_BUILD_PARALLEL_LEVEL}")
if(NOT opt_make_j)
  set(opt_make_j ${ncpu})
endif()

# and for number tests running in parallel
if(NOT "$ENV{MAX_CPUS_USE}" STREQUAL "")
  message(STATUS "limiting parallel running of tests to $ENV{MAX_CPUS_USE} tests.")
  set(ncpu "$ENV{MAX_CPUS_USE}")
endif()

# Windows static MT build cannot run tests fully parallel.....
if(WIN32 AND "$ENV{BUILDOPTIONS}" MATCHES "static" AND "$ENV{THREAD}" MATCHES "MT")
  math(EXPR ncpu "${ncpu}/2")
  message(STATUS "limiting parallel running of tests to ${ncpu}.")
endif()

 
#---Build Configuration--------------------------------------------------------
#---Generator
get_configuration_tag(tag)

if(NOT WIN32)
  set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
  set(CTEST_BUILD_COMMAND "make -s -i -j${opt_make_j}")
else(WIN32)
  if(tag MATCHES x64)
    set(win64 " Win64")
  endif()
  # be4 adding a new generator, make sure that cmake knows about this....
  if(tag MATCHES vc145)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 18 2026")
  elseif(tag MATCHES vc144)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 17 2022")
  elseif(tag MATCHES vc143)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 17 2022")
  elseif(tag MATCHES vc142)
    set(CTEST_CMAKE_GENERATOR "Visual Studio 16 2019")
  else()
    set(CTEST_CMAKE_GENERATOR "NMake Makefiles")
  endif()
endif()

#---Build Type and Name
list(APPEND CMAKE_CONFIGURATION_TYPES Release Debug RelWithDebInfo MinSizeRel TestRelease FullRelWithDebInfo)

if("$ENV{BUILDTYPE}" STREQUAL "" OR "$ENV{BUILDTYPE}" STREQUAL "RelWithDebInfo")
  set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
  set(CTEST_BUILD_NAME ${tag})
elseif("${CMAKE_CONFIGURATION_TYPES}" MATCHES "$ENV{BUILDTYPE}")
  set(CTEST_BUILD_CONFIGURATION "$ENV{BUILDTYPE}")
  set(CTEST_BUILD_NAME ${tag}-$ENV{BUILDTYPE})
else()
  set(CTEST_BUILD_CONFIGURATION "RelWithDebInfo")
  set(CTEST_BUILD_NAME ${tag}-$ENV{BUILDTYPE})
endif()

if(NOT "$ENV{THREAD}" STREQUAL "")
  set(CTEST_BUILD_NAME ${CTEST_BUILD_NAME}-$ENV{THREAD})
endif()

if(NOT "$ENV{BUILDOPTIONS}" STREQUAL "")
  set(CTEST_BUILD_NAME ${CTEST_BUILD_NAME}-$ENV{BUILDOPTIONS})
endif()

set(CTEST_CONFIGURATION_TYPE "${CTEST_BUILD_CONFIGURATION}")

#---Custom CTest settings---------------------------------------------------
set(CTEST_CUSTOM_TESTS_IGNORE test19 test29 test39 test49 test47)
set(CTEST_CUSTOM_MAXIMUM_FAILED_TEST_OUTPUT_SIZE "100000")
set(CTEST_CUSTOM_MAXIMUM_PASSED_TEST_OUTPUT_SIZE "10000")
if(WIN32)
  set(CTEST_CUSTOM_TESTS_IGNORE ${CTEST_CUSTOM_TESTS_IGNORE} example-ext-geometry-olap)
  set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
    "Ranlux64Engine.+: warning C4293:"
	  "SystemOfUnits.+: warning C4005: 'pascal'"
    ": warning LNK4221:")
elseif(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "aarch64" AND CMAKE_SYSTEM_NAME MATCHES "Linux")
   message(STATUS "g4macros.cmake: disable warning for std::pair(double,double)...")
   set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
      "when C[+]+17 is enabled changed to match C[+]+14 in GCC 10.1")
else()
  set(CTEST_CUSTOM_WARNING_EXCEPTION ${CTEST_CUSTOM_WARNING_EXCEPTION}
    "warning: ignoring return value of"
    "clang: warning: argument unused" 
    "include/xercesc/util/regx/Token.hpp"
    "note: variable tracking size limit exceeded with -fvar-tracking-assignments"
	  "warning: assuming signed overflow does not occur when assuming that ")
endif()

set(CTEST_TEST_TIMEOUT "$ENV{CTEST_TIMEOUT}")
if(NOT CTEST_TEST_TIMEOUT)
  set(CTEST_TEST_TIMEOUT 1500)
endif()

 
#---Set Runtime environment-------------------------------------------------
if(WIN32) 
  if(NOT CTEST_CMAKE_GENERATOR MATCHES Makefiles)
    set(_cfg /${CTEST_BUILD_CONFIGURATION})
  endif()
  set(ENV{PATH} "${CTEST_BINARY_DIRECTORY}/BuildProducts/${_cfg}/bin;$ENV{PATH}")
endif()

#---Workaround network issues on macOS machines-----------------------------
if(APPLE)
  set(ENV{G4GDML_SCHEMA_FILE} "${CTEST_SOURCE_DIRECTORY}/source/persistency/gdml/schema/gdml.xsd")
endif()
