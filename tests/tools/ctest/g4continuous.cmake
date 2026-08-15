cmake_minimum_required(VERSION 3.16...3.27 FATAL_ERROR)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---mark continuous in build name-------------------------------------------
set(CTEST_BUILD_NAME c_$ENV{gitlabUserName}_$ENV{MergeRequestId}_${CTEST_BUILD_NAME})

#---Addional CTest settings-------------------------------------------------
#---Jenkins will merge the current master into a MR, so use the MR commit---
set(CTEST_UPDATE_VERSION_OVERRIDE "$ENV{MergeRequestLastCommit}")

#---Clean or not clean the Binary-------------------------------------------
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/is_docker)
  ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

#---CDash Retry Upload Settings---------------------------------------------
set(G4NIGHTLY_RETRY_ARGS
  RETRY_COUNT 5 
  RETRY_DELAY 120)

#---CTest commands----------------------------------------------------------
ctest_start(Continuous)

ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})
ctest_configure(BUILD ${CTEST_BINARY_DIRECTORY}
  SOURCE ${CTEST_SOURCE_DIRECTORY}
  OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_submit(PARTS Update Configure ${G4NIGHTLY_RETRY_ARGS})

ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
ctest_submit(PARTS Build ${G4NIGHTLY_RETRY_ARGS})

if("$ENV{BUILDOPTIONS}" MATCHES "GranularCheck")
  ctest_test(PARALLEL_LEVEL ${ncpu}
    INCLUDE_LABEL "Continuous"
    EXCLUDE "^test"
    INCLUDE "^validate")
else()
  ctest_test(PARALLEL_LEVEL ${ncpu}
    INCLUDE_LABEL "Nightly|Continuous"
    EXCLUDE "largeN$"
    INCLUDE "^test|^source|^validate$")
endif()

ctest_submit(PARTS Test Coverage MemCheck Notes ExtraFiles ${G4NIGHTLY_RETRY_ARGS})
