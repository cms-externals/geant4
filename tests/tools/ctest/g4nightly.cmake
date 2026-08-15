cmake_minimum_required(VERSION 3.16...3.27 FATAL_ERROR)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---Addional CTest settings-------------------------------------------------
#---Enable Memcheck builds if requested
if(DEFINED ENV{WITH_MEMCHECK})
  set(WITH_MEMCHECK TRUE)
  # Always valgrind for now (NB: should also check platform, but that's nominally
  # done through the WITH_MEMCHECK env var
  find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)
  message(STATUS "running with memcheck command ${CTEST_MEMORYCHECK_COMMAND}" )
  set(CTEST_MEMORYCHECK_TYPE "Valgrind")

  # NB: Because tests are lauched through a CMake driver script, must use
  # --trace-children to check the actual test executable
  set(CTEST_MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --track-origins=yes")

  if(DEFINED ENV{VALGRIND_ROOTSUPP})
    message(STATUS "using root supression file $ENV{VALGRIND_ROOTSUPP}")
    set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE "$ENV{VALGRIND_ROOTSUPP}")
  endif()
endif()

# - Filter tests that timeout when run in BoundsCheck
if("$ENV{BUILDOPTIONS}" MATCHES "BoundsCheck")
  set(BOUNDSCHECK_EXCLUDE_ARG EXCLUDE_LABEL "ExcludeFromBoundsCheck")
endif()

#---Clean or not clean the Binary-------------------------------------------
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/is_docker )
   ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

#---CDash Retry Upload Settings---------------------------------------------
set(G4NIGHTLY_RETRY_ARGS
  RETRY_COUNT 5 
  RETRY_DELAY 120)

#---CTest commands----------------------------------------------------------
ctest_start("Nightly" TRACK $ENV{CdashTrack})

ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})
ctest_configure(BUILD ${CTEST_BINARY_DIRECTORY}
  SOURCE  ${CTEST_SOURCE_DIRECTORY}
  OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_submit(PARTS Update Configure ${G4NIGHTLY_RETRY_ARGS})

ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY} NUMBER_ERRORS G4BUILD_NUMBER_ERRORS)
ctest_submit(PARTS Build ${G4NIGHTLY_RETRY_ARGS})

# Temporarily here to test behaviour in Continuous builds.
# Eventually only in Nightlies
# On Windows and with no build errors, deploy any downloaded datasets
if(WIN32 AND (G4BUILD_NUMBER_ERRORS EQUAL 0))
  execute_process(COMMAND ${CMAKE_COMMAND} --install ${CTEST_BINARY_DIRECTORY} --component Data
    COMMAND_ECHO STDOUT
    RESULT_VARIABLE G4DEPLOYDATA_RES
    ERROR_VARIABLE G4DEPLOYDATA_ERR
  )
  if(G4DEPLOYDATA_ERR)
    message(STATUS "Deployment of datasets on Windows failed: ${G4DEPLOYDATA_ERR}")
    # TODO: Likely need a cleanup mechanism in case of errors, but pure "uninstall" won't
    # work since we don't install everything.
  endif()
endif()  

if(WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  # Only build (not run) the example tests that memcheck will exercise
  # This avoids double-running tests in this specific build which are exercised in other configurations
  ctest_test(PARALLEL_LEVEL ${ncpu}
    INCLUDE "example-(bas|ext).*-build$"
    EXCLUDE_LABEL "ExcludeFromMemCheck")
else()
  ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE_LABEL "Nightly" ${BOUNDSCHECK_EXCLUDE_ARG})
endif()
ctest_submit(PARTS Test ${G4NIGHTLY_RETRY_ARGS})

if(WITH_MEMCHECK AND CTEST_MEMORYCHECK_COMMAND)
  message(STATUS "memory checking using ${CTEST_MEMORYCHECK_COMMAND}" )
  # Increase default timeout, and set env var G4TestDriver.cmake uses for launching tests
  set(CTEST_TEST_TIMEOUT 7000)
  set(ENV{CTEST_TIMEOUT} ${CTEST_TEST_TIMEOUT})
  # Valgrind may need additional flags, e.g. --fair-sched=yes, or may need to limit
  # number of threads to get good performance. For now, follow old GNUmake practice
  # and only run with a Serial RunManager.
  # Take simplest option as used in old GNUmake builds
  set(ENV{G4FORCE_RUN_MANAGER_TYPE} Serial)
  ctest_memcheck(PARALLEL_LEVEL ${ncpu}
    INCLUDE "example-(bas|ext)"
    EXCLUDE "-build$"
    EXCLUDE_LABEL "ExcludeFromMemCheck")
endif()

ctest_submit(PARTS Coverage MemCheck Notes ExtraFiles ${G4NIGHTLY_RETRY_ARGS})
