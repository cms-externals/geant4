#----------------------------------------------------------------------------
# Re-set target output directories to ensure targets are output to per-test
# subdirectories
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

foreach(_c ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER ${_c} _c_upper)
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_c_upper})
  unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_c_upper})
  unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_c_upper})
endforeach()

#----------------------------------------------------------------------------
# Add core test projects
# - Code/Library structure checks
add_subdirectory(code_checks)

# - Code/Library Unit Tests 
add_subdirectory(ctests_source)

# - Core system tests only if we're not testing a custom library modularization
#   scheme. Geant4Config.cmake does not yet auto-generate the list of targets.
if(__GEANT4_LIBRARY_DEFINITION_FILE)
  return()
endif()

# - Basic/Extended/Advanced Examples test builds and runs
add_subdirectory(ctests_examples)

#----------------------------------------------------------------------------
# Add in-detail test projects
# - Integration tests build and run
add_subdirectory(ctests_integration)

# - Verification tests build and run
add_subdirectory(ctests_verification)
