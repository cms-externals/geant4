# - Find needed modules/packages
find_package(Python3 3.9 QUIET COMPONENTS Interpreter)
find_package(Git QUIET)

#-------------------------------------------------------------------------------------------------
# - Missing/overspecified sources
#-------------------------------------------------------------------------------------------------
geant4_add_test(validate-sources
  COMMAND ${CMAKE_COMMAND} -P geant4_validate_sources.cmake
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  LABELS Continuous)

#-------------------------------------------------------------------------------------------------
# - Source code module checks
#-------------------------------------------------------------------------------------------------
if(Python3_FOUND)
  # Module checker command
  set(G4MODULE_CHECK_CMD ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/cmake/Modules/geant4_module_check.py -v -db ${PROJECT_BINARY_DIR}/G4ModuleInterfaceMap.csv)

  # Tests for the script itself
  # Note that output is only checked at a very basic level as it is dependent on the module composition
  geant4_add_test(validate-module-check-list
    COMMAND ${G4MODULE_CHECK_CMD} --list
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PASSREGEX "G4globman;G4intercoms"
    LABELS Continuous)
  geant4_add_test(validate-module-check-source
    COMMAND ${G4MODULE_CHECK_CMD} --source G4globman
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PASSREGEX "${PROJECT_SOURCE_DIR}/source/global/management"
    LABELS Continuous)
  geant4_add_test(validate-module-check-interface
    COMMAND ${G4MODULE_CHECK_CMD} --interface G4globman
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    PASSREGEX "G4String.hh;G4Version.hh"
    LABELS Continuous)

  # TODO: this fails on Windows:
  # G4Test stderr:
  # No module provides header 'G4String.hh'
  # System specific as passes on all other platforms
  if(NOT WIN32)
    geant4_add_test(validate-module-check-provides
      COMMAND ${G4MODULE_CHECK_CMD} --provides G4String.hh
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      PASSREGEX "G4globman"
      LABELS Continuous)
  endif()

  # Only run this test on standard builds where we know expected library name
  if(NOT __GEANT4_LIBRARY_DEFINITION_FILE)
    geant4_add_test(validate-module-check-library
      COMMAND ${G4MODULE_CHECK_CMD} --library G4globman
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      PASSREGEX "G4global"
      LABELS Continuous)
  endif()

  # Regression: check compiled and header-only modules can have usage requirements queried
  geant4_add_test(validate-module-check-usage-full
    COMMAND ${G4MODULE_CHECK_CMD} --usage-requirements G4globman
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LABELS Continuous)
  geant4_add_test(validate-module-check-usage-interface
    COMMAND ${G4MODULE_CHECK_CMD} --usage-requirements G4hepgeometry
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LABELS Continuous)

  # Actual module Validation tests
  geant4_add_test(validate-no-module-cycles
    COMMAND ${G4MODULE_CHECK_CMD} --find-cycles
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LABELS Continuous Nightly)

  # Windows checks have some issues with file encoding, e.g.
  # File "C:\Program Files\Python39\lib\encodings\cp1252.py", line 23, in decode
  #  return codecs.charmap_decode(input,self.errors,decoding_table)[0]
  #  UnicodeDecodeError: 'charmap' codec can't decode byte 0x9d in position 8133: character maps to <undefined>
  if(NOT WIN32)
    geant4_add_test(validate-module-consistency
      COMMAND ${G4MODULE_CHECK_CMD} --find-inconsistencies
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      LABELS Continuous Nightly)
  endif()
endif()


#-------------------------------------------------------------------------------------------------
# - Test against direct use of getenv for data directories
#-------------------------------------------------------------------------------------------------
if(UNIX AND Git_FOUND)
  geant4_add_test(validate-g4-find-datadir
    COMMAND /bin/sh -c "${GIT_EXECUTABLE} --git-dir=${PROJECT_SOURCE_DIR}/.git --work-tree=${PROJECT_SOURCE_DIR} grep 'getenv.*(.*G4.*DATA.*)' -- source/ >&2 || true"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    LABELS Continuous)
endif()

#-------------------------------------------------------------------------------------------------
# - Test for C++ Standard Features supported by current compiler
#-------------------------------------------------------------------------------------------------
# We *don't* use geant4_add_test because that can only build one target
# and not obvious that can be used on a whole project for VS (no "all" equivalent?)
add_test(NAME validate-cxx-features COMMAND ${CMAKE_CTEST_COMMAND}
  --build-and-test ${CMAKE_CURRENT_SOURCE_DIR}/checkCXXFeatures ${CMAKE_CURRENT_BINARY_DIR}/checkCXXFeatures
  --build-generator ${CMAKE_GENERATOR}
  --build-makeprogram ${CMAKE_MAKE_PROGRAM}
  --build-project checkCXXFeatures
  --build-config $<CONFIGURATION>
  --build-noclean
  --build-options --no-warn-unused-cli -DCCF_BUILD_FAILING_CHECKS=ON -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
  --test-command checkCXXFeatures )
# Runs on all
set_property(TEST validate-cxx-features PROPERTY LABELS Nightly Continuous)
