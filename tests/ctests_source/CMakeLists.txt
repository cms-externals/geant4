#-----------------------------------------------------------------------
# Per-category tests
# Major reason for having these separate from the toolkit is that it's
# easier to implement the needed module->library transform linking here.
# tests/ is processed after source/ so all module grouping is defined
# by that point.
# Otherwise, would have to do a little chicken/egg processing
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# - Setup Dependencies
#-----------------------------------------------------------------------
# Gtest chosen as very common/familiar library
find_package(GTest QUIET)

# If we don't find GTest, add it as an external (so, this would eventually go into externals or
# as part of the testing setup)
if(GTest_FOUND)
  message(STATUS "Found GTest at ${GTest_DIR}")
else()
  message(STATUS "Using FetchContent to import GoogleTest locally")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    # v1.14.0
    URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip)
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  # Don't install gtest itself
  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  FetchContent_MakeAvailable(googletest)
endif()

#-----------------------------------------------------------------------
# G4tools tests can require MPI
# - At present, this is only for examples, but we set things up here
#   for convenience and future proofing it G4mpi becomes a toolkit component.
# Relies on MPI being found _before_ being called
function(_geant4addunittest_mpi_cmd outvar np test_exe)
  if(MPI_FOUND)
    if(np GREATER 1)
      set(_cmd "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} "${np}"
        ${MPIEXEC_PREFLAGS} ${test_exe} ${MPIEXEC_POSTFLAGS} ${ARGN})
    else()
      set(_cmd ${test_exe} ${ARGN})
    endif()
  else()
    set(_cmd "${test_exe}" ${ARGN})
  endif()
  set(${outvar} "${_cmd}" PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# - API for setting up/declaring GTest-based unit tests
#-----------------------------------------------------------------------
# - Temporary location, eventually to be refactored into G4TestAPI module
#[=======================================================================[.rst:
.. command:: geant4_setup_unit_tests

  Set dependencies for the unit tests in the current CMakeLists file. It
  should be called before any tests are defined

    geant4_setup_unit_tests(
      [LINK_LIBRARIES <target> [<target> ...]]
      )

    ``LINK_LIBRARIES``
      Link test executables to these libraries by default

#]=======================================================================]
function(geant4_setup_unit_tests)
  cmake_parse_arguments(PARSE
    ""
    ""
    "LINK_LIBRARIES"
    ${ARGN}
  )

  # Set the directory relative to the "test" directory
  file(RELATIVE_PATH PARSE_DIR "${PROJECT_SOURCE_DIR}/tests/ctests_source" "${CMAKE_CURRENT_SOURCE_DIR}")

  # Set special variables
  foreach(_var LINK_LIBRARIES DIR)
    set(GEANT4_UNIT_TEST_${_var} "${PARSE_${_var}}" PARENT_SCOPE)
  endforeach()
endfunction()

#[=======================================================================[.rst:
.. command:: geant4_add_unit_test

  Add a C++ GoogleTest test::

    geant4_add_unit_test(
      <filename>
      [LINK_LIBRARIES lib1 [...]]
      [NP n1]
      [SOURCES file1 [...]]
      )

    ``<filename>``
      Test source file name relative to current directory

    ``LINK_LIBRARIES``
      Extra libraries to link to. By default, unit tests will link against the
      package's current library.

    ``NP``
      The number of processes to launch via MPI for this unit test. Only relevant
      if the test needs MPI.

    ``SOURCES``
      Additional sources to compile into the test executable

  The target and test name are derived from the filename (without extension)
  and the directory path of that file relative to `ctests_source`.

#]=======================================================================]
function(geant4_add_unit_test SOURCE_FILE)
  cmake_parse_arguments(PARSE
    ""
    ""
    "LINK_LIBRARIES;NP;SOURCES"
    ${ARGN}
  )
  if(PARSE_UNPARSED_ARGUMENTS)
    message(SEND_ERROR "Unknown keywords given to geant4_add_unit_test(): "
            "\"${PARSE_UNPARSED_ARGUMENTS}\"")
  endif()

  if(NOT PARSE_NP)
    set(PARSE_NP 1)
  endif()

  get_filename_component(_BASENAME "${SOURCE_FILE}" NAME_WE)
  get_filename_component(_BASEDIR "${SOURCE_FILE}" DIRECTORY)
  if(_BASEDIR)
    set(_BASENAME "${_BASEDIR}/${_BASENAME}")
  endif()

  set(_TARGET "${_BASENAME}")
  set(_OUTPUT_NAME "${_TARGET}")
  
  if(GEANT4_UNIT_TEST_DIR)
    set(_TARGET "${GEANT4_UNIT_TEST_DIR}/${_TARGET}")
    set(_PREFIX "${GEANT4_UNIT_TEST_DIR}/")
  endif()
  string(REGEX REPLACE "[^a-zA-Z0-9_]+" "_" _TARGET "${_TARGET}")
  string(REGEX REPLACE "[^a-zA-Z0-9_]+" "_" _OUTPUT_NAME "${_OUTPUT_NAME}")

  add_executable(${_TARGET} "${SOURCE_FILE}" ${PARSE_SOURCES})
  set_property(TARGET ${_TARGET} PROPERTY OUTPUT_NAME "${_OUTPUT_NAME}")
  geant4_test_link_libraries(${_TARGET} 
    PRIVATE 
      ${GEANT4_UNIT_TEST_LINK_LIBRARIES}
      ${PARSE_LINK_LIBRARIES}
      GTest::gtest_main
  )

  set(_EXE_NAME "$<TARGET_FILE:${_TARGET}>")
  set(_TEST_NAME "source/${_PREFIX}${_BASENAME}")
  _geant4addunittest_mpi_cmd(_test_cmd "${PARSE_NP}" "${_EXE_NAME}")

  add_test(NAME ${_TEST_NAME} COMMAND ${_test_cmd})
  set_tests_properties(${_TEST_NAME} PROPERTIES LABELS "Nightly;Continuous")
  
  if(WIN32)
    # Requires CMake 3.22, which is fine in CI, but needs improvement
    set_tests_properties(${_TEST_NAME}
      PROPERTIES 
        ENVIRONMENT_MODIFICATION "PATH=path_list_prepend:$<SHELL_PATH:$<TARGET_FILE_DIR:GTest::gtest>>"
    )
  endif()  
endfunction()

#-----------------------------------------------------------------------
# - Recurse into per-category tests
#-----------------------------------------------------------------------
add_subdirectory(externals)
add_subdirectory(global)
add_subdirectory(intercoms)
add_subdirectory(materials)
