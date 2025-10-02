# We only need to run these tests on a subset of platforms
# UNIX only...
if(WIN32)
  return()
endif()
# ... except for macOS
if(APPLE)
  return()
endif()
# ... Shared libs only ...
if(NOT BUILD_SHARED_LIBS)
  return()
endif()
# ... without VecGeom solids ...
if(GEANT4_USE_USOLIDS)
  return()
endif()
# ... only for GCC > 7 ...
if((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8))
  return()
endif()
# ..without BoundsCheck...
if(BOUNDSCHECK)
  return()
endif()

# GNUMake system relies on its own environment setup script,
# so wrap any command in a script that sources this before
# exex-ing the command
configure_file(gnumake_wrapper.sh.in gnumake_wrapper.sh @ONLY)
set(GMB_WRAPPER "${CMAKE_CURRENT_BINARY_DIR}/gnumake_wrapper.sh")

# Basic test setup only, create parallel binary tree to avoid
# clashes with CMake tests
set(GMB_BASIC_SOURCE_DIR ${PROJECT_SOURCE_DIR}/examples/basic)
set(GMB_BASIC_BINARY_DIR ${PROJECT_BINARY_DIR}/examples/gnumake)
file(MAKE_DIRECTORY ${GMB_BASIC_BINARY_DIR})

set(GMB_COMMON_TESTS
  B1
  B2/B2a
  B2/B2b
  B3/B3a
  B3/B3b
  B4/B4a
  B4/B4b
  B4/B4c
  B4/B4d
  B5)

foreach(_gmb_test ${GMB_COMMON_TESTS})
  get_filename_component(_gmb_test_tag ${_gmb_test} NAME)

  set(_gmb_test_build "gnumake-bas-${_gmb_test_tag}-build")
  set(_gmb_test_run "gnumake-bas-${_gmb_test_tag}-run")

  set(_gmb_test_src_dir ${GMB_BASIC_SOURCE_DIR}/${_gmb_test})
  set(_gmb_test_bin_dir ${PROJECT_BINARY_DIR}/examples/gnumake/${_gmb_test})
  file(MAKE_DIRECTORY ${_gmb_test_bin_dir})

  set(_gmb_extra_arg)
  if(${_gmb_test} MATCHES "^B4")
    set(_gmb_extra_arg "-m")
  endif()
  string(SUBSTRING ${_gmb_test} 0 2 _gmb_test_macro)
  set(_gmb_test_macro ${_gmb_extra_arg} ${_gmb_test_src_dir}/example${_gmb_test_macro}.in)

  geant4_add_test(${_gmb_test_build} COMMAND ${GMB_WRAPPER} make -C ${_gmb_test_src_dir}
    WORKING_DIRECTORY ${_gmb_test_bin_dir})

  # GNUmake builds in parallel can potentially clash over creation of the G4TMP/G4SYSTEM directory
  # so we force the build part of the test to run in serial using a resource lock
  set_tests_properties(${_gmb_test_build} PROPERTIES RESOURCE_LOCK G4SystemMkDir)

  geant4_add_test(${_gmb_test_run} COMMAND ${GMB_WRAPPER} example${_gmb_test_tag} ${_gmb_test_macro}
    WORKING_DIRECTORY ${_gmb_test_bin_dir}
    DEPENDS ${_gmb_test_build})
endforeach()
