# Locate Pythia6 library
# in a directory defined via PYTHIA6 environment variable
#
# Defines:
#  PYTHIA6_FOUND
#  PYTHIA6_LIBRARIES
set(TEST_PYTHIA6_ROOT_DIR  "" ${PYTHIA6_ROOT_DIR})
if(TEST_PYTHIA6_ROOT_DIR STREQUAL "")
  if(DEFINED ENV{PYTHIA6_ROOT_DIR})
    set(PYTHIA6_ROOT_DIR  $ENV{PYTHIA6_ROOT_DIR})
  else()
    set(PYTHIA6_ROOT_DIR  "/usr")
  endif()
endif()

find_library(PYTHIA6_LIBRARY
  NAMES Pythia6 pythia6 pythia6-$ENV{PYTHIA6_VERSION}
  HINTS $ENV{PYTHIA6} $ENV{PYTHIA6}/lib $ENV{PYTHIA6}/lib64 ${PYTHIA6_ROOT_DIR}/lib ${PYTHIA6_ROOT_DIR}/lib64)

set(PYTHIA6_LIBRARIES ${PYTHIA6_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set PYTHIA6_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Pythia6 DEFAULT_MSG PYTHIA6_LIBRARIES)

mark_as_advanced(PYTHIA6_FOUND PYTHIA6_LIBRARIES)
