#-----------------------------------------------------------------------------#
# - Setup
#-----------------------------------------------------------------------------#
# Always link to ZLIB and EXPAT as these are always present
geant4_setup_unit_tests(LINK_LIBRARIES Geant4::G4tools ${G4ZLIB_LIBRARIES} ${G4EXPAT_LIBRARIES})

# Quick and dirty way for all tests to use "detail" headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Optional target for tests needing HDF5
add_library(G4ToolsHDF5 INTERFACE)
target_compile_definitions(G4ToolsHDF5
  INTERFACE
    $<$<BOOL:${HDF5_FOUND}>:TOOLS_USE_HDF5>
)
target_link_libraries(G4ToolsHDF5
  INTERFACE
    $<$<BOOL:${HDF5_FOUND}>:${HDF5_LIBRARIES}>
)

#-----------------------------------------------------------------------------#
# - Tests
#-----------------------------------------------------------------------------#
geant4_add_unit_test(G4toolsCompile.test.cc)
geant4_add_unit_test(G4tools_get_lines.test.cc)

# Needs MPI which is only in an example at present, so we find but don't
# require it. These tests currently have a runtime issue when built with
# the Intel compiler, so are also disabled there temporarily.
find_package(MPI QUIET COMPONENTS CXX)
if(MPI_FOUND AND NOT (CMAKE_CXX_COMPILER_ID MATCHES "^Intel.*"))
  geant4_add_unit_test(mpi_histos.cc
    LINK_LIBRARIES MPI::MPI_CXX
    NP 3
  )
  geant4_add_unit_test(mpi_ntuple.cc
    LINK_LIBRARIES MPI::MPI_CXX
    NP 3
  )
endif()
geant4_add_unit_test(mt_ntuple_ref.cc)
geant4_add_unit_test(mt_ntuple_seq.cc)
geant4_add_unit_test(mt_ntuple.cc)

# - The combined "utest" tests
set(_g4tools_utest_sources
  utest/axml.cc
  utest/columns.cc
  utest/csv_histo.cc
  utest/csv_ntuple.cc
  utest/histo.cc
  utest/root_histo.cc
  utest/root_ntuple.cc
  $<$<BOOL:${HDF5_FOUND}>:
    utest/hdf5_histo.cc
    utest/hdf5_ntuple.cc
    utest/hdf5.cc
  >
    # These require external data files, to be looked at later as the
    # files are quite large.
    #utest/lego_offscreen_zb.cc
    #utest/viewplot_ttf.cc
    #utest/viewplot.cc 
)

geant4_add_unit_test(utest.cc 
  SOURCES ${_g4tools_utest_sources}
  LINK_LIBRARIES G4ToolsHDF5
)
