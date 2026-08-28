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
# - Core
geant4_add_unit_test(G4toolsCompile.test.cc)
geant4_add_unit_test(lina/vec.test.cc)
geant4_add_unit_test(lina/box3.test.cc)
geant4_add_unit_test(lina/line.test.cc)
geant4_add_unit_test(lina/plane.test.cc)
geant4_add_unit_test(lina/clip.test.cc)
geant4_add_unit_test(lina/mat4.test.cc)
geant4_add_unit_test(lina/qrot.test.cc)

geant4_add_unit_test(buf2lines.test.cc)
geant4_add_unit_test(charp_out.test.cc)
geant4_add_unit_test(get_lines.test.cc)
geant4_add_unit_test(num2s.test.cc)
geant4_add_unit_test(rcmp.test.cc)
geant4_add_unit_test(S_STRING.test.cc)
geant4_add_unit_test(schar.test.cc)
geant4_add_unit_test(sout.test.cc)
geant4_add_unit_test(sprintf.test.cc)
geant4_add_unit_test(srep.test.cc)
geant4_add_unit_test(sto.test.cc)
geant4_add_unit_test(strip.test.cc)
geant4_add_unit_test(tokenize.test.cc)
geant4_add_unit_test(tos.test.cc)
geant4_add_unit_test(touplow.test.cc)

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

#-----------------------------------------------------------------------------#
# - Micro Benchmarks
#-----------------------------------------------------------------------------#
if(GEANT4_USE_benchmark)
  # We don't add these as a test just yet
  add_executable(sto_bench sto.bench.cc)
  target_link_libraries(sto_bench benchmark::benchmark Geant4::G4tools)
endif()
