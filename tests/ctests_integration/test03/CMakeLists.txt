cmake_minimum_required(VERSION 3.16...3.21)
project(test03)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executables, one per output type (CSV, HDF5, ROOT, XML, GENERIC)
add_executable(test03-Csv test03.cc ${sources} ${headers})
target_compile_definitions(test03-Csv PRIVATE TEST_ANALYSIS_CSV)
target_include_directories(test03-Csv PRIVATE include)
target_link_libraries(test03-Csv PRIVATE ${Geant4_LIBRARIES})

if(Geant4_hdf5_FOUND)
  add_executable(test03-Hdf5 test03.cc ${sources} ${headers})
  target_compile_definitions(test03-Hdf5 PRIVATE TEST_ANALYSIS_HDF5)
  target_include_directories(test03-Hdf5 PRIVATE include)
  target_link_libraries(test03-Hdf5 PRIVATE ${Geant4_LIBRARIES})
endif()

add_executable(test03-Root test03.cc ${sources} ${headers})
target_compile_definitions(test03-Root PRIVATE TEST_ANALYSIS_ROOT)
target_include_directories(test03-Root PRIVATE include)
target_link_libraries(test03-Root PRIVATE ${Geant4_LIBRARIES})

add_executable(test03-Xml test03.cc ${sources} ${headers})
target_compile_definitions(test03-Xml PRIVATE TEST_ANALYSIS_XML)
target_include_directories(test03-Xml PRIVATE include)
target_link_libraries(test03-Xml PRIVATE ${Geant4_LIBRARIES})

add_executable(test03-Generic test03.cc ${sources} ${headers})
target_compile_definitions(test03-Generic PRIVATE TEST_ANALYSIS_GENERIC)
target_include_directories(test03-Generic PRIVATE include)
target_link_libraries(test03-Generic PRIVATE ${Geant4_LIBRARIES})

# Copy scripts to build area
set(TEST_ANALYSIS_SCRIPTS
  testAnalysis.in
  gui.mac
  init.mac
  init_vis.mac
  h1.mac
  h2.mac
  h3.mac
  p1.mac
  p2.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${TEST_ANALYSIS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/macros/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# Define the tests, if available
if(GEANT4_ENABLE_TESTING)
  # CSV
  geant4_add_test(test03-Csv-write
    COMMAND test03-Csv -m run2.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test03-Csv-read
    COMMAND test03-Csv -m run2.mac -read on
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test03-Csv-write)

  # HDF5
  if(TARGET test03-Hdf5)
    geant4_add_test(test03-Hdf5-write
      COMMAND test03-Hdf5  -m run2.mac
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Hdf5-read
      COMMAND test03-Hdf5  -m run2.mac -read on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS test03-Hdf5-write)
  endif()

  # ROOT
  if(GEANT4_BUILD_MULTITHREADED)
    geant4_add_test(test03-Root-write
      COMMAND test03-Root -m run2.mac -merge on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-read
      COMMAND test03-Root -m run2.mac -merge on -read on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS test03-Root-write)
  else()
    geant4_add_test(test03-Root-write
      COMMAND test03-Root -m run2.mac
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-read
      COMMAND test03-Root -m run2.mac -read on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS test03-Root-write)
  endif()

  if(GEANT4_BUILD_MULTITHREADED)
    geant4_add_test(test03-Root-writeCw
      COMMAND test03-Root -m run2.mac -o testAnalysisCw -merge on -rowWise off -rowMode off
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-writeRwextra
      COMMAND test03-Root -m run2.mac -o testAnalysisRw -merge on -rowWise on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-write-mt-merge1
      COMMAND test03-Root -m run2.mac -o testAnalysisNf1 -merge on -nfiles 1
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-write-mt-merge2
      COMMAND test03-Root -m run2.mac -o testAnalysisNf2 -merge on -nfiles 2
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-write-mt-merge-max
      COMMAND test03-Root -m run2.mac -o testAnalysisNfmax -merge on -nfiles max
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Root-write-mt-no-merge
      COMMAND test03-Root -m run2.mac -o testAnalysisNoMerge -merge off
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()

  # XML
  geant4_add_test(test03-Xml-write
    COMMAND test03-Xml -m run2.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test03-Xml-read
    COMMAND test03-Xml -m run2.mac -read on
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test03-Xml-write)

  # Generic - ROOT output
  if(GEANT4_BUILD_MULTITHREADED)
    geant4_add_test(test03-Generic-write
      COMMAND test03-Generic -m run2.mac -merge on -o testAnalysisG -plot on -ascii on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
    geant4_add_test(test03-Generic-read
      COMMAND test03-Generic -m run2.mac -merge on -read on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS test03-Generic-write)
  else()
    geant4_add_test(test03-Generic-write
      COMMAND test03-Generic -m run2.mac -o testAnalysisG -plot on -ascii on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
    geant4_add_test(test03-Generic-read
      COMMAND test03-Generic -m run2.mac -read on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS test03-Generic-write)
  endif()

  if(GEANT4_BUILD_MULTITHREADED)
    geant4_add_test(test03-Generic-writeCw
      COMMAND test03-Generic -m run2.mac -o testAnalysisGCw -merge on -rowWise off -rowMode off
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Generic-writeRwextra
      COMMAND test03-Generic -m run2.mac -o testAnalysisGRw -merge on -rowWise on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Generic-write-mt-merge1
      COMMAND test03-Generic -m run2.mac -o testAnalysisGNf1 -merge on -nfiles 1
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Generic-write-mt-merge2
      COMMAND test03-Generic -m run2.mac -o testAnalysisGNf2 -merge on -nfiles 2
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Generic-write-mt-merge-max
      COMMAND test03-Generic -m run2.mac -o testAnalysisGNfmax -merge on -nfiles max
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(test03-Generic-write-mt-no-merge
      COMMAND test03-Generic -m run2.mac -o testAnalysisGNoMerge -merge off
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()
endif()
