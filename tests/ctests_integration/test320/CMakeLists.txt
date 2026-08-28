cmake_minimum_required(VERSION 3.16...3.21)
project(test320)

# Find Geant4 package, activating all available UI and Vis drivers by default
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test320 test320.cc ${sources} ${headers})
target_include_directories(test320 PRIVATE include)
target_link_libraries(test320 PRIVATE ${Geant4_LIBRARIES})

# Copy all macros to the build directory
set(EXAMPLEB1_SCRIPTS
  clean-files.sh
  init_vis.mac
  run-csv.mac
  run-hdf5.mac
  run-root.mac
  run-root-merge.mac
  run-xml.mac
  run.mac
  vis.mac
)

foreach(_macro ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_macro}
    ${PROJECT_BINARY_DIR}/${_macro}
    COPYONLY
    )
endforeach()

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test320-csv
    COMMAND test320 -m run-csv.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  if(Geant4_hdf5_FOUND AND
    (NOT GEANT4_TEST_ENVIRONMENT MATCHES "G4RUN_MANAGER_TYPE=TBB"))
    geant4_add_test(test320-hdf5
      COMMAND test320 -m run-hdf5.mac
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()

  geant4_add_test(test320-root
    COMMAND test320 -m run-root.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  if(GEANT4_BUILD_MULTITHREADED)
    geant4_add_test(test320-root-merge
      COMMAND test320 -m run-root-merge.mac -merge on
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()

  geant4_add_test(test320-xml
    COMMAND test320 -m run-xml.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
