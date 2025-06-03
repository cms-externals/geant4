cmake_minimum_required(VERSION 3.16...3.21)
project(test70)

# Vis only if testing is disabled
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" OFF)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test70 test70.cc ${sources} ${headers})
target_include_directories(test70 PRIVATE include)
target_link_libraries(test70 PRIVATE ${Geant4_LIBRARIES})

#-------------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
foreach(_script long.mac short-mt.mac short.mac)
  configure_file(${_script} ${PROJECT_BINARY_DIR}/ COPYONLY)
endforeach()

# ----------------------------------------------------------------------------------------------------
# --- Test definitions
# ----------------------------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test70
    COMMAND test70 -mac ${CMAKE_CURRENT_SOURCE_DIR}/short.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    LABELS Nightly Continuous)

  geant4_add_test(test70-mt
    COMMAND test70 -mt 2 -mac ${CMAKE_CURRENT_SOURCE_DIR}/short-mt.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    LABELS Nightly Continuous)

  # ----------------------------------------------------------------------------------------------------
  # --- Nightly ONLY
  # ----------------------------------------------------------------------------------------------------

  #geant4_add_test(test70-long
  #  COMMAND test70 -mac ${CMAKE_CURRENT_SOURCE_DIR}/long.mac
  #  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  #  LABELS Nightly)

  #geant4_add_test(test70-long-mt
  #  COMMAND test70 -mt 10 -mac ${CMAKE_CURRENT_SOURCE_DIR}/long.mac
  #  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  #  LABELS Nightly)
endif()
