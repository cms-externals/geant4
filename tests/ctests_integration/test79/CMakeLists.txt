#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(test79)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build test79. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TESTtest79_SCRIPTS
  init_vis.mac
  vis.mac
  run.mac
  Si_planes_111.dat
  )

foreach(_script ${TESTtest79_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(test79 test79.cc ${sources} ${headers})
target_include_directories(test79 PRIVATE include)
target_link_libraries(test79 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test79
    COMMAND test79 run.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
