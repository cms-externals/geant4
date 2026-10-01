#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(radiobiology)

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
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(radiobiology radiobiology.cc ${sources} ${headers})
target_link_libraries(radiobiology ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy macros, input data, and analysis scripts to the build directory. This
# allows the executable and post-processing tools to run from that directory.
#
file(GLOB macros RELATIVE ${PROJECT_SOURCE_DIR} macro/*.mac)
file(GLOB_RECURSE data RELATIVE ${PROJECT_SOURCE_DIR} data/*)
file(GLOB analysis RELATIVE ${PROJECT_SOURCE_DIR} analysis/*.py)

foreach(_file ${macros} ${data} ${analysis})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_file}
    ${PROJECT_BINARY_DIR}/${_file}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS radiobiology DESTINATION bin)
