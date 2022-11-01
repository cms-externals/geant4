
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.12...3.20)

if(G4TARGET)

  if(NOT DEFINED G4TARGET)
    message(FATAL_ERROR "G4TARGET variable is NOT defined.  It must be set using cmake -DG4TARGET=...")
  else()  
    message(STATUS  "Building project '" ${G4TARGET} "' in geometry/Solids/CSG/test")
  endif()

# Option 1 - create unit tests one by one
  project(${G4TARGET})

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
# Setup include directory for this project
#
  include(${Geant4_USE_FILE})
  include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
  add_executable(${G4TARGET} ${G4TARGET}.cc)
  target_link_libraries(${G4TARGET} ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
  install(TARGETS ${G4TARGET} DESTINATION bin)

else()
# Variant for creating all unit tests
#  Enabled by  GEANT4_BUILD_TESTS  ( tbc )
  message(STATUS "Configuring unit tests in geometry/ navigation/ tests")

  if(GEANT4_BUILD_TESTS)

    GEANT4_ADD_UNIT_TESTS(test*.cc
                      INCLUDE_DIRS geometry/management/include
                                   geometry/volumes/include
                                   geometry/navigation/include
                                   geometry/solids/CSG/include
                                   geometry/solids/Boolean/include
                                   geometry/solids/specific/include
                                   geometry/magneticfield/include
                                   global/management/include
                                   global/HEPRandom/include
                                   global/HEPNumerics/include
                                   global/HEPGeometry/include
                                   graphics_reps/include
                                   run/include
                                   particles/management/include
                                   particles/bosons/include
                                   particles/leptons/include
                                   particles/hadrons/mesons/include
                                   particles/hadrons/barions/include
                                   particles/hadrons/ions/include
                                   materials/include
                                   tracking/include
                                   track/include
                                   intercoms/include
                                   interfaces/basic/include
                                   interfaces/common/include
                                   ${CLHEP_INCLUDE_DIRS}
                      LIBRARIES G4geometry G4particles G4interfaces G4run)
  else()
    message(STATUS "Tests for geometry/ navigation/  are not built.")
  endif()		    
endif()