cmake_minimum_required(VERSION 3.16...3.27)
set(name CaTS)
project(${name} VERSION 0.1.0)

option(WITH_ROOT "Build example with ROOT" ON)
set(WITH_ENGINE "geant4" CACHE STRING
    "Simulation engine: opticks, celeritas, or geant4")
string(TOLOWER "${WITH_ENGINE}" WITH_ENGINE)
set_property(CACHE WITH_ENGINE PROPERTY VALUE "${WITH_ENGINE}")
set(_valid_engines opticks celeritas geant4)
if(NOT WITH_ENGINE IN_LIST _valid_engines)
  message(FATAL_ERROR
    "WITH_ENGINE must be one of: opticks, celeritas, geant4 (got '${WITH_ENGINE}')")
endif()
message(STATUS "WITH_ENGINE: ${WITH_ENGINE}")
option(G4ANALYSIS_USE "Build example with Analysis" ON)

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
# CaTS requires shared libraries
#
if(NOT Geant4_shared_FOUND)
  message(FATAL_ERROR "CaTS must use shared libraries")
endif()

if(WITH_ENGINE STREQUAL "opticks")
    # Add Opticks CMake modules to the search path
    list(APPEND CMAKE_MODULE_PATH "$ENV{OPTICKS_HOME}/cmake/Modules")
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
    set(OPTICKS_PREFIX "$ENV{OPTICKS_PREFIX}")

    # Print debug information about CMAKE_MODULE_PATH and OPTICKS_HOME
    message(STATUS "*************************************************************")
    message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
    message(STATUS "OPTICKS_HOME:  $ENV{OPTICKS_HOME}")
    message(STATUS "*************************************************************")

    include(OpticksBuildOptions)
    find_package(G4CX REQUIRED CONFIG)
    find_package(U4 REQUIRED CONFIG)
    message(STATUS "     G4OPTICKS: ${G4OPTICKS}")
else()
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/modules")
  include(CaTSCXXFlags)
endif()

if(WITH_ENGINE STREQUAL "celeritas")
  find_package(Celeritas REQUIRED)
  message(STATUS "     Celeritas_DIR: ${Celeritas_DIR}")
endif()

if(WITH_ROOT)
  find_package(ROOT REQUIRED)
# ROOT version 6 required
   if(ROOT_FOUND)
     STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
     if(NOT VERSION6MATCH)
       message(FATAL_ERROR "${name} requires ROOT 6")
     endif()
   endif()
# Include ROOT's CMake functions for dictionary generation
#  since root6.20, the file is renamed and included by default, so include
#  only when we find the *old* name
   if(EXISTS "${ROOT_DIR}/modules/RootNewMacros.cmake")
     include("${ROOT_DIR}/modules/RootNewMacros.cmake")
   endif()
endif()

if (WITH_ENGINE STREQUAL "opticks")
    file(GLOB OpRHDRS ${PROJECT_SOURCE_DIR}/include/Opticks/*.hh)
    file(GLOB OpSRCS ${PROJECT_SOURCE_DIR}/src/Opticks/*.cc)
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
if(WITH_ROOT)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
   include_directories(${PROJECT_SOURCE_DIR}/include
     ${Geant4_INCLUDE_DIR}
     ${ROOT_INCLUDE_DIRS})
#----------------------------------------------------------------------------
# Generate dictionaries, add ROOT libraries properties
#
   REFLEX_GENERATE_DICTIONARY(CaTSClasses include/CaTSClasses.hh SELECTION xml/selection.xml)
   add_library(CaTSClassesDict SHARED ${sources} CaTSClasses.cxx)
   set(libsuffix .so)
   set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})
   set_target_properties(CaTSClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
   target_link_libraries(CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
  include_directories(${PROJECT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIR})
endif()

file(GLOB detectors ${PROJECT_SOURCE_DIR}/gdml/*.gdml)
file(GLOB schemas ${PROJECT_SOURCE_DIR}/gdml/*.xsd)
file(GLOB _scripts_all ${PROJECT_SOURCE_DIR}/scripts/*)
set(scripts)
foreach(_s IN LISTS _scripts_all)
    if(NOT IS_DIRECTORY ${_s})
        list(APPEND scripts ${_s})
    endif()
endforeach()
file(GLOB macros ${PROJECT_SOURCE_DIR}/macros/*.mac)
file(GLOB runscripts ${PROJECT_SOURCE_DIR}/scripts/run.sh  ${PROJECT_SOURCE_DIR}/scripts/check.sh)

add_executable(${name} ${name}.cc ${sources} ${headers})
target_compile_definitions(${name} PRIVATE WITH_ROOT)
target_link_libraries(${name} PRIVATE CaTSClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
if(WITH_ENGINE STREQUAL "opticks")
  target_sources(${name} PRIVATE ${OpRHDRS} ${OpSRCS})
  target_compile_definitions(${name} PRIVATE WITH_G4OPTICKS)
  target_link_libraries(${name} PRIVATE Opticks::G4CX Opticks::U4)
elseif(WITH_ENGINE STREQUAL "celeritas")
  target_compile_definitions(${name} PRIVATE WITH_CELERITAS)
  target_link_libraries(${name} PRIVATE Celeritas::celeritas Celeritas::G4)
else()
  target_sources(${name} PRIVATE include/SteppingAction.hh)
endif()
set(install_targets ${name})

if(WITH_ROOT)
   add_executable(readPhotonHits readPhotonHits.cc  ${sources} ${headers})
   target_link_libraries(readPhotonHits CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
   add_executable(readCalorimeterHits  readCalorimeterHits.cc   ${sources} ${headers})
   target_link_libraries(readCalorimeterHits CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
   add_executable(readMscHits readMscHits.cc  ${sources} ${headers})
   target_link_libraries(readMscHits CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
   add_executable(readDRCalorimeterHits  readDRCalorimeterHits.cc   ${sources} ${headers})
   target_link_libraries(readDRCalorimeterHits CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
   add_executable(validateCeleritas validateCeleritas.cc  ${sources} ${headers})
   target_link_libraries(validateCeleritas CaTSClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} )
   link_directories( ${ROOT_LIBRARY_DIR} )
endif()

install(TARGETS ${install_targets} DESTINATION bin)
install(FILES ${detectors} ${schemas}  ${scripts} ${macros}  DESTINATION bin)
install(FILES ${detectors} DESTINATION bin/gdml)
install(FILES ${macros} DESTINATION bin/macros)
install(PROGRAMS ${runscripts}  DESTINATION bin)


if(WITH_ENGINE STREQUAL "opticks")
    # Environment variables for Opticks
    # Sets the maximum number of photons to simulate on the GPU
    # Selects which GPU to use
    # Enables or disables debugging
    install(FILES OptickEnv DESTINATION bin)
endif ()

if(WITH_ROOT)
  install(TARGETS CaTSClassesDict DESTINATION bin)
  install(TARGETS readPhotonHits DESTINATION bin)
  install(TARGETS readCalorimeterHits DESTINATION bin)
  install(TARGETS readDRCalorimeterHits DESTINATION bin)
  install(TARGETS readMscHits DESTINATION bin)
  install(TARGETS validateCeleritas DESTINATION bin)
  install(FILES ${PROJECT_BINARY_DIR}/CaTSClasses_rdict.pcm DESTINATION bin)
endif()

message(STATUS "G4ANALYSIS_USE: ${G4ANALYSIS_USE}")
message(STATUS "WITH_ROOT:      ${WITH_ROOT}")
if(WITH_ROOT)
  message(STATUS "     ROOTSYS: $ENV{ROOTSYS}")
endif(WITH_ROOT)
