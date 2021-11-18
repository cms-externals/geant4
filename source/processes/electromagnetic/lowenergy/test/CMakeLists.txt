#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project($ENV{G4TARGET})

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executables, and link it to the Geant4 libraries
#
add_executable($ENV{G4TARGET} $ENV{G4TARGET}.cc ${sources} ${headers})
target_link_libraries($ENV{G4TARGET} ${Geant4_LIBRARIES} )

install(TARGETS $ENV{G4TARGET} DESTINATION bin)

