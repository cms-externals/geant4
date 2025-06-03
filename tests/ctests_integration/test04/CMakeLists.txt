cmake_minimum_required(VERSION 3.16...3.21)
project(test04)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test04
  test04.cc
  src/Tst04ActionInitialization.cc
  src/Tst04PrimaryGeneratorAction.cc
  src/Tst04DetectorConstruction.cc
  src/Tst04RunAction.cc
  src/Tst04DetectorMessenger.cc
  src/Tst04SteppingAction.cc
  src/Tst04PhysicsList.cc
  src/Tst04TargetParameterisation.cc
  )
target_include_directories(test04 PRIVATE include)
target_link_libraries(test04 PRIVATE ${Geant4_LIBRARIES})

add_executable(test04_hadronic_exerciser test04_hadronic_exerciser.cc)
target_link_libraries(test04_hadronic_exerciser PRIVATE ${Geant4_LIBRARIES} )

#-----------------------------------------------------------------------
# Check if user provides data locally, otherwise
# download the required dataset and update testing environment to use it
#
if(NOT ${G4PARTICLEHPDATA})
  set (TENDL_VERSION 1.3.2)

  if(EXISTS ${GEANT4_INSTALL_DATADIR}/G4TENDL${TENDL_VERSION})
    set(G4PARTICLEHPDATA ${GEANT4_INSTALL_DATADIR}/G4TENDL${TENDL_VERSION})
  else()
    include(ExternalProject)
    ExternalProject_Add(test04-data
      PREFIX Data/tmp
      SOURCE_DIR Data/g4tendl
      URL http://cern.ch/geant4/support/source/G4TENDL.${TENDL_VERSION}.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND "")
   set(G4PARTICLEHPDATA ${PROJECT_BINARY_DIR}/Data/g4tendl)
  endif()

  list(APPEND GEANT4_TEST_ENVIRONMENT "G4PARTICLEHPDATA=${G4PARTICLEHPDATA}")
endif()

#---Test definitions----------------------------------------------------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test04_hadronic_exerciser
    COMMAND test04_hadronic_exerciser
    OUTPUT  test04.in)
  geant4_add_test(test04
    COMMAND test04 test04.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS test04_hadronic_exerciser
    TIMEOUT 2700)
endif()
