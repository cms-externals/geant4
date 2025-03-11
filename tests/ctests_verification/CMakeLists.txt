# - Only process if we have a `verification` subdirectory in the project dir
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/verification")
  return()
endif()

# - Pre-check for packages that some tests have a hard requirement on
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})
find_package(ROOT QUIET)

#-------------------------------------------------------------------------------------------------
#---Simplified Calorimeter and other PhysicsChecks tests----------------------------------------------------
#-------------------------------------------------------------------------------------------------
include(${PROJECT_SOURCE_DIR}/verification/SimplifiedCalorimeter/physval/CTestDefinitions.txt OPTIONAL)

# - Build and run tests from verification
message("CTESTS FOR verification Processsing")

if(EXISTS ${PROJECT_SOURCE_DIR}/verification/electromagnetic)
  set(VERIFICATION_EM
    attenuation
	  EmSamplingCalo
	  FluctTest
	  MscHanson
	  MSCL3
    MscTest1
	  TestDEDX2
	  tileatlas
	  Tracker
	  zmumu)

    foreach(name IN LISTS VERIFICATION_EM)
      geant4_add_test(em-testing-suite-${name}
        COMMAND ${PROJECT_BINARY_DIR}/verification/electromagnetic/${name}/${name}
                ${PROJECT_SOURCE_DIR}/verification/electromagnetic/${name}/${name}.in
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/verification/electromagnetic/${name}
        BINARY_DIR ${PROJECT_BINARY_DIR}/verification/electromagnetic/${name}
        BUILD ${name}
        ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
    endforeach()
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/verification/G4Med_suite)
  set(G4MED_TEST
    CCCStest
    )
    foreach(name IN LISTS G4MED_TEST)
      geant4_add_test(g4med-suite-${name}
        COMMAND ${PROJECT_BINARY_DIR}/verification/G4Med_suite/${name}/${name}
                ${PROJECT_SOURCE_DIR}/verification/G4Med_suite/${name}/${name}.in
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/verification/G4Med_suite/${name}
        BINARY_DIR ${PROJECT_BINARY_DIR}/verification/G4Med_suite/${name}
        BUILD ${name}
        ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
    endforeach()
endif()

if(EXISTS ${PROJECT_SOURCE_DIR}/verification/hadronic)
  if(ROOT_FOUND)
    geant4_add_test(verification-hadronic-FragTest
      COMMAND ${PROJECT_BINARY_DIR}/verification/hadronic/FragTest/FragTest
              ${PROJECT_SOURCE_DIR}/verification/hadronic/FragTest/FragTest.in
      SOURCE_DIR ${PROJECT_SOURCE_DIR}/verification/hadronic/FragTest
      BINARY_DIR ${PROJECT_BINARY_DIR}/verification/hadronic/FragTest
      BUILD FragTest
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  else()
      message(STATUS "G4 VERIFICATION TESTS: ROOT package not found --> FragTest disabled ")
  endif()
endif()
