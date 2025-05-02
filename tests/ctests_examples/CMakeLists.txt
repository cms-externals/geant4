# - Pre-check for packages that some tests have a hard requirement on
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
find_package(AIDA QUIET)
find_package(VecGeom QUIET)
find_package(Pythia6 QUIET)
find_package(Pythia8 QUIET)

# - Set up variable to check if in BoundsCheck build using CMAKE_CXX_FLAGS
# All code sharing stl containers MUST be compiled with this option, root persistency cannot be used
string(REGEX MATCH "-D_GLIBCXX_DEBUG" BOUNDSCHECK ${CMAKE_CXX_FLAGS})

# Generate custom excludes file
configure_file(CTestCustom.cmake.in ${PROJECT_BINARY_DIR}/CTestCustom.cmake @ONLY)

#-------------------------------------------------------------------------------------------------
#---Example Coding--------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
if(NOT WIN32)
  option(GEANT4_TEST_EXAMPLE_GUIDELINES "Enable tests of examples for coding guideline violations" OFF)
  mark_as_advanced(GEANT4_TEST_EXAMPLE_GUIDELINES)
  if(GEANT4_TEST_EXAMPLE_GUIDELINES)
    geant4_add_test(example-bas-guidelines
      COMMAND ${CMAKE_CURRENT_LIST_DIR}/check_basic_examples.sh --find-inconsistencies
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      LABELS Continuous Nightly)

    geant4_add_test(example-ext-guidelines
      COMMAND ${CMAKE_CURRENT_LIST_DIR}/check_extended_examples.sh --find-inconsistencies
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      LABELS Continuous Nightly)
  endif()
endif()

#-------------------------------------------------------------------------------------------------
#---Basic Examples: Build and Run-----------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
# - Base directories
set(EX_BASIC_SOURCE_DIR ${PROJECT_SOURCE_DIR}/examples/basic)
set(EX_BASIC_BINARY_DIR ${PROJECT_BINARY_DIR}/examples/basic)

geant4_add_test(example-bas-b1
  COMMAND ${EX_BASIC_BINARY_DIR}/B1/exampleB1 ${EX_BASIC_SOURCE_DIR}/B1/exampleB1.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B1
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B1
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b2a
  COMMAND ${EX_BASIC_BINARY_DIR}/B2/B2a/exampleB2a ${EX_BASIC_SOURCE_DIR}/B2/B2a/exampleB2.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B2/B2a
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B2/B2a
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b2b
  COMMAND ${EX_BASIC_BINARY_DIR}/B2/B2b/exampleB2b ${EX_BASIC_SOURCE_DIR}/B2/B2b/exampleB2.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B2/B2b
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B2/B2b
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b3a
  COMMAND ${EX_BASIC_BINARY_DIR}/B3/B3a/exampleB3a ${EX_BASIC_SOURCE_DIR}/B3/B3a/exampleB3.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B3/B3a
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B3/B3a
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b3b
  COMMAND ${EX_BASIC_BINARY_DIR}/B3/B3b/exampleB3b ${EX_BASIC_SOURCE_DIR}/B3/B3b/exampleB3.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B3/B3b
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B3/B3b
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b4a
  COMMAND ${EX_BASIC_BINARY_DIR}/B4/B4a/exampleB4a -m ${EX_BASIC_SOURCE_DIR}/B4/B4a/exampleB4.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B4/B4a
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B4/B4a
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b4b
  COMMAND ${EX_BASIC_BINARY_DIR}/B4/B4b/exampleB4b -m ${EX_BASIC_SOURCE_DIR}/B4/B4b/exampleB4.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B4/B4b
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B4/B4b
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b4c
  COMMAND ${EX_BASIC_BINARY_DIR}/B4/B4c/exampleB4c -m ${EX_BASIC_SOURCE_DIR}/B4/B4c/exampleB4.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B4/B4c
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B4/B4c
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b4d
  COMMAND ${EX_BASIC_BINARY_DIR}/B4/B4d/exampleB4d -m ${EX_BASIC_SOURCE_DIR}/B4/B4d/exampleB4.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B4/B4d
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B4/B4d
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-bas-b5
  COMMAND ${EX_BASIC_BINARY_DIR}/B5/exampleB5 ${EX_BASIC_SOURCE_DIR}/B5/exampleB5.in
  SOURCE_DIR ${EX_BASIC_SOURCE_DIR}/B5
  BINARY_DIR ${EX_BASIC_BINARY_DIR}/B5
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# Check that basic examples build with GNUmake system
include(GNUMakeTests.cmake)

#-------------------------------------------------------------------------------------------------
#---Extended Examples: Build and Run--------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
# - Base directories
set(EX_EXTEND_SOURCE_DIR ${PROJECT_SOURCE_DIR}/examples/extended)
set(EX_EXTEND_BINARY_DIR ${PROJECT_BINARY_DIR}/examples/extended)

foreach(_i 01 03)
  geant4_add_test(example-ext-analysis-anaex${_i}-build
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/analysis/AnaEx${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  foreach(_j csv root xml)
    geant4_add_test(example-ext-analysis-anaex${_i}-${_j}
      COMMAND ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx${_i}/AnaEx${_i}
              ${EX_EXTEND_SOURCE_DIR}/analysis/AnaEx${_i}/AnaEx${_i}-${_j}.in
      WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx${_i}
      DEPENDS example-ext-analysis-anaex${_i}-build
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()

  if(Geant4_hdf5_FOUND)
    geant4_add_test(example-ext-analysis-anaex${_i}-hdf5
      COMMAND ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx${_i}/AnaEx${_i}
              ${EX_EXTEND_SOURCE_DIR}/analysis/AnaEx${_i}/AnaEx${_i}-hdf5.in
      WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx${_i}
      DEPENDS example-ext-analysis-anaex${_i}-build
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()
endforeach()

if(ROOT_FOUND)
  geant4_add_test(example-ext-analysis-anaex02
    COMMAND ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx02/AnaEx02 ${EX_EXTEND_SOURCE_DIR}/analysis/AnaEx02/AnaEx02.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/analysis/AnaEx02
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/analysis/AnaEx02
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-ext-analysis-b1con
  COMMAND ${EX_EXTEND_BINARY_DIR}/analysis/B1Con/exampleB1Con ${EX_EXTEND_SOURCE_DIR}/analysis/B1Con/exampleB1Con.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/analysis/B1Con
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/analysis/B1Con
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i 01 02 03)
  geant4_add_test(example-ext-biasing-b${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/biasing/B${_i}/exampleB${_i}
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/biasing/B${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/biasing/B${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

foreach(_i 01 02 03 04 05 06 07)
  geant4_add_test(example-ext-biasing-gb${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/biasing/GB${_i}/exampleGB${_i}
            ${EX_EXTEND_SOURCE_DIR}/biasing/GB${_i}/exampleGB${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/biasing/GB${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/biasing/GB${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-biasing-reversemc01
  COMMAND ${EX_EXTEND_BINARY_DIR}/biasing/ReverseMC01/exampleRMC01
          ${EX_EXTEND_SOURCE_DIR}/biasing/ReverseMC01/run_adjoint_simulation_electron.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/biasing/ReverseMC01
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/biasing/ReverseMC01
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-common
  COMMAND ${EX_EXTEND_BINARY_DIR}/common/exampleCommon
          ${EX_EXTEND_SOURCE_DIR}/common/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/common
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/common
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18)
  geant4_add_test(example-ext-electromagnetic-testem${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm${_i}/TestEm${_i}
            ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm${_i}/TestEm${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-electromagnetic-testem0-DirectAccess
  COMMAND ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm0/DirectAccess
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm0
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm0
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-electromagnetic-testem3-geom
  COMMAND ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm3/TestEm3
          ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm3/geom.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm3
  DEPENDS example-ext-electromagnetic-testem3
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-electromagnetic-testem5-geom
  COMMAND ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm5/TestEm5
          ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm5/geom.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm5
  DEPENDS example-ext-electromagnetic-testem5
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-electromagnetic-testem5-pixe_ANSTO
  COMMAND ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm5/TestEm5
          ${EX_EXTEND_SOURCE_DIR}/electromagnetic/TestEm5/pixe_ANSTO.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/electromagnetic/TestEm5
  DEPENDS example-ext-electromagnetic-testem5
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-errorpropagation-errProp
  COMMAND ${EX_EXTEND_BINARY_DIR}/errorpropagation/errProp/errProp
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/errorpropagation/errProp
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/errorpropagation/errProp
  BUILD
  ENVIRONMENT
    G4ERROR_TARGET=PLANE_SURFACE
    G4ERROR_MODE=BACKWARDS
    G4ERROR_PROP=UNTIL_TARGET
    ${GEANT4_TEST_ENVIRONMENT} )

geant4_add_test(example-ext-eventgenerator-userPrimaryGenerator
  COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/userPrimaryGenerator/userPrimaryGenerator
          ${EX_EXTEND_SOURCE_DIR}/eventgenerator/userPrimaryGenerator/userPrimaryGenerator.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/userPrimaryGenerator
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/userPrimaryGenerator
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# -----
# - GPS Example Tests
geant4_add_test(example-ext-eventgenerator-exgps
  COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps/exgps
          ${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/exgps_batch.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# - Additional exgps tests
geant4_add_test(gps-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps
  BUILD
  LABELS PhysicsChecks)

set(_wd ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps/stattest)
file(MAKE_DIRECTORY ${_wd})

set(_wd2 ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps/stattest-full)
file(MAKE_DIRECTORY ${_wd2})

configure_file(
	${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/stattest/x-ray-spectrum-geant4-180kV.dat
	${_wd2}/x-ray-spectrum-geant4-180kV.dat
	COPYONLY)

geant4_add_test(gps-run-simple
  COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps/exgps
	        ${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/exgps_MT.in
  WORKING_DIRECTORY ${_wd}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS gps-build
  LABELS PhysicsChecks)

geant4_add_test(gps-run
  COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/exgps/exgps
	        ${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/stattest/run.mac
  WORKING_DIRECTORY ${_wd2}
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  DEPENDS gps-build
  LABELS PhysicsChecks)

find_package(StatTest QUIET)
if(STATTEST_FOUND)
	stattest_add_test(gps-simple-checkOutput
    G4TEST gps-run-simple
    CONFIG "${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/stattest/conf-simple.qa"
    REFERENCE "${GEANT4_TEST_REFERENCES_URL}/examples/extended/eventgenerator/exgps/reference-V10-01-ref-01-simple.root"
    INPUT "${_wd}/exgps_t*.root"
    IMG "${_wd}/gps-simple.pdf"
    LABELS PhysicsChecks)

  stattest_add_test(gps-checkOutput
	  G4TEST gps-run
    CONFIG "${EX_EXTEND_SOURCE_DIR}/eventgenerator/exgps/stattest/conf.qa"
    REFERENCE "${GEANT4_TEST_REFERENCES_URL}/examples/extended/eventgenerator/exgps/reference-V10-01-ref-01.root"
    INPUT "${_wd2}/exgps_t*.root"
    IMG "${_wd2}/gps.pdf"
    LABELS PhysicsChecks)
endif()
# -----

geant4_add_test(example-ext-eventgenerator-particleGun-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/particleGun
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/particleGun
  BUILD)

foreach(_i 1 2 3 4)
  geant4_add_test(example-ext-eventgenerator-particleGun-r${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/particleGun/particleGun
            ${EX_EXTEND_SOURCE_DIR}/eventgenerator/particleGun/run${_i}.mac
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/eventgenerator/particleGun
    DEPENDS example-ext-eventgenerator-particleGun-build
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

if(Pythia6_FOUND)
  geant4_add_test(example-ext-eventgenerator-pythia-decayer6-build
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/pythia/decayer6
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/decayer6
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(example-ext-eventgenerator-pythia-decayer6
    COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/decayer6/pythia6_decayer
            ${EX_EXTEND_SOURCE_DIR}/eventgenerator/pythia/decayer6/pythia6_decayer.in
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/decayer6
    DEPENDS example-ext-eventgenerator-pythia-decayer6-build
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

if(Pythia8_FOUND)
  geant4_add_test(example-ext-eventgenerator-pythia-py8decayer-build
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/eventgenerator/pythia/py8decayer
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/py8decayer
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(example-ext-eventgenerator-pythia-py8decayer
    COMMAND ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/py8decayer/pythia8_decayer
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/eventgenerator/pythia/py8decayer
    DEPENDS example-ext-eventgenerator-pythia-py8decayer-build
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-ext-exoticphysics-monopole
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/monopole/monopole
          -m ${EX_EXTEND_SOURCE_DIR}/exoticphysics/monopole/monopole.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/monopole
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/monopole
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-dmparticle
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/dmparticle/dmparticle
          ${EX_EXTEND_SOURCE_DIR}/exoticphysics/dmparticle/dmparticle.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/dmparticle
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/dmparticle
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# -----
# - Phonon example is relatively complex, needs work
#Check if we are in a MT specific build
set(_MTsuffix "")
if(GEANT4_BUILD_MULTITHREADED)
  # use MATCHES to match also AppleClang
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(_MTsuffix "-clang-MT")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(_MTsuffix "-MT")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(_MTsuffix "-clang-MT")
  elseif(MSVC)
    set(_MTsuffix "-win-MT")
  endif()
endif()
set(_compareScript ${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon/ComparePhonon.cmake)
set(_causticRef ${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon/caustic.out)
set(_timingRef ${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon/timing.out)

#Check if specific OS references are available, if not use default
if(EXISTS ${_causticRef}${_MTsuffix})
  set(_causticRef ${_causticRef}${_MTsuffix})
else()
  message(STATUS "ctests_examples: example-ext-exoticphysics-phonons, cannot find caustic.out${_MTsuffix} references, using default caustic.out")
endif()

if(EXISTS ${_timingRef}${_MTsuffix})
  set(_timingRef ${_timingRef}${_MTsuffix})
else()
  message(STATUS "ctests_examples: example-ext-exoticphysics-phonons, cannot find timing.out${_MTsuffix} references, using default timing.out")
endif()

# Split tests into build, run, then comparison
# 1. Build
geant4_add_test(example-ext-exoticphysics-phonons-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/phonon
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# GF: disable temporarily
if(FALSE)
  # 2. Run
  geant4_add_test(example-ext-exoticphysics-phonons-run
    COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/phonon/XGeBox
      ${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon/run.in
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/exoticphysics/phonon
    DEPENDS example-ext-exoticphysics-phonons-build
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      G4LATTICEDATA=${EX_EXTEND_SOURCE_DIR}/exoticphysics/phonon/CrystalMaps)

  # 3. Compare
  geant4_add_test(example-ext-exoticphysics-phonons-timecomp
    DEPENDS example-ext-exoticphysics-phonons-run
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/exoticphysics/phonon
    COMMAND ${CMAKE_COMMAND} -Dtest_file=timing.ssv -Dreference_file=${_timingRef} -P ${_compareScript})

  geant4_add_test(example-ext-exoticphysics-phonons-causticcomp
    DEPENDS example-ext-exoticphysics-phonons-run
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/exoticphysics/phonon
    COMMAND ${CMAKE_COMMAND} -Dtest_file=caustic.ssv -Dreference_file=${_causticRef} -P ${_compareScript})
endif()
# -----

geant4_add_test(example-ext-exoticphysics-channeling-ch0
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch0/ch0
		      ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch0/2009_PLB680_129.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/channeling/ch0
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch0
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-channeling-ch1
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch1/ch1
	  ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch1/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/channeling/ch1
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch1
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-channeling-ch2
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch2/ch2
	  ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch2/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/channeling/ch2
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-channeling-ch3
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch3/ch3
          ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch3/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/channeling/ch3
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/channeling/ch3
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-saxs-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/saxs
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/saxs
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-saxs
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/saxs/saxs
          ${EX_EXTEND_SOURCE_DIR}/exoticphysics/saxs/saxs.in
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/exoticphysics/saxs
  DEPENDS example-ext-exoticphysics-saxs-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-exoticphysics-saxs-slits
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/saxs/saxs
          ${EX_EXTEND_SOURCE_DIR}/exoticphysics/saxs/saxs_slits.in
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/exoticphysics/saxs
  DEPENDS example-ext-exoticphysics-saxs-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i 01 02 03 05)
  geant4_add_test(example-ext-field-field${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/field/field${_i}/field${_i}
            ${EX_EXTEND_SOURCE_DIR}/field/field${_i}/field${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/field/field${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/field/field${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

foreach(_i 04 06)
  geant4_add_test(example-ext-field-field${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/field/field${_i}/field${_i}
            -m ${EX_EXTEND_SOURCE_DIR}/field/field${_i}/field${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/field/field${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/field/field${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-field-blinetracer
  COMMAND ${EX_EXTEND_BINARY_DIR}/field/BlineTracer/testBlineTracer
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/field/BlineTracer
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/field/BlineTracer
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

if(GEANT4_USE_G3TOG4)
  geant4_add_test(example-ext-g3tog4-clgeometry
    COMMAND ${EX_EXTEND_BINARY_DIR}/g3tog4/clGeometry/clGeometry
            ${EX_EXTEND_SOURCE_DIR}/g3tog4/clGeometry/data/testmodel.dat
            -m ${EX_EXTEND_SOURCE_DIR}/g3tog4/clGeometry/clGeometry.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/g3tog4/clGeometry
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/g3tog4/clGeometry
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-ext-geometry-transforms
  COMMAND ${EX_EXTEND_BINARY_DIR}/geometry/transforms/transforms
          ${EX_EXTEND_SOURCE_DIR}/geometry/transforms/transforms.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/geometry/transforms
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/geometry/transforms
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

if(VecGeom_FOUND AND NOT BOUNDSCHECK)
  geant4_add_test(example-ext-geometry-vecGeomNav
    COMMAND ${EX_EXTEND_BINARY_DIR}/geometry/vecGeomNavigation/exampleVecGeomNav
            ${EX_EXTEND_SOURCE_DIR}/geometry/vecGeomNavigation/vecgeomNav.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/geometry/vecGeomNavigation
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/geometry/vecGeomNavigation
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
endif()


foreach(_i 00 01 02 03 04 05 06 07 08 10)
  geant4_add_test(example-ext-hadronic-hadr${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/Hadr${_i}/Hadr${_i}
            ${EX_EXTEND_SOURCE_DIR}/hadronic/Hadr${_i}/hadr${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/Hadr${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/Hadr${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-hadronic-hadr09
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/Hadr09/Hadr09
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/Hadr09
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/Hadr09
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-FissionFragment
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/FissionFragment/FissionFragment
	  -i ${EX_EXTEND_SOURCE_DIR}/hadronic/FissionFragment/batch.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/FissionFragment
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/FissionFragment
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-FlukaCern-ProcessLevel-CrossSection
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/FlukaCern/ProcessLevel/CrossSection/HadronNucleusXS
	  ${EX_EXTEND_SOURCE_DIR}/hadronic/FlukaCern/ProcessLevel/CrossSection/all_XS.in FTFP_BERT_HP
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/FlukaCern/ProcessLevel/CrossSection
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/FlukaCern/ProcessLevel/CrossSection
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-FlukaCern-ProcessLevel-FinalState
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/FlukaCern/ProcessLevel/FinalState/HadNucIneEvents --test
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/FlukaCern/ProcessLevel/FinalState
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/FlukaCern/ProcessLevel/FinalState
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-NeutronSource
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/NeutronSource/NeutronSource
          ${EX_EXTEND_SOURCE_DIR}/hadronic/NeutronSource/neutronSource.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/NeutronSource
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/NeutronSource
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-ParticleFluence-Sphere
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Sphere/Sphere
          ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Sphere/all_together.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Sphere
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Sphere
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-ParticleFluence-ConcentricSpheres
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/ConcentricSpheres/ConcentricSpheres
          ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/ConcentricSpheres/all_together.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/ConcentricSpheres
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/ConcentricSpheres
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-ParticleFluence-Layer
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Layer/Layer
          ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Layer/all_together.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Layer
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Layer
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-hadronic-ParticleFluence-Calo
  COMMAND ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Calo/Calo
          ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Calo/atlasHec.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/hadronic/ParticleFluence/Calo
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/hadronic/ParticleFluence/Calo
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dicom
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/DICOM/DICOM
          ${EX_EXTEND_SOURCE_DIR}/medical/DICOM/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/DICOM
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/DICOM
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dicom2
  PRECMD ${CMAKE_COMMAND}
         -DDICOM_DIR=${EX_EXTEND_BINARY_DIR}/medical/DICOM
         ${EX_EXTEND_SOURCE_DIR}/medical/DICOM2
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/DICOM2/dicom2
          ${EX_EXTEND_SOURCE_DIR}/medical/DICOM2/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/DICOM2
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/DICOM2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/medical/DICOM2
  DEPENDS example-ext-medical-dicom)

# dicom2 relies on a library built by dicom. This is fine on rpath-based platforms,
# but Windows builds rely on PATH to locate DLLs so we update the environment
# We apply the env change to the test *directly* because geant4_add_test/driver can't
# handle path-like env vars properly
if(WIN32)
  if(CMAKE_GENERATOR MATCHES Makefiles)
    set(dcfg)
  else()
    set(dcfg $<CONFIG>)
  endif()
  # Trick on Windows: PATH separator is semi-colon, but CMake uses this internally as
  # a list separator. Must escape ; *3* times to properly pass it through to CTest environment!
  # - https://stackoverflow.com/questions/59862894/how-do-i-make-a-list-in-cmake-with-the-semicolon-value
  # - https://gitlab.kitware.com/cmake/cmake/-/issues/17946
  file(TO_NATIVE_PATH "${EX_EXTEND_BINARY_DIR}/medical/DICOM/${dcfg}" DICOMPATH)
  string(REPLACE ";" "\\\;" DICOM2PATH "${DICOMPATH};$ENV{PATH}")
  set_property(TEST example-ext-medical-dicom2 APPEND PROPERTY
    ENVIRONMENT PATH=${DICOM2PATH} DICOMPATH="${DICOMPATH}")
endif()

geant4_add_test(example-ext-medical-electronScattering
  PRECMD  ${CMAKE_COMMAND} -E make_directory results
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/electronScattering/electronScattering
          ${EX_EXTEND_SOURCE_DIR}/medical/electronScattering/electronScattering.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/electronScattering
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/electronScattering
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-electronScattering2
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/electronScattering2/electronScattering2
          ${EX_EXTEND_SOURCE_DIR}/medical/electronScattering2/electronScattering2.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/electronScattering2
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/electronScattering2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-radiobiology
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/radiobiology/radiobiology
          ${EX_EXTEND_SOURCE_DIR}/medical/radiobiology/macro/defaultMacro.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/radiobiology
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/radiobiology
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


geant4_add_test(example-ext-medical-fanoCavity
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/fanoCavity/fanoCavity
          ${EX_EXTEND_SOURCE_DIR}/medical/fanoCavity/fanoCavity.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/fanoCavity
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/fanoCavity
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-fanoCavity2
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/fanoCavity2/fanoCavity2
          ${EX_EXTEND_SOURCE_DIR}/medical/fanoCavity2/fanoCavity2.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/fanoCavity2
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/fanoCavity2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-GammaTherapy
  COMMAND ${EX_EXTEND_BINARY_DIR}/medical/GammaTherapy/GammaTherapy
          ${EX_EXTEND_SOURCE_DIR}/medical/GammaTherapy/GammaTherapy.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/GammaTherapy
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/medical/GammaTherapy
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-optical-opnovice-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/optical/OpNovice
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-optical-opnovice-opnovice
  COMMAND ${EX_EXTEND_BINARY_DIR}/optical/OpNovice/OpNovice
          -m ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice/OpNovice.in
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/optical/OpNovice
  DEPENDS example-ext-optical-opnovice-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-optical-opnovice-opgdml
  COMMAND ${EX_EXTEND_BINARY_DIR}/optical/OpNovice/OpNovice
          -g ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice/gdml/NoviceExample.gdml
          -m ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice/OpNovice.in
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/optical/OpNovice
  DEPENDS example-ext-optical-opnovice-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-optical-opnovice-optPhoton
  COMMAND ${EX_EXTEND_BINARY_DIR}/optical/OpNovice/OpNovice
          -m ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice/optPhoton.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/optical/OpNovice
  DEPENDS example-ext-optical-opnovice-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-optical-opnovice2-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice2
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/optical/OpNovice2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i boundary coated complexRindex electron fresnel OpNovice2 scint_by_particle wls)
  geant4_add_test(example-ext-optical-opnovice2-${_i}
   COMMAND ${EX_EXTEND_BINARY_DIR}/optical/OpNovice2/OpNovice2
           ${EX_EXTEND_SOURCE_DIR}/optical/OpNovice2/${_i}.mac
   WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/optical/OpNovice2
   DEPENDS example-ext-optical-opnovice2-build
   ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-optical-lxe-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/optical/LXe
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/optical/LXe
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i cerenkov LXe photon wls)
  geant4_add_test(example-ext-optical-lxe-${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/optical/LXe/LXe
            ${EX_EXTEND_SOURCE_DIR}/optical/LXe/${_i}.mac
    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/optical/LXe
    DEPENDS example-ext-optical-lxe-build
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-optical-wls-electron
  COMMAND ${EX_EXTEND_BINARY_DIR}/optical/wls/wls
          ${EX_EXTEND_SOURCE_DIR}/optical/wls/electron.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/optical/wls
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/optical/wls
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i 1 2 3 a)
  geant4_add_test(example-ext-parameterisations-gflash${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/parameterisations/gflash/gflash${_i}/ExGflash${_i}
            ${EX_EXTEND_BINARY_DIR}/parameterisations/gflash/gflash${_i}/test.mac
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parameterisations/gflash/gflash${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parameterisations/gflash/gflash${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()


geant4_add_test(example-ext-parameterisations-par01
  COMMAND ${EX_EXTEND_BINARY_DIR}/parameterisations/Par01/examplePar01
          ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par01/examplePar01.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par01
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parameterisations/Par01
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-parameterisations-par03
  COMMAND ${EX_EXTEND_BINARY_DIR}/parameterisations/Par03/examplePar03
          -m ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par03/examplePar03.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par03
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parameterisations/Par03
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-parameterisations-par04
  COMMAND ${EX_EXTEND_BINARY_DIR}/parameterisations/Par04/examplePar04
          -m ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par04/examplePar04.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par04
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parameterisations/Par04
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

if(GEANT4_USE_GDML)
  geant4_add_test(example-ext-parameterisations-par02
    COMMAND ${EX_EXTEND_BINARY_DIR}/parameterisations/Par02/examplePar02
            ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par02/examplePar02.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parameterisations/Par02
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parameterisations/Par02
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(example-ext-persistency-gdml-g01
    PRECMD  ${CMAKE_COMMAND} -E remove -f g01.gdml
    COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G01/load_gdml
            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G01/solids.gdml
            g01.gdml
            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G01/g01.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G01
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G01
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

#if(NOT WIN32)
#  # Repeat g01 using local schema only
#  geant4_add_test(example-ext-persistency-gdml-g01-schema
#    PRECMD  ${CMAKE_COMMAND} -E remove -f g01-schema.gdml
#    COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G01/load_gdml
#            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G01/solids.gdml
#            g01-schema.gdml
#            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G01/g01.in
#    WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G01
#    DEPENDS example-ext-persistency-gdml-g01
#    PASSREGEX "G4GDML: Loaded alternative schema URI: ${PROJECT_SOURCE_DIR}/source/persistency/gdml/schema/gdml.xsd"
#    ENVIRONMENT G4GDML_SCHEMA_FILE=${PROJECT_SOURCE_DIR}/source/persistency/gdml/schema/gdml.xsd ${GEANT4_TEST_ENVIRONMENT})
#endif()

  geant4_add_test(example-ext-persistency-gdml-g04
    COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G04/gdml_det
            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G04/auxiliary.gdml
            ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G04/g04.mac
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/persistency/gdml/G04
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/persistency/gdml/G04
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

#P01/02 require ROOT 6
if(ROOT_FOUND AND Geant4_shared_FOUND AND NOT BOUNDSCHECK)
  string(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
  if(VERSION6MATCH)
    geant4_add_test(example-ext-persistency-p01-write
      COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/P01/exampleP01
              ${EX_EXTEND_SOURCE_DIR}/persistency/P01/run.mac
      SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/persistency/P01
      BINARY_DIR ${EX_EXTEND_BINARY_DIR}/persistency/P01
      BUILD
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

    geant4_add_test(example-ext-persistency-p01-read
      COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/P01/readHits hits.root
      SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/persistency/P01
      BINARY_DIR ${EX_EXTEND_BINARY_DIR}/persistency/P01
      BUILD
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS example-ext-persistency-p01-write)
   endif()
endif()

geant4_add_test(example-ext-persistency-p03
  PRECMD  ${CMAKE_COMMAND} -E copy ${EX_EXTEND_SOURCE_DIR}/persistency/P03/g4geom.txt g4geom.txt
  COMMAND ${EX_EXTEND_BINARY_DIR}/persistency/P03/textGeom
          ${EX_EXTEND_SOURCE_DIR}/persistency/P03/batch.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/persistency/P03
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/persistency/P03
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-polarisation-pol01
  COMMAND ${EX_EXTEND_BINARY_DIR}/polarisation/Pol01/Pol01
          ${EX_EXTEND_SOURCE_DIR}/polarisation/Pol01/pol01.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/polarisation/Pol01
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/polarisation/Pol01
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-physicslists-factory
  COMMAND ${EX_EXTEND_BINARY_DIR}/physicslists/factory/factory
          -m ${EX_EXTEND_SOURCE_DIR}/physicslists/factory/run.mac
          -p QGSP_BERT
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/physicslists/factory
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/physicslists/factory
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-physicslists-factory-environment
  COMMAND ${EX_EXTEND_BINARY_DIR}/physicslists/factory/factory
          -m ${EX_EXTEND_SOURCE_DIR}/physicslists/factory/run.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/physicslists/factory
  DEPENDS example-ext-physicslists-factory
  ENVIRONMENT PHYSLIST=FTFP_BERT_EMV ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-physicslists-genericPL
  COMMAND ${EX_EXTEND_BINARY_DIR}/physicslists/genericPL/genericPL
          -m ${EX_EXTEND_SOURCE_DIR}/physicslists/genericPL/run.mac
          -p FTFP_BERT.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/physicslists/genericPL
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/physicslists/genericPL
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-physicslists-extensibleFactory
  COMMAND ${EX_EXTEND_BINARY_DIR}/physicslists/extensibleFactory/extensibleFactory
           -m ${EX_EXTEND_SOURCE_DIR}/physicslists/extensibleFactory/run.mac
           -p MySpecialPhysList
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/physicslists/extensibleFactory
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/physicslists/extensibleFactory
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-radioactivedecay-activation
  COMMAND ${EX_EXTEND_BINARY_DIR}/radioactivedecay/Activation/Activation
          ${EX_EXTEND_SOURCE_DIR}/radioactivedecay/Activation/Activation.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/radioactivedecay/Activation
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/radioactivedecay/Activation
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

foreach(_i 01 02)
  geant4_add_test(example-ext-radioactivedecay-rdecay${_i}
    COMMAND ${EX_EXTEND_BINARY_DIR}/radioactivedecay/rdecay${_i}/rdecay${_i}
            ${EX_EXTEND_SOURCE_DIR}/radioactivedecay/rdecay${_i}/rdecay${_i}.in
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/radioactivedecay/rdecay${_i}
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/radioactivedecay/rdecay${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-runandevent-re01
  PRECMD  ${CMAKE_COMMAND} -E copy ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE01/pythia_event.data pythia_event.data
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE01/RE01
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE01/sample.in
  ERROR sample.err
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE01
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE01
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re02
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02/RE02
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE02/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE02
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re02-run3
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02/RE02
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE02/run3.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02
  DEPENDS example-ext-runandevent-re02
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re02-run4
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02/RE02
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE02/run4.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE02
  DEPENDS example-ext-runandevent-re02
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re03
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE03/RE03
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE03/run.mac
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE03
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE03
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

file(WRITE ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE04/score.ref "G4VScoringMesh::DrawColorChart(): no visualization system\n")
geant4_add_test(example-ext-runandevent-re04
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE04/RE04
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE04/score.mac
  POSTCMD ${CMAKE_COMMAND} -E compare_files score.err score.ref
  ERROR score.err
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE04
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE04
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  LABELS Nightly Continuous)

geant4_add_test(example-ext-runandevent-re05
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE05/exampleRE05
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE05/exampleRE05.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE05
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE05
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re06
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE06/exampleRE06
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE06/exampleRE06.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE06
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE06
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re07-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE07
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07
  BUILD)

geant4_add_test(example-ext-runandevent-re07-processes
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07/exampleRE07
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE07/processes.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07
  DEPENDS example-ext-runandevent-re07-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re07-tracking
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07/exampleRE07
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE07/tracking.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07
  DEPENDS example-ext-runandevent-re07-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-runandevent-re07-specialized
  COMMAND ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07/exampleRE07
          ${EX_EXTEND_SOURCE_DIR}/runAndEvent/RE07/specialized.mac
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/runAndEvent/RE07
  DEPENDS example-ext-runandevent-re07-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


geant4_add_test(example-ext-visualization-uservisaction-build
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/visualization/userVisAction
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/visualization/userVisAction
  BUILD)

geant4_add_test(example-ext-visualization-uservisaction-uva
  COMMAND ${EX_EXTEND_BINARY_DIR}/visualization/userVisAction/userVisAction
          ${EX_EXTEND_SOURCE_DIR}/visualization/userVisAction/userVisAction.in
  WORKING_DIRECTORY ${EX_EXTEND_BINARY_DIR}/visualization/userVisAction
  DEPENDS example-ext-visualization-uservisaction-build
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-visualization-movies
  COMMAND ${EX_EXTEND_BINARY_DIR}/visualization/movies/movies
          -m ${EX_EXTEND_SOURCE_DIR}/visualization/movies/movies.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/visualization/movies
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/visualization/movies
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


geant4_add_test(example-ext-exoticphysics-ExUCN
  COMMAND ${EX_EXTEND_BINARY_DIR}/exoticphysics/ucn/ExUCN
          -m ${EX_EXTEND_SOURCE_DIR}/exoticphysics/ucn/ExUCN.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/exoticphysics/ucn
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/exoticphysics/ucn
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

set(DNA_EXAMPLE_SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/medical/dna)
set(DNA_EXAMPLE_BIN_DIR ${EX_EXTEND_BINARY_DIR}/medical/dna)

geant4_add_test(example-ext-medical-dna-pdb4dna
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/pdb4dna/pdb4dna
          -mac ${DNA_EXAMPLE_SOURCE_DIR}/pdb4dna/pdb4dna.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/pdb4dna
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/pdb4dna
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-wholeNuclearDNA
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/wholeNuclearDNA/wholeNuclearDNA
          -mac ${DNA_EXAMPLE_SOURCE_DIR}/wholeNuclearDNA/wholeNuclearDNA.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/wholeNuclearDNA
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/wholeNuclearDNA
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-dnadamage1
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/dnadamage1/dnadamage1
          ${DNA_EXAMPLE_SOURCE_DIR}/dnadamage1/dnadamage1.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/dnadamage1
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/dnadamage1
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-dnadamage2
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/dnadamage2/dnadamage2
          ${DNA_EXAMPLE_SOURCE_DIR}/dnadamage2/dnadamage2.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/dnadamage2
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/dnadamage2
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-dnaphysics
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/dnaphysics/dnaphysics
          ${DNA_EXAMPLE_SOURCE_DIR}/dnaphysics/dnaphysics.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/dnaphysics
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/dnaphysics
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-microdosimetry
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/microdosimetry/microdosimetry
          ${DNA_EXAMPLE_SOURCE_DIR}/microdosimetry/microdosimetry.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/microdosimetry
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/microdosimetry
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-microyz
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/microyz/microyz
          ${DNA_EXAMPLE_SOURCE_DIR}/microyz/microyz.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/microyz
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/microyz
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-microprox
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/microprox/microprox
          ${DNA_EXAMPLE_SOURCE_DIR}/microprox/microprox.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/microprox
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/microprox
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-radial
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/radial/radial
          ${DNA_EXAMPLE_SOURCE_DIR}/radial/radial.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/radial
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/radial
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-svalue
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/svalue/svalue
          ${DNA_EXAMPLE_SOURCE_DIR}/svalue/svalue.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/svalue
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/svalue
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-splitting
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/splitting/splitting
          ${DNA_EXAMPLE_SOURCE_DIR}/splitting/run.mac
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/splitting
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/splitting
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-wvalue
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/wvalue/wvalue
          ${DNA_EXAMPLE_SOURCE_DIR}/wvalue/wvalue.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/wvalue
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/wvalue
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-slowing
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/slowing/slowing
          ${DNA_EXAMPLE_SOURCE_DIR}/slowing/slowing.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/slowing
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/slowing
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-spower
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/spower/spower
          ${DNA_EXAMPLE_SOURCE_DIR}/spower/spower.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/spower
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/spower
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-neuron
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/neuron/neuron
          -mac ${DNA_EXAMPLE_SOURCE_DIR}/neuron/neuron.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/neuron
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/neuron
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-range
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/range/range
          ${DNA_EXAMPLE_SOURCE_DIR}/range/range.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/range
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/range
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-clustering
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/clustering/clustering -mac
          ${DNA_EXAMPLE_SOURCE_DIR}/clustering/runOneEvent.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/clustering
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/clustering
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-icsd
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/icsd/icsd
          ${DNA_EXAMPLE_SOURCE_DIR}/icsd/icsd.mac
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/icsd
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/icsd
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-mfp
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/mfp/mfp
          ${DNA_EXAMPLE_SOURCE_DIR}/mfp/mfp.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/mfp
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/mfp
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-AuNP
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/AuNP/AuNP
          ${DNA_EXAMPLE_SOURCE_DIR}/AuNP/AuNP.mac
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/AuNP
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/AuNP
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-scavenger
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/scavenger/scavenger
          ${DNA_EXAMPLE_SOURCE_DIR}/scavenger/beam.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/scavenger
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/scavenger
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-jetcounter
        COMMAND ${DNA_EXAMPLE_BIN_DIR}/jetcounter/jetcounter run
        SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/jetcounter
        BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/jetcounter
        BUILD
        ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-ext-medical-dna-UHDR
        COMMAND ${DNA_EXAMPLE_BIN_DIR}/UHDR/UHDR
                ${DNA_EXAMPLE_SOURCE_DIR}/UHDR/beam.in
        SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/UHDR
        BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/UHDR
        BUILD
        ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


foreach(_i 1 2 3 4 5)
  # flag -mt always ON if G4MULTITHREADED is not defined, it will
  # run in sequential mode
  geant4_add_test(example-ext-medical-dna-chem${_i}
    COMMAND ${DNA_EXAMPLE_BIN_DIR}/chem${_i}/chem${_i}
            -mac ${DNA_EXAMPLE_SOURCE_DIR}/chem${_i}/beam.in
    SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/chem${_i}
    BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/chem${_i}
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

geant4_add_test(example-ext-medical-dna-chem6
  COMMAND ${DNA_EXAMPLE_BIN_DIR}/chem6/chem6
          ${DNA_EXAMPLE_SOURCE_DIR}/chem6/beam.in
  SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/chem6
  BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/chem6
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

set(THREADSAFE_SCORERS_ENVIRON ${GEANT4_TEST_ENVIRONMENT})
if(NOT GEANT4_BUILD_MULTITHREADED)
  list(APPEND THREADSAFE_SCORERS_ENVIRON G4FORCE_RUN_MANAGER_TYPE=Serial)
endif()

geant4_add_test(example-ext-parallel-threadsafe-scorers
  COMMAND ${EX_EXTEND_BINARY_DIR}/parallel/ThreadsafeScorers/ts_scorers
          ${EX_EXTEND_SOURCE_DIR}/parallel/ThreadsafeScorers/ts_scorers.in
  SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parallel/ThreadsafeScorers
  BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parallel/ThreadsafeScorers
  BUILD
  ENVIRONMENT ${THREADSAFE_SCORERS_ENVIRON})

# - MPI Examples
find_package(MPI QUIET)
if(MPI_CXX_FOUND)
  set(G4MPI_BINDIR "${EX_EXTEND_BINARY_DIR}/parallel/MPI")
  set(G4MPI_TEST_NAME_STUB "example-ext-parallel-mpi")
  set(G4MPI_BUILD_TEST "${G4MPI_TEST_NAME_STUB}-core" )

  # - Build the MPI library and example
  geant4_add_test(${G4MPI_BUILD_TEST}
    COMMAND ${CMAKE_COMMAND} -E echo "G4mpi build complete"
    SOURCE_DIR ${EX_EXTEND_SOURCE_DIR}/parallel/MPI
    BINARY_DIR ${EX_EXTEND_BINARY_DIR}/parallel/MPI
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    LABELS MPI)

  # - Test exMPI01
  set(G4MPI_EX01_BINDIR ${G4MPI_BINDIR}/exMPI01)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex01-sequential
    COMMAND ${G4MPI_EX01_BINDIR}/exMPI01 run.mac
    WORKING_DIRECTORY ${G4MPI_EX01_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex01
    COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX01_BINDIR}/exMPI01 run.mac
    WORKING_DIRECTORY ${G4MPI_EX01_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)

  # - Test exMPI02
  # Will only be available if ROOT is
  if(ROOT_FOUND)
    set(G4MPI_EX02_BINDIR ${G4MPI_BINDIR}/exMPI02)

    geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex02-sequential
      COMMAND ${G4MPI_EX02_BINDIR}/exMPI02 run.mac
      WORKING_DIRECTORY ${G4MPI_EX02_BINDIR}
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      DEPENDS ${G4MPI_BUILD_TEST}
      LABELS MPI)

    geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex02
      COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX02_BINDIR}/exMPI02 run.mac
      WORKING_DIRECTORY ${G4MPI_EX02_BINDIR}
		  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
		  DEPENDS ${G4MPI_BUILD_TEST}
		  LABELS MPI)
  else()
    message(STATUS "G4 Examples: MPI example exMPI02 disabled: ROOT not found")
  endif()

  # - Test exMPI03
  set(G4MPI_EX03_BINDIR ${G4MPI_BINDIR}/exMPI03)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex03-sequential
    COMMAND ${G4MPI_EX03_BINDIR}/exMPI03 run.mac
    WORKING_DIRECTORY ${G4MPI_EX03_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex03
    COMMAND ${MPIEXEC} -n 2 ${G4MPI_EX03_BINDIR}/exMPI03 run.mac
    WORKING_DIRECTORY ${G4MPI_EX03_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)

  # - Test exMPI04
  set(G4MPI_EX04_BINDIR ${G4MPI_BINDIR}/exMPI04)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex04-sequential
    COMMAND ${G4MPI_EX04_BINDIR}/exMPI04 run.mac
    WORKING_DIRECTORY ${G4MPI_EX04_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)

  geant4_add_test(${G4MPI_TEST_NAME_STUB}-ex04
    COMMAND ${MPIEXEC} -n 3 ${G4MPI_EX04_BINDIR}/exMPI04 run.mac
    WORKING_DIRECTORY ${G4MPI_EX04_BINDIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    DEPENDS ${G4MPI_BUILD_TEST}
    LABELS MPI)
endif()

#-------------------------------------------------------------------------------------------------
#---Advanced Examples: Build and Run--------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
set(EX_ADVANCED_SOURCE_DIR ${PROJECT_SOURCE_DIR}/examples/advanced)
set(EX_ADVANCED_BINARY_DIR ${PROJECT_BINARY_DIR}/examples/advanced)

geant4_add_test(example-adv-STCyclotron
  COMMAND ${EX_ADVANCED_BINARY_DIR}/STCyclotron/STCyclotron
          ${EX_ADVANCED_SOURCE_DIR}/STCyclotron/run.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/STCyclotron
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/STCyclotron
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-air_shower
  COMMAND ${EX_ADVANCED_BINARY_DIR}/air_shower/Ultra
          ${EX_ADVANCED_SOURCE_DIR}/air_shower/UltraMacro.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/air_shower
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/air_shower
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-amsEcal
  COMMAND ${EX_ADVANCED_BINARY_DIR}/amsEcal/amsEcal
          ${EX_ADVANCED_SOURCE_DIR}/amsEcal/run1.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/amsEcal
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/amsEcal
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-brachytherapy
  COMMAND ${EX_ADVANCED_BINARY_DIR}/brachytherapy/Brachy
          ${EX_ADVANCED_SOURCE_DIR}/brachytherapy/FlexiSourceMacro.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/brachytherapy
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/brachytherapy
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-ICRP110Phantoms
  COMMAND ${EX_ADVANCED_BINARY_DIR}/ICRP110_HumanPhantoms/ICRP110phantoms
          ${EX_ADVANCED_SOURCE_DIR}/ICRP110_HumanPhantoms/male_head.in
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/ICRP110_HumanPhantoms
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/ICRP110_HumanPhantoms
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-ICRP145Phantoms
  COMMAND ${EX_ADVANCED_BINARY_DIR}/ICRP145_HumanPhantoms/ICRP145phantoms
          -m ${EX_ADVANCED_SOURCE_DIR}/ICRP145_HumanPhantoms/example.in
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/ICRP145_HumanPhantoms
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/ICRP145_HumanPhantoms
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# ICRP145/ICRP110 appear to consume significant memory resource, so prevent them running
# concurrently with a resource lock.
# TODO: ICRP145 may suffer from a memory leak as it is observed to consume increasing amounts
# of memory during a run. To be investigated further.
#set_tests_properties(example-adv-ICRP110Phantoms example-adv-ICRP145Phantoms
#  PROPERTIES RESOURCE_LOCK G4LargeMemoryUse)

# Temporarily revert to use RUN_SERIAL on ICRP145 to determine if this can run in isolation at all
set_tests_properties(example-adv-ICRP145Phantoms PROPERTIES RUN_SERIAL TRUE)

geant4_add_test(example-adv-fastAerosol
  COMMAND ${EX_ADVANCED_BINARY_DIR}/fastAerosol/fastAerosol
          ${EX_ADVANCED_SOURCE_DIR}/fastAerosol/test.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/fastAerosol
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/fastAerosol
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

# ChargeExchangeMC has a mandatory GDML dependency
if(GEANT4_USE_GDML)
  geant4_add_test(example-adv-ChargeExchangeMC
    COMMAND ${EX_ADVANCED_BINARY_DIR}/ChargeExchangeMC/ChargeExchangeMC
            -p ${EX_ADVANCED_BINARY_DIR}/ChargeExchangeMC/preinit.mac
            -m ${EX_ADVANCED_BINARY_DIR}/ChargeExchangeMC/batch.mac
		SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/ChargeExchangeMC
    BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/ChargeExchangeMC
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )

  geant4_add_test(example-adv-gorad
    COMMAND ${EX_ADVANCED_BINARY_DIR}/gorad/gorad
            ${EX_ADVANCED_SOURCE_DIR}/gorad/test.mac
    SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/gorad
    BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/gorad
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-adv-composite_calorimeter
  COMMAND ${EX_ADVANCED_BINARY_DIR}/composite_calorimeter/CompositeCalorimeter
          ${EX_ADVANCED_SOURCE_DIR}/composite_calorimeter/test.g4mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/composite_calorimeter
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/composite_calorimeter
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
              CCAL_CONFPATH=./dataconf
              CCAL_SENSITIVECONF=g4testbeamhcal96.conf
              CCAL_GEOMETRYCONF=testbeamhcal96.conf
              CCAL_GLOBALPATH=./dataglobal
              CCAL_GEOMPATH=./datageom
              CCAL_VISPATH=./datavis)

geant4_add_test(example-adv-eRosita
  COMMAND ${EX_ADVANCED_BINARY_DIR}/eRosita/eRosita
          ${EX_ADVANCED_SOURCE_DIR}/eRosita/eRosita.in
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/eRosita
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/eRosita
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-gammaknife
  COMMAND ${EX_ADVANCED_BINARY_DIR}/gammaknife/gammaknife
          ${EX_ADVANCED_SOURCE_DIR}/gammaknife/simulation.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/gammaknife
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/gammaknife
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-gammaray_telescope
  COMMAND ${EX_ADVANCED_BINARY_DIR}/gammaray_telescope/GammaRayTel
          ${EX_ADVANCED_SOURCE_DIR}/gammaray_telescope/gammaraytel.in
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/gammaray_telescope
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/gammaray_telescope
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-microbeam
  COMMAND ${EX_ADVANCED_BINARY_DIR}/microbeam/microbeam
          ${EX_ADVANCED_SOURCE_DIR}/microbeam/microbeam.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/microbeam
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/microbeam
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-microelectronics
  COMMAND ${EX_ADVANCED_BINARY_DIR}/microelectronics/microelectronics
          ${EX_ADVANCED_SOURCE_DIR}/microelectronics/microelectronics.in
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/microelectronics
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/microelectronics
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

if(GEANT4_USE_SYSTEM_CLHEP)
  geant4_add_test(example-adv-nanobeam
    COMMAND ${EX_ADVANCED_BINARY_DIR}/nanobeam/nanobeam
            ${EX_ADVANCED_SOURCE_DIR}/nanobeam/default.mac
    SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/nanobeam
    BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/nanobeam
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-adv-hadrontherapy
  COMMAND ${EX_ADVANCED_BINARY_DIR}/hadrontherapy/hadrontherapy
          ${EX_ADVANCED_SOURCE_DIR}/hadrontherapy/batch.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/hadrontherapy
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/hadrontherapy
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-human_phantom
  COMMAND ${EX_ADVANCED_BINARY_DIR}/human_phantom/phantom batch.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/human_phantom
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/human_phantom
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-purging_magnet
  COMMAND ${EX_ADVANCED_BINARY_DIR}/purging_magnet/PurgMag MacroTesting.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/purging_magnet
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/purging_magnet
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-exp_microdosimetry
  COMMAND ${EX_ADVANCED_BINARY_DIR}/exp_microdosimetry/exp_microdosimetry
          ${EX_ADVANCED_SOURCE_DIR}/exp_microdosimetry/run.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/exp_microdosimetry
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/exp_microdosimetry
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-xray_fluorescence
  COMMAND ${EX_ADVANCED_BINARY_DIR}/xray_fluorescence/XrayFluo livermore.mac 1
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/xray_fluorescence
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/xray_fluorescence
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} XRAYDATA=. )

geant4_add_test(example-adv-xray_telescope
  COMMAND ${EX_ADVANCED_BINARY_DIR}/xray_telescope/XrayTel test.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/xray_telescope
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/xray_telescope
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-xray_TESdetector
  COMMAND ${EX_ADVANCED_BINARY_DIR}/xray_TESdetector/XrayTESdetector run01.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/xray_TESdetector
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/xray_TESdetector
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-xray_SiliconPoreOptics
  COMMAND ${EX_ADVANCED_BINARY_DIR}/xray_SiliconPoreOptics/XraySiliconPoreOptics run01.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/xray_SiliconPoreOptics
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/xray_SiliconPoreOptics
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-iort_therapy
  COMMAND ${EX_ADVANCED_BINARY_DIR}/iort_therapy/iort_therapy
			    ${EX_ADVANCED_BINARY_DIR}/iort_therapy/batch.mac
	SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/iort_therapy
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/iort_therapy
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-medical_linac
  COMMAND ${EX_ADVANCED_BINARY_DIR}/medical_linac/medical_linac
          ${EX_ADVANCED_BINARY_DIR}/medical_linac/run.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/medical_linac
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/medical_linac
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


if(ROOT_FOUND)
  geant4_add_test(example-adv-doiPET
    COMMAND ${EX_ADVANCED_BINARY_DIR}/doiPET/doiPET
            ${EX_ADVANCED_BINARY_DIR}/doiPET/run.mac
    SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/doiPET
    BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/doiPET
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()


geant4_add_test(example-adv-underground_physics
	COMMAND ${EX_ADVANCED_BINARY_DIR}/underground_physics/DMX
	        ${EX_ADVANCED_BINARY_DIR}/underground_physics/gamma_1000.mac
	SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/underground_physics
	BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/underground_physics
	BUILD
	ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-underground_physics-alpha
	COMMAND ${EX_ADVANCED_BINARY_DIR}/underground_physics/DMX
	        ${EX_ADVANCED_BINARY_DIR}/underground_physics/alpha.mac
  WORKING_DIRECTORY ${EX_ADVANCED_BINARY_DIR}/underground_physics
	DEPENDS example-adv-underground_physics
	ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-lArCal
  COMMAND ${EX_ADVANCED_BINARY_DIR}/lAr_calorimeter/lArCal
	        ${EX_ADVANCED_BINARY_DIR}/lAr_calorimeter/run1.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/lAr_calorimeter
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}//lAr_calorimeter
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


geant4_add_test(example-adv-HGCal_testbeam_electron
  COMMAND ${EX_ADVANCED_BINARY_DIR}/HGCal_testbeam/HGCal_testbeam
          ${EX_ADVANCED_SOURCE_DIR}/HGCal_testbeam/run.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/HGCal_testbeam
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/HGCal_testbeam
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

if(ROOT_FOUND AND Geant4_shared_FOUND AND NOT BOUNDSCHECK)
  geant4_add_test(example-adv-CaTS
    COMMAND ${EX_ADVANCED_BINARY_DIR}/CaTS/CaTS
    -g ${EX_ADVANCED_SOURCE_DIR}/CaTS/gdml/simpleLArTPC.gdml
    -m  ${EX_ADVANCED_SOURCE_DIR}/CaTS/macros/test_G4.mac
    SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/CaTS
    BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/CaTS
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()

geant4_add_test(example-adv-eFLASH_radiotherapy
  COMMAND ${EX_ADVANCED_BINARY_DIR}/eFLASH_radiotherapy/flash
          ${EX_ADVANCED_SOURCE_DIR}/eFLASH_radiotherapy/run.mac
  SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/eFLASH_radiotherapy
  BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/eFLASH_radiotherapy
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

geant4_add_test(example-adv-stim_pixe_tomography
        COMMAND ${EX_ADVANCED_BINARY_DIR}/stim_pixe_tomography/stim_pixe_tomography -s ${EX_ADVANCED_SOURCE_DIR}/stim_pixe_tomography/test.in
        SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/stim_pixe_tomography
        BINARY_DIR ${EX_ADVANCED_BINARY_DIR}/stim_pixe_tomography
        BUILD
        ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

set(ADV_DNA_EXAMPLE_SOURCE_DIR ${EX_ADVANCED_SOURCE_DIR}/dna)
set(ADV_DNA_EXAMPLE_BIN_DIR ${EX_ADVANCED_BINARY_DIR}/dna)

geant4_add_test(example-adv-dna-moleculardna
  COMMAND ${ADV_DNA_EXAMPLE_BIN_DIR}/moleculardna/molecular
  -m ${ADV_DNA_EXAMPLE_SOURCE_DIR}/moleculardna/cylinders.mac
  SOURCE_DIR ${ADV_DNA_EXAMPLE_SOURCE_DIR}/moleculardna
  BINARY_DIR ${ADV_DNA_EXAMPLE_BIN_DIR}/moleculardna
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})


geant4_add_test(example-adv-dna-dsbandrepair
  COMMAND ${ADV_DNA_EXAMPLE_BIN_DIR}/dsbandrepair/dsbandrepair
          ${ADV_DNA_EXAMPLE_BIN_DIR}/dsbandrepair/dsbandrepair.in
  SOURCE_DIR ${ADV_DNA_EXAMPLE_SOURCE_DIR}/dsbandrepair
  BINARY_DIR ${ADV_DNA_EXAMPLE_BIN_DIR}/dsbandrepair
  BUILD
  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(example-adv-dna-cellularPhantom
    COMMAND ${ADV_DNA_EXAMPLE_BIN_DIR}/cellularPhantom/cellularPhantom
            ${ADV_DNA_EXAMPLE_BIN_DIR}/cellularPhantom/run.mac
    SOURCE_DIR ${ADV_DNA_EXAMPLE_SOURCE_DIR}/cellularPhantom
    BINARY_DIR ${ADV_DNA_EXAMPLE_BIN_DIR}/cellularPhantom
    BUILD
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
