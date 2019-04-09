cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CLHEP_RNG_Tests C CXX)

FIND_PACKAGE(Geant4 REQUIRED)
INCLUDE(${Geant4_USE_FILE})

SET(RNG_TYPES HepJamesRandom MixMaxRng RanluxEngine
              Ranlux64Engine DualRand RanecuEngine
              MTwistEngine RanshiEngine)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_LIST_DIR})

# used for where the binary will end up
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CONFIG)

SET(G4RNG_LIBS )
SET(G4RNG_LIB_MATCH G4clhep G4global)
FOREACH(_LIB ${Geant4_LIBRARIES})
    FOREACH(_MATCH ${G4RNG_LIB_MATCH})
        IF("${_LIB}" MATCHES "${_MATCH}")
            LIST(APPEND G4RNG_LIBS ${_LIB})
        ENDIF()
    ENDFOREACH()
ENDFOREACH()

MESSAGE(STATUS "testRngMT: linking ${G4RNG_LIBS}")

FOREACH(RNG ${RNG_TYPES})

    # create the executable
    ADD_EXECUTABLE(${RNG} ${RNG}.cc RngUnitTest.hh)

    # link to Geant4 libs
    TARGET_LINK_LIBRARIES(${RNG} ${G4RNG_LIBS})

    # add the dependencies
    ADD_DEPENDENCIES(${RNG} ${G4RNG_LIBS})

    # add the test
    INSTALL(TARGETS ${RNG} DESTINATION bin)

ENDFOREACH(RNG ${RNG_TYPES})
