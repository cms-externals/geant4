cmake_minimum_required(VERSION 3.16...3.21)
project(test38)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#   0 = those compiled and added directly to main program
#   1 = separate library, but linked at build time
#   2 = separate library, but loaded only at run time
#
# ComandLineParser taken from test70

set(sources0
  ${PROJECT_SOURCE_DIR}/src/CommandLineParser.cc
  ${PROJECT_SOURCE_DIR}/include/CommandLineParser.hh
  ${PROJECT_SOURCE_DIR}/include/MyPL0.hh
  ${PROJECT_SOURCE_DIR}/include/MyPL0.icc)

set(sources1
  ${PROJECT_SOURCE_DIR}/src/MyStuff1.cc
  ${PROJECT_SOURCE_DIR}/include/MyPL1.hh
  ${PROJECT_SOURCE_DIR}/include/MyPL1.icc)

set(sources2
  ${PROJECT_SOURCE_DIR}/src/MyStuff2.cc
  ${PROJECT_SOURCE_DIR}/src/G4NewDecayPhysics.cc
  ${PROJECT_SOURCE_DIR}/src/G4NewExoticPhysics.cc
  ${PROJECT_SOURCE_DIR}/include/MyPL2.hh
  ${PROJECT_SOURCE_DIR}/include/MyPL2.icc
  ${PROJECT_SOURCE_DIR}/include/MyNSPL3.hh
	${PROJECT_SOURCE_DIR}/include/MyNSPL3.icc
	${PROJECT_SOURCE_DIR}/include/MyNSPL4.hh
	${PROJECT_SOURCE_DIR}/include/MyNSPL4.icc)

#----------------------------------------------------------------------------
# Make libraries
set(T38PRELINK "")
set(T38TRY_EXTRAS 0)

if(Geant4_shared_FOUND)
  set(T38BUILDTYPE "SHARED")

  # for now do NOT try to prelink lib1 on windows as MYGEANT4_LINK_LIBRARY
  # doesn't yet make the appropriate .lib to go along with the .dll
  # when it does this should be 'if(T38TRY_EXTRAS)
  if(NOT WIN32 AND T38TRY_EXTRAS)
    set(T38PRELINK g4plftest38_lib1 ${CMAKE_DL_LIBS})
  else()
    set(T38PRELINK ${CMAKE_DL_LIBS})
    set(T38SPECIALFLAG "-DG4TEST38_WIN_NO_PRELINK1")
  endif()

else()
  set(T38BUILDTYPE "STATIC")

  if(T38TRY_EXTRAS)
    set(T38PRELINK g4plftest38_lib1 g4plftest38_lib2)
  endif()
endif()

if(T38TRY_EXTRAS)
  add_library(g4plftest38_lib1 ${T38BUILDTYPE} ${sources1})
  target_include_directories(g4plftest38_lib1 PUBLIC include)
  target_link_libraries(g4plftest38_lib1 PUBLIC ${Geant4_LIBRARIES})
  set_target_properties(g4plftest38_lib1 PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)

  add_library(g4plftest38_lib2 ${T38BUILDTYPE} ${sources2})
  target_include_directories(g4plftest38_lib2 PUBLIC include)
  target_link_libraries(g4plftest38_lib2 PUBLIC ${Geant4_LIBRARIES})
  set_target_properties(g4plftest38_lib2 PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# - Build the main exe
add_executable(g4plfactory_test38 g4plfactory_test38.cc ${sources0})
target_include_directories(g4plfactory_test38 PRIVATE include)
target_compile_definitions(g4plfactory_test38
  PRIVATE
    G4TEST38_${T38BUILDTYPE}=1
    $<IF:$<BOOL:${T38TRY_EXTRAS}>,G4TEST38_TRY_EXTRAS=1,G4TEST38_SKIP_EXTRAS=1>
    ${T38SPECIALFLAG}
)
target_link_libraries(g4plfactory_test38 PRIVATE ${T38PRELINK} ${Geant4_LIBRARIES})

#---Test definition-----------------------------------------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test38
    COMMAND g4plfactory_test38 --defaults -f -F -c -r 2 -V 0 -v 0
	  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
	  LABELS Nightly Continuous PhysicsChecks)
endif()