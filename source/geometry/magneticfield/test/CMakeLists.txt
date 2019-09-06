cmake_minimum_required(VERSION 3.1)

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

   add_subdirectory(field07)
#add_subdirectory(RkTestRun)

   find_package(Boost COMPONENTS) 

   if(Boost_FOUND)
       add_subdirectory(boost) 
   endif()

   add_subdirectory(NTST)

   if(EXISTS ${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/test/TestInterpolation)  
     add_subdirectory(TestInterpolation)
   endif()

endif()   