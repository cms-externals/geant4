cmake_minimum_required(VERSION 3.16...3.21)
project(test202)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test202 test202.cc ${sources} ${headers})
target_include_directories(test202 PRIVATE include)
target_link_libraries(test202 PRIVATE ${Geant4_LIBRARIES})

configure_file(testVisSceneAdd.mac testVisSceneAdd.mac COPYONLY) # test expects file to be local
configure_file(testVisSet.mac testVisSet.mac COPYONLY) # test expects file to be local
configure_file(vis.mac vis.mac COPYONLY) # test expects file to be local
configure_file(run0.mac run0.mac COPYONLY) # this file local too
configure_file(OGLIX.mac OGLIX.mac COPYONLY) # this file local too
configure_file(xvfb.sh xvfb.sh COPYONLY) # this file local too

# OGLIX-reference directory
configure_file(OGLIX-reference/visTest-DrawVolume_0.eps OGLIX-reference/visTest-DrawVolume_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-VisAdd_0.eps OGLIX-reference/visTest-VisAdd_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-VisSet_0.eps OGLIX-reference/visTest-VisSet_0.eps  COPYONLY) # this file local too
configure_file(OGLIX-reference/visTest-BeamOn_0.eps OGLIX-reference/visTest-BeamOn_0.eps  COPYONLY) # this file local too

if(GEANT4_ENABLE_TESTING)
  find_package(Xvfb QUIET)

  if(XVFB_FOUND)
    message(STATUS "G4 TESTS: found Xvfb --> run test202")

    geant4_add_test(test202-OGLIX
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/xvfb.sh $<TARGET_FILE:test202> OGLIX
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      TIMEOUT 1500
      ERROR expected_errors.txt)
  endif()
endif()
