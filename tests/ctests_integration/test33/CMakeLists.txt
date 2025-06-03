cmake_minimum_required(VERSION 3.16...3.21)
project(test33)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test33 test33.cc ${sources} ${headers})
target_include_directories(test33 PRIVATE include)
target_link_libraries(test33 PRIVATE ${Geant4_LIBRARIES})

set(Test33_SCRIPTS
    SamplingChange.mac
    MassGeo_TimedApp.mac
    ParallelGeo_TimedApp.mac
  )

foreach(_script ${Test33_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test33
    COMMAND test33 ${CMAKE_CURRENT_SOURCE_DIR}/MassGeo_TimedApp.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  geant4_add_test(test33-1
    COMMAND test33 ${CMAKE_CURRENT_SOURCE_DIR}/ParallelGeo_TimedApp.mac
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
