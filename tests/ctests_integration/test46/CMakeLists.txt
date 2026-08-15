cmake_minimum_required(VERSION 3.16...3.21)
project(test46)
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test46 test46.cc ${sources} ${headers})
target_include_directories(test46 PRIVATE include)
target_link_libraries(test46 PRIVATE ${Geant4_LIBRARIES})

configure_file(
  ${PROJECT_SOURCE_DIR}/vis.mac
  ${PROJECT_BINARY_DIR}/vis.mac
  COPYONLY)

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test46
    COMMAND test46 ${CMAKE_CURRENT_SOURCE_DIR}/test46.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
