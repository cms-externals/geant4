cmake_minimum_required(VERSION 3.16...3.21)
project(test64)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test64 test64.cc ${sources} ${headers})
target_include_directories(test64 PRIVATE include)
target_link_libraries(test64 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  # Bertini tests
  geant4_add_test(test64-bert
    COMMAND test64 ${CMAKE_CURRENT_SOURCE_DIR}/bert.in
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  	LABELS Nightly Continuous)

  # Fritiof tests
  geant4_add_test(test64-ftfp
    COMMAND test64 ${CMAKE_CURRENT_SOURCE_DIR}/ftfp.in
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  	LABELS Nightly Continuous)
endif()
