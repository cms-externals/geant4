cmake_minimum_required(VERSION 3.16...3.21)
project(test73)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. --> test73 disabled")
  return()
endif()

# Locate sources and headers for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test73 test73.cc ${sources} ${headers})
target_compile_definitions(test73 PRIVATE G4ANALYSIS_USEROOT)
target_include_directories(test73 PRIVATE include ${ROOT_INCLUDE_DIRS})
target_link_libraries(test73 PRIVATE ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#---Test definition-----------------------------------------
#Prepare running are:
# 1- create working directories
# 2- Copy macros in the directories
# 3- Define tests (since we run in subdirs we can run in parallel)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/local ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT)
file(GLOB MACROS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.mac xmacro.C)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/local)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX)
file(COPY ${MACROS} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT)

#The following tests are enabled only in PhysicsChecks since they
#Will fail for physics reasons: the MSC model is not adeguate for pi-
#When new G4 developments are such that they will not fail anymore
#they can be promoted to nightlies adding the Nichtly and Continous
#LABELS

if(GEANT4_ENABLE_TESTING)
  #geant4_add_test(test73-local
  #  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test73 run_all.mac local
  #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/local
  #    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
  #    LABELS PhysicsChecks)

  geant4_add_test(test73-FTFP_BERT_EMX
    COMMAND test73 run_all.mac FTFP_BERT_EMX
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT_EMX
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      LABELS PhysicsChecks)

  geant4_add_test(test73-FTFP_BERT
    COMMAND test73 run_all.mac FTFP_BERT
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/FTFP_BERT
      ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
      LABELS PhysicsChecks)
endif()
