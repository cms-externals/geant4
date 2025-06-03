cmake_minimum_required(VERSION 3.16...3.21)
project(test68)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(Python3 QUIET COMPONENTS Interpreter)
if(NOT Python3_FOUND)
  message(STATUS "Python3 interpreter not found !")
  return()
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test68 test68.cc ${sources} ${headers})
target_include_directories(test68 PRIVATE include)
target_link_libraries(test68 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)

  Set(Use_Timeout 3600)
  set(Use_Labels_Selective Nightly Continuous)
  if(WIN32)
     set(Use_Timeout 7200)
     set(Use_Labels_Selective Nightly)
     message(" On WIndows: increase TIMEOUT to ${Use_Timeout}")
  endif()
      
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/start.rndm ${CMAKE_CURRENT_BINARY_DIR}/start.rndm COPYONLY)

  # NB: These DO NOT work at present because the python scripts assume the exe in the path
  # Need updating to take the full exe path as the argument

  geant4_add_test(test68-reproducibility-k0lFeSci
	  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_k0lFeSci.py $<CONFIG>
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    TIMEOUT ${Use_Timeout}
	  LABELS Nightly Continuous)

  geant4_add_test(test68-reproducibility-pimCuLAr
	  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_pimCuLAr.py $<CONFIG>
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    TIMEOUT ${Use_Timeout}
	  LABELS ${Use_Labels_Selective})

  geant4_add_test(test68-reproducibility-nPbSci
	  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_nPbSci.py $<CONFIG>
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    TIMEOUT ${Use_Timeout}
	  LABELS ${Use_Labels_Selective})

  geant4_add_test(test68-reproducibility-pWLAr
	  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_pWLAr.py $<CONFIG>
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    TIMEOUT ${Use_Timeout}
	  LABELS ${Use_Labels_Selective})

  geant4_add_test(test68-reproducibility-kmPbWO4
	  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/reproducibility_kmPbWO4.py $<CONFIG>
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    TIMEOUT ${Use_Timeout}
	  LABELS ${Use_Labels_Selective})
endif()
