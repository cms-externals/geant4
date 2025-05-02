cmake_minimum_required(VERSION 3.16...3.21)
project(test67)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(ROOT QUIET)
if(ROOT_FOUND)
  # copy .root files to build area
  file(GLOB _refs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.root test67*.in)
  file(COPY ${_refs} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else()
  message(STATUS "G4 TESTS: ROOT package not found. --> test67 Analysis disabled")
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(test67 test67.cc ${sources} ${headers})
target_compile_definitions(test67 PRIVATE $<$<BOOL:${ROOT_FOUND}>:G4_USE_ROOT>)
target_include_directories(test67 PRIVATE include $<$<BOOL:${ROOT_FOUND}>:${ROOT_INCLUDE_DIRS}>)
target_link_libraries(test67 PRIVATE ${Geant4_LIBRARIES} "$<$<BOOL:${ROOT_FOUND}>:${ROOT_LIBRARIES}>")

if(GEANT4_ENABLE_TESTING)
	#---Test definitions----------------------------------------------------------------------------------
	geant4_add_test(test67
	  COMMAND test67 ${CMAKE_CURRENT_SOURCE_DIR}/test67.in
	  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
		LABELS Nightly Continuous PhysicsChecks)

	#Add special tests for PhysicsList group
	if(ROOT_FOUND)
	  #Define a test that builds the executable: needed independent since tests after will be
	  #executed in parallel. This ensure correct dependency treatment
	  #Loop on all different test types
	  foreach(_confs test67_livermore test67_stdopt3 test67_stdopt0)
		  #Use for testname same string as macro name (test names need to be unique)
		  geant4_add_test(${_confs}
			  COMMAND test67 ${CMAKE_CURRENT_SOURCE_DIR}/${_confs}.in
			  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
			  LABELS PhysicsChecks)
	   endforeach()
	endif()
endif()
