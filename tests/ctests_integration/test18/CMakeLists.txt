# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(test18)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test18
  test18.cc
  src/Nucleus.cc
  src/RadioactiveDecayGun.cc
  src/RadioactiveDecayGunmessenger.cc
  src/Tst18ActionInitialization.cc
  src/Tst18EventAction.cc
  src/Tst18EventActionMessenger.cc
  src/Tst18GeometryConstruction.cc
  src/Tst18PhysicsList.cc
  src/Tst18PrimaryGeneratorAction.cc
  src/Tst18RunAction.cc
  src/Tst18RunActionMessenger.cc
  src/Tst18StepCut.cc
  src/Tst18SteppingAction.cc
  src/UIcmdWithNucleusAndUnit.cc)
target_include_directories(test18 PRIVATE include)
target_link_libraries(test18 PRIVATE ${Geant4_LIBRARIES})

#--Test definitions -------
if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test18
    COMMAND test18 ${CMAKE_CURRENT_SOURCE_DIR}/test18.large_N.in
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
