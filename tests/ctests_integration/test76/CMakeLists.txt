cmake_minimum_required(VERSION 3.16...3.21)
project(test76)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(test76 test76.cc include/testing.hh)
target_include_directories(test76 PRIVATE include)
target_link_libraries(test76 PRIVATE ${Geant4_LIBRARIES})

if(GEANT4_ENABLE_TESTING)
  geant4_add_test(test76-statistics
    COMMAND test76
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
    LABELS Nightly Continuous
    TIMEOUT 600) # 10 minutes
endif()
