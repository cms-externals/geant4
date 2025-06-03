cmake_minimum_required(VERSION 3.16...3.27 FATAL_ERROR)

#---Common Geant4 CTest script----------------------------------------------
include(${CTEST_SCRIPT_DIRECTORY}/g4common.cmake)

#---Addional CTest settings-------------------------------------------------

#---Clean or not clean the Binary-------------------------------------------
if(NOT EXISTS ${CTEST_BINARY_DIRECTORY}/is_docker)
  ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})
endif()

#---CTest commands----------------------------------------------------------
ctest_start("Release" TRACK "Release")

ctest_update(SOURCE ${CTEST_SOURCE_DIRECTORY})
ctest_configure(BUILD ${CTEST_BINARY_DIRECTORY} 
  SOURCE ${CTEST_SOURCE_DIRECTORY}
  OPTIONS "${CTEST_CONFIG_OPTIONS}")
ctest_submit(PARTS Update Configure)

ctest_read_custom_files("${CTEST_BINARY_DIRECTORY}")
ctest_build(BUILD ${CTEST_BINARY_DIRECTORY})
ctest_submit(PARTS Build)

ctest_test(PARALLEL_LEVEL ${ncpu} INCLUDE_LABEL "Nightly")
ctest_submit(PARTS Test Coverage MemCheck Notes ExtraFiles)
