#----------------------------------------------------------------------------
# Re-set target output directories to ensure targets are output to per-test
# subdirectories
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

foreach(_c ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER ${_c} _c_upper)
  unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${_c_upper})
  unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${_c_upper})
  unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${_c_upper})
endforeach()

#----------------------------------------------------------------------------
# Add test projects
#
add_subdirectory(test00)
add_subdirectory(test01)
add_subdirectory(test02)
add_subdirectory(test03)
add_subdirectory(test04)
add_subdirectory(test05)
add_subdirectory(test06)
add_subdirectory(test07)
add_subdirectory(test08)
add_subdirectory(test09)

add_subdirectory(test10)
add_subdirectory(test11)
add_subdirectory(test12)
add_subdirectory(test13)
add_subdirectory(test14)
add_subdirectory(test15)
add_subdirectory(test16)
add_subdirectory(test17)
add_subdirectory(test18)
# ? Doesn't seem to be part of CTesting (uses code in test23)...
add_subdirectory(test19)

add_subdirectory(test22)
# ? Doesn't seem to be part of CTesting
add_subdirectory(test23)
add_subdirectory(test24)
add_subdirectory(test25)
add_subdirectory(test27)
add_subdirectory(test28)

add_subdirectory(test30)
add_subdirectory(test31)
add_subdirectory(test32)
add_subdirectory(test320)
add_subdirectory(test33)
add_subdirectory(test34)
add_subdirectory(test35)
add_subdirectory(test36)
add_subdirectory(test37)
add_subdirectory(test38)

add_subdirectory(test40)
add_subdirectory(test41)
add_subdirectory(test44)
add_subdirectory(test45)
add_subdirectory(test46)
# ? Doesn't seem to be part of CTesting (uses code in test23)...
add_subdirectory(test47)
add_subdirectory(test48)

add_subdirectory(test50)
add_subdirectory(test52)
add_subdirectory(test53)
add_subdirectory(test54)
add_subdirectory(test55)
add_subdirectory(test56)
add_subdirectory(test58)

add_subdirectory(test60)
add_subdirectory(test61)
add_subdirectory(test62)
add_subdirectory(test64)
add_subdirectory(test65)
add_subdirectory(test67)
add_subdirectory(test68)
add_subdirectory(test69)

add_subdirectory(test70)
add_subdirectory(test73)
add_subdirectory(test74)
# ? Doesn't seem to be part of CTesting (uses code in test23)...
add_subdirectory(test75)
add_subdirectory(test76)
add_subdirectory(test77)
add_subdirectory(test78)
add_subdirectory(test79)

add_subdirectory(test202)

# Only relevant if building in MT mode
if(GEANT4_BUILD_MULTITHREADED)
  add_subdirectory(test400)
endif()
