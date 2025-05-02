cmake_minimum_required(VERSION 3.16...3.21)
project(test06)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

if(Geant4_builtin_zlib_FOUND)
	if(HAVE_OFF64_T)
	   add_definitions(-D_LARGEFILE64_SOURCE=1 -D_FILE_OFFSET_BITS=64)
	endif()
  
  # Intel/Clang may warn about -Wdeprecated-non-prototype, but per https://github.com/madler/zlib/issues/633
  # we suppress this warning if the compiler supports the flag
  include(CheckCCompilerFlag)
  check_c_compiler_flag("-Wno-deprecated-non-prototype" TEST06_NEEDS_DNP)
  
	set(G4ZLIB_TARGET "$<IF:$<BOOL:${Geant4_shared_FOUND}>,Geant4::G4zlib,Geant4::G4zlib-static>")
	add_executable(test06-zlib-exam example.c)
  target_link_libraries(test06-zlib-exam PRIVATE ${G4ZLIB_TARGET})
  add_executable(test06-zlib-mgz minigzip.c)
  target_link_libraries(test06-zlib-mgz PRIVATE ${G4ZLIB_TARGET})

  if(TEST06_NEEDS_DNP)
    target_compile_options(test06-zlib-exam PRIVATE "-Wno-deprecated-non-prototype")
    target_compile_options(test06-zlib-mgz PRIVATE "-Wno-deprecated-non-prototype")
  endif()

  #---Test definitions----------------------------------------------------------------------------------
  if(GEANT4_ENABLE_TESTING)
    geant4_add_test(test06-zlib-exam COMMAND test06-zlib-exam)
    geant4_add_test(test06-zlib-mgz
      COMMAND test06-zlib-mgz -d foo.gz
      DEPENDS test06-zlib-exam
      POSTCMD  ${CMAKE_COMMAND} -E compare_files foo ${PROJECT_SOURCE_DIR}/mgz.txt)
  endif()
endif()
