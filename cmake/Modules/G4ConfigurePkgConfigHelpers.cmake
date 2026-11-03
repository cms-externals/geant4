#.rst:
# G4ConfigurePkgConfigHelpers
# ---------------------------
#
# This module configures and installs pkg-config scripts and the
# geant4-config program to help clients compile and link against
# the Geant4 libraries.
#
# The geant4-config script provides an sh based interface to provide
# information on the Geant4 installation, including installation prefix,
# version number, compiler and linker flags.
#
# The script is generated from a template file and then installed to the
# known bindir as an executable.
#
# Paths are always hardcoded in the build tree version as this is never
# intended to be relocatable.
# The Install Tree script uses self-location based on that in
# {root,clhep}-config is the install itself is relocatable, otherwise
# absolute paths are encoded.
#
#

#-----------------------------------------------------------------
# License and Disclaimer
#
# The  Geant4 software  is  copyright of the Copyright Holders  of
# the Geant4 Collaboration.  It is provided  under  the terms  and
# conditions of the Geant4 Software License,  included in the file
# LICENSE and available at  http://cern.ch/geant4/license .  These
# include a list of copyright holders.
#
# Neither the authors of this software system, nor their employing
# institutes,nor the agencies providing financial support for this
# work  make  any representation or  warranty, express or implied,
# regarding  this  software system or assume any liability for its
# use.  Please see the license in the file  LICENSE  and URL above
# for the full disclaimer and the limitation of liability.
#
# This  code  implementation is the result of  the  scientific and
# technical work of the GEANT4 collaboration.
# By using,  copying,  modifying or  distributing the software (or
# any work based  on the software)  you  agree  to acknowledge its
# use  in  resulting  scientific  publications,  and indicate your
# acceptance of all terms of the Geant4 Software license.
#
#-----------------------------------------------------------------


#-----------------------------------------------------------------------
# function get_system_include_dirs
#          return list of directories our C++ compiler searches
#          by default.
#
#          The idea comes from CMake's inbuilt technique to do this
#          for the Eclipse and CodeBlocks generators, but we implement
#          our own function because the CMake functionality is internal
#          so we can't rely on it.
function(get_system_include_dirs _dirs)
  # Only for GCC, Clang and Intel
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES GNU OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Clang OR "${CMAKE_CXX_COMPILER_ID}" MATCHES Intel)
    # Proceed
    file(WRITE "${PROJECT_BINARY_DIR}/CMakeFiles/g4dummy" "\n")

    # Save locale, them to "C" english locale so we can parse in English
    set(_orig_lc_all      $ENV{LC_ALL})
    set(_orig_lc_messages $ENV{LC_MESSAGES})
    set(_orig_lang        $ENV{LANG})

    set(ENV{LC_ALL}      C)
    set(ENV{LC_MESSAGES} C)
    set(ENV{LANG}        C)

    execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ARG1} -v -E -x c++ -dD g4dummy
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/CMakeFiles
      ERROR_VARIABLE _cxxOutput
      OUTPUT_VARIABLE _cxxStdout
      )

    file(REMOVE "${PROJECT_BINARY_DIR}/CMakeFiles/g4dummy")

    # Parse and extract search dirs
    set(_resultIncludeDirs )
    if( "${_cxxOutput}" MATCHES "> search starts here[^\n]+\n *(.+ *\n) *End of (search) list" )
      string(REGEX MATCHALL "[^\n]+\n" _includeLines "${CMAKE_MATCH_1}")
      foreach(nextLine ${_includeLines})
        string(REGEX REPLACE "\\(framework directory\\)" "" nextLineNoFramework "${nextLine}")
        string(STRIP "${nextLineNoFramework}" _includePath)
        list(APPEND _resultIncludeDirs "${_includePath}")
      endforeach()
    endif()

    # Restore original locale
    set(ENV{LC_ALL}      ${_orig_lc_all})
    set(ENV{LC_MESSAGES} ${_orig_lc_messages})
    set(ENV{LANG}        ${_orig_lang})

    set(${_dirs} ${_resultIncludeDirs} PARENT_SCOPE)
  else()
    set(${_dirs} "" PARENT_SCOPE)
  endif()
endfunction()

# Convert CMake library lists to shell linker flags. Preserve full paths
# so that versioned libraries and macOS SDK stubs select the same library
# as the Geant4 build. CMake's LOCATION_<CONFIG> property selects the configured
# location of imported targets, including Release/Debug mappings.
function(_g4_configure_library_flags _output)
  set(_flags)
  set(_selection general)
  foreach(_library IN LISTS ARGN)
    if(_library MATCHES "^(debug|optimized|general)$")
      set(_selection "${_library}")
      continue()
    endif()
    if((_selection STREQUAL "debug" AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
        OR (_selection STREQUAL "optimized" AND CMAKE_BUILD_TYPE STREQUAL "Debug"))
      set(_selection general)
      continue()
    endif()
    set(_selection general)

    if(TARGET "${_library}")
      get_target_property(_type "${_library}" TYPE)
      if(_type STREQUAL "INTERFACE_LIBRARY")
        get_target_property(_interface "${_library}" INTERFACE_LINK_LIBRARIES)
        if(_interface)
          _g4_configure_library_flags(_interface_flags ${_interface})
          list(APPEND _flags "${_interface_flags}")
        endif()
        continue()
      endif()
      string(TOUPPER "${CMAKE_BUILD_TYPE}" _config)
      if(WIN32 AND _type STREQUAL "SHARED_LIBRARY")
        # DLLs are runtime artifacts. Link clients against the import
        # library, respecting the package's configuration mappings.
        get_property(_mapped TARGET "${_library}"
          PROPERTY MAP_IMPORTED_CONFIG_${_config} SET)
        if(_mapped)
          get_target_property(_configs "${_library}" MAP_IMPORTED_CONFIG_${_config})
        else()
          get_target_property(_configs "${_library}" IMPORTED_CONFIGURATIONS)
          if(NOT _configs)
            set(_configs "")
          endif()
          set(_configs "${_config};;${_configs}")
        endif()
        set(_import_library "")
        if(_mapped AND _configs STREQUAL "")
          get_target_property(_import_library "${_library}" IMPORTED_IMPLIB)
        endif()
        foreach(_import_config IN LISTS _configs)
          set(_property IMPORTED_IMPLIB)
          if(NOT _import_config STREQUAL "")
            string(TOUPPER "${_import_config}" _import_config)
            string(APPEND _property "_${_import_config}")
          endif()
          get_target_property(_import_library "${_library}" ${_property})
          if(_import_library)
            break()
          endif()
        endforeach()
        if(NOT _import_library)
          message(FATAL_ERROR "No import library found for ${_library} (${CMAKE_BUILD_TYPE})")
        endif()
        set(_library "${_import_library}")
      elseif(CMAKE_BUILD_TYPE)
        get_target_property(_library "${_library}" LOCATION_${_config})
      else()
        get_target_property(_library "${_library}" LOCATION)
      endif()
    endif()
    if(_library MATCHES "^(.*)/([^/]+)\\.framework(/.*)?$")
      list(APPEND _flags "-F${CMAKE_MATCH_1} -framework ${CMAKE_MATCH_2}")
    elseif(IS_ABSOLUTE "${_library}" OR _library MATCHES "^[A-Za-z]:[/\\\\]"
        OR _library MATCHES "^-")
      list(APPEND _flags "${_library}")
    else()
      list(APPEND _flags "-l${_library}")
    endif()
  endforeach()
  string(JOIN " " _flags ${_flags})
  set(${_output} "${_flags}" PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------
# Only create script if we have a global library build...
#
#if(NOT GEANT4_BUILD_GRANULAR_LIBS AND UNIX)
if(NOT GEANT4_BUILD_GRANULAR_LIBS)
  # Get implicit search paths
  get_system_include_dirs(_cxx_compiler_dirs)

  # Setup variables needed for expansion in configuration file
  # - C++ Filesystem, if needed
  if(GEANT4_CXX_FILESYSTEM_LIBRARY)
    set(G4_LINK_CXX_FILESYSTEM "-l${GEANT4_CXX_FILESYSTEM_LIBRARY}")
  endif()

  # - Static libs
  if(BUILD_STATIC_LIBS)
    set(G4_BUILTWITH_STATICLIBS "yes")
  else()
    set(G4_BUILTWITH_STATICLIBS "no")
  endif()

  # - Multithreading
  if(GEANT4_BUILD_MULTITHREADED)
    set(G4_BUILTWITH_MULTITHREADING "yes")
  else()
    set(G4_BUILTWITH_MULTITHREADING "no")
  endif()

  # - PHP_AS_HP
  if(GEANT4_BUILD_PHP_AS_HP)
    set(G4_BUILTWITH_PHP_AS_HP "yes")
  else()
    set(G4_BUILTWITH_PHP_AS_HP "no")
  endif()

  # - Smart Stack
  if(GEANT4_USE_SMARTSTACK)
    set(G4_BUILTWITH_SMARTSTACK "yes")
  else()
    set(G4_BUILTWITH_SMARTSTACK "no")
  endif()

  # - CLHEP
  if(GEANT4_USE_SYSTEM_CLHEP)
    set(G4_BUILTWITH_CLHEP "no")
    #inc path
    get_filename_component(G4_SYSTEM_CLHEP_INCLUDE_DIR "${CLHEP_INCLUDE_DIR}" ABSOLUTE)

    _g4_configure_library_flags(G4_SYSTEM_CLHEP_LIBRARIES ${CLHEP_LIBRARIES})
  else()
    set(G4_BUILTWITH_CLHEP "yes")
  endif()

  # - EXPAT
  if(GEANT4_USE_SYSTEM_EXPAT)
    set(G4_BUILTWITH_EXPAT "no")
  else()
    set(G4_BUILTWITH_EXPAT "yes")
  endif()

  # - ZLIB
  if(GEANT4_USE_SYSTEM_ZLIB)
    set(G4_BUILTWITH_ZLIB "no")
  else()
    set(G4_BUILTWITH_ZLIB "yes")
  endif()

  # - PTL
  if(GEANT4_USE_SYSTEM_PTL)
    set(G4_BUILTWITH_PTL "no")
  else()
    set(G4_BUILTWITH_PTL "yes")
  endif()

  # - GDML
  if(GEANT4_USE_GDML)
    set(G4_BUILTWITH_GDML "yes")
    set(G4_XERCESC_INCLUDE_DIRS ${XercesC_INCLUDE_DIR})
    list(REMOVE_DUPLICATES G4_XERCESC_INCLUDE_DIRS)
    if(_cxx_compiler_dirs)
      list(REMOVE_ITEM G4_XERCESC_INCLUDE_DIRS ${_cxx_compiler_dirs})
    endif()

    set(G4_XERCESC_CFLAGS )
    foreach(_dir ${G4_XERCESC_INCLUDE_DIRS})
      set(G4_XERCESC_CFLAGS "${G4_XERCESC_CFLAGS} -I${_dir}")
    endforeach()
  else()
    set(G4_BUILTWITH_GDML "no")
  endif()

  # - G3ToG4
  if(GEANT4_USE_G3TOG4)
    set(G4_BUILTWITH_G3TOG4 "yes")
  else()
    set(G4_BUILTWITH_G3TOG4 "no")
  endif()

  # - USolids
  if(GEANT4_USE_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
    set(G4_BUILTWITH_USOLIDS "yes")
    get_target_property(G4_USOLIDS_INCLUDE_DIRS VecGeom::vecgeom INTERFACE_INCLUDE_DIRECTORIES)
    list(REMOVE_DUPLICATES G4_USOLIDS_INCLUDE_DIRS)
    foreach(_dir ${G4_USOLIDS_INCLUDE_DIRS})
      set(G4_USOLIDS_CFLAGS "${G4_USOLIDS_CFLAGS} -I${_dir}")
    endforeach()
    # NB: should ALSO account for VecGeom having compile_options, but
    # this is better handled through proper pkg-config support
  else()
    set(G4_BUILTWITH_USOLIDS "no")
  endif()

  # - Freetype
  if(GEANT4_USE_FREETYPE)
    set(G4_BUILTWITH_FREETYPE "yes")
  else()
    set(G4_BUILTWITH_FREETYPE "no")
  endif()

  # - HDF5
  if(GEANT4_USE_HDF5)
    set(G4_BUILTWITH_HDF5 "yes")
  else()
    set(G4_BUILTWITH_HDF5 "no")
  endif()

  # - Qt
  if(GEANT4_USE_QT)
    set(G4_BUILTWITH_QT "yes")
    set(_qtcomps OpenGL Widgets Gui Core)
    if(QT_VERSION_MAJOR GREATER 5)
      list(PREPEND _qtcomps OpenGLWidgets)
    endif()

    set(G4_QT_INCLUDE_DIRS )
    foreach(_qtc ${_qtcomps})
      list(APPEND G4_QT_INCLUDE_DIRS ${Qt${QT_VERSION_MAJOR}${_qtc}_INCLUDE_DIRS})
    endforeach()

    list(REMOVE_DUPLICATES G4_QT_INCLUDE_DIRS)
    if(_cxx_compiler_dirs)
      list(REMOVE_ITEM G4_QT_INCLUDE_DIRS ${_cxx_compiler_dirs})
    endif()

    set(G4_QT_CFLAGS )
    foreach(_dir ${G4_QT_INCLUDE_DIRS})
      set(G4_QT_CFLAGS "${G4_QT_CFLAGS} -I${_dir}")
    endforeach()

    if(APPLE AND G4_QT_CFLAGS MATCHES "QtCore\.framework")
      set(G4_QT_CFLAGS "${G4_QT_CFLAGS} -F${G4QTLIBPATH}")
    endif()
  else()
    set(G4_BUILTWITH_QT "no")
  endif()

  # - Motif
  if(GEANT4_USE_XM)
    set(G4_BUILTWITH_MOTIF "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_MOTIF "no")
  endif()

  # - RayTracerX
  if(GEANT4_USE_RAYTRACER_X11)
    set(G4_BUILTWITH_RAYTRACERX11 "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_RAYTRACERX11 "no")
  endif()

  # - OpenGL X11
  if(GEANT4_USE_OPENGL_X11)
    set(G4_BUILTWITH_OPENGLX11 "yes")
    set(G4_CONFIG_NEEDS_X11 TRUE)
  else()
    set(G4_BUILTWITH_OPENGLX11 "no")
  endif()

  # - OpenInventor
  if(GEANT4_USE_INVENTOR)
    set(G4_BUILTWITH_INVENTOR "yes")
  else()
    set(G4_BUILTWITH_INVENTOR "no")
  endif()

  # - VTK
  if(GEANT4_USE_VTK)
    set(G4_BUILTWITH_VTK "yes")
  else()
    set(G4_BUILTWITH_VTK "no")
  endif()

  # If we have a module that uses X11, We have to play with the X11
  # paths to get a clean set suitable for inclusion
  if(G4_CONFIG_NEEDS_X11)
    set(_raw_x11_includes ${X11_INCLUDE_DIR})
    list(REMOVE_DUPLICATES _raw_x11_includes)
    if(_cxx_compiler_dirs)
      list(REMOVE_ITEM _raw_x11_includes ${_cxx_compiler_dirs})
    endif()
    set(G4_X11_CFLAGS )
    foreach(_p ${_raw_x11_includes})
      set(G4_X11_CFLAGS "-I${_p} ${G4_X11_CFLAGS}")
    endforeach()
  endif()

  # External dependencies must also be available to shell-based clients.
  # Keep GUI dependencies separate for the --without-gui queries and put
  # their libraries after the Geant4 libraries which depend on them.
  set(_external_includes)
  set(_external_libraries)
  set(_gui_includes)
  set(_gui_libraries)
  if(GEANT4_USE_SYSTEM_EXPAT)
    list(APPEND _external_includes ${EXPAT_INCLUDE_DIRS})
    list(APPEND _external_libraries ${EXPAT_LIBRARIES})
  endif()
  if(GEANT4_USE_SYSTEM_ZLIB)
    list(APPEND _external_includes ${ZLIB_INCLUDE_DIRS})
    list(APPEND _external_libraries ${ZLIB_LIBRARIES})
  endif()
  if(GEANT4_USE_SYSTEM_PTL)
    get_target_property(_ptl_includes PTL::ptl INTERFACE_INCLUDE_DIRECTORIES)
    if(_ptl_includes)
      list(APPEND _external_includes ${_ptl_includes})
    endif()
    list(APPEND _external_libraries PTL::ptl)
  endif()
  if(GEANT4_USE_TBB)
    list(APPEND _external_includes ${TBB_INCLUDE_DIRS})
    list(APPEND _external_libraries TBB::tbb)
  endif()
  if(GEANT4_USE_GDML)
    list(APPEND _external_libraries ${XercesC_LIBRARIES})
  endif()
  if(GEANT4_USE_USOLIDS OR GEANT4_USE_PARTIAL_USOLIDS)
    list(APPEND _external_libraries VecGeom::vecgeom)
  endif()
  if(GEANT4_USE_HDF5)
    list(APPEND _external_includes ${HDF5_INCLUDE_DIRS})
    list(APPEND _external_libraries ${HDF5_C_LIBRARIES})
  endif()
  if(GEANT4_USE_FREETYPE)
    list(APPEND _external_includes ${FREETYPE_INCLUDE_DIRS})
    list(APPEND _external_libraries ${FREETYPE_LIBRARIES})
  endif()
  if(GEANT4_USE_VTK)
    list(APPEND _gui_includes ${VTK_INCLUDE_DIRS})
    list(APPEND _gui_libraries ${VTK_LIBRARIES})
  endif()
  if(GEANT4_USE_INVENTOR)
    if(GEANT4_USE_INVENTOR_QT)
      list(APPEND _gui_libraries SoQt::SoQt)
    elseif(GEANT4_USE_INVENTOR_XT)
      list(APPEND _gui_libraries SoXt::SoXt)
    elseif(GEANT4_USE_INVENTOR_WIN)
      list(APPEND _gui_libraries SoWin::SoWin)
    endif()
    list(APPEND _gui_libraries Coin::Coin)
  endif()
  if(GEANT4_USE_QT)
    foreach(_qtc ${_qtcomps})
      list(APPEND _gui_libraries Qt${QT_VERSION_MAJOR}::${_qtc})
    endforeach()
  endif()
  if(GEANT4_USE_XM)
    list(APPEND _gui_includes ${MOTIF_INCLUDE_DIR})
    list(APPEND _gui_libraries ${MOTIF_LIBRARIES})
  endif()
  if(GEANT4_USE_INVENTOR OR GEANT4_USE_QT OR GEANT4_USE_OPENGL_X11
      OR GEANT4_USE_XM OR GEANT4_USE_OPENGL_WIN32)
    list(APPEND _gui_includes ${OPENGL_INCLUDE_DIR})
    list(APPEND _gui_libraries ${OPENGL_LIBRARIES})
    if(APPLE AND G4_CONFIG_NEEDS_X11)
      list(APPEND _gui_includes ${XQuartzGL_INCLUDE_DIR})
      list(APPEND _gui_libraries ${XQuartzGL_gl_LIBRARY})
    endif()
  endif()
  if(G4_CONFIG_NEEDS_X11)
    list(APPEND _gui_libraries ${X11_Xmu_LIB} ${X11_Xt_LIB}
      ${X11_Xext_LIB} ${X11_SM_LIB} ${X11_ICE_LIB} ${X11_X11_LIB})
  endif()

  # Config-mode packages such as VTK and Coin expose includes on their
  # imported targets rather than in the older <Package>_INCLUDE_DIRS vars.
  foreach(_library IN LISTS _gui_libraries)
    if(TARGET "${_library}")
      get_target_property(_includes "${_library}" INTERFACE_INCLUDE_DIRECTORIES)
      if(_includes)
        list(APPEND _gui_includes ${_includes})
      endif()
    endif()
  endforeach()

  foreach(_kind external gui)
    string(TOUPPER "${_kind}" _upper_kind)
    list(REMOVE_DUPLICATES _${_kind}_includes)
    if(_cxx_compiler_dirs)
      list(REMOVE_ITEM _${_kind}_includes ${_cxx_compiler_dirs})
    endif()
    set(G4_${_upper_kind}_CFLAGS)
    foreach(_dir IN LISTS _${_kind}_includes)
      # Imported packages may wrap their includes in BUILD_INTERFACE.
      # Both scripts describe clients of this configured dependency tree.
      string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+)>" "\\1" _dir "${_dir}")
      string(GENEX_STRIP "${_dir}" _dir)
      if(_dir)
        string(APPEND G4_${_upper_kind}_CFLAGS " -I${_dir}")
      endif()
    endforeach()
    _g4_configure_library_flags(G4_${_upper_kind}_LIBS ${_${_kind}_libraries})
  endforeach()

  # Configure the script
  # - BUILD TREE
  # Ouch, the include path will be LONG, but at least we always have
  # absolute paths...
  set(GEANT4_CONFIG_SELF_LOCATION "# BUILD TREE IS NON-RELOCATABLE")
  set(GEANT4_CONFIG_INSTALL_PREFIX "${PROJECT_BINARY_DIR}")
  set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
  # NB: this only works for *single* mode generators. With multimode
  # generators, which mode to use is not clear...
  set(GEANT4_CONFIG_LIBDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

  get_property(__geant4_buildtree_include_dirs GLOBAL PROPERTY
    GEANT4_BUILDTREE_INCLUDE_DIRS)

  foreach(_dir ${__geant4_buildtree_include_dirs})
    set(GEANT4_CONFIG_INCLUDE_DIRS "${GEANT4_CONFIG_INCLUDE_DIRS} \\
    ${_dir}")
  endforeach()

  # - Data
  geant4_export_datasets(BUILD GEANT4_CONFIG_DATASET_DESCRIPTIONS)

  # Configure the build tree script
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4-config.in
    ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
    @ONLY
    )

  file(COPY
    ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/geant4-config
    DESTINATION ${PROJECT_BINARY_DIR}
    FILE_PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    )

  # - Install Tree
  # Much easier :-)
  # Non-Relocatable case...
  if(CMAKE_INSTALL_IS_NONRELOCATABLE)
    # Hardcoded paths
    set(GEANT4_CONFIG_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
    set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
    set(GEANT4_CONFIG_LIBDIR "${CMAKE_INSTALL_FULL_LIBDIR}")
    set(GEANT4_CONFIG_INCLUDE_DIRS "${CMAKE_INSTALL_FULL_INCLUDEDIR}/Geant4")
  else()
    # Calculate base of self contained install based on relative path from
    # CMAKE_INSTALL_FULL_BINDIR to CMAKE_INSTALL_PREFIX.
    file(RELATIVE_PATH _bin_to_prefix ${CMAKE_INSTALL_FULL_BINDIR} ${CMAKE_INSTALL_PREFIX})
    # Strip any trailing path separators just for neatness.
    string(REGEX REPLACE "[/\\]$" "" _bin_to_prefix "${_bin_to_prefix}")

    set(GEANT4_CONFIG_INSTALL_PREFIX "$scriptloc/${_bin_to_prefix}")
    set(GEANT4_CONFIG_INSTALL_EXECPREFIX \"\")
    set(GEANT4_CONFIG_LIBDIR "\${prefix}/${CMAKE_INSTALL_LIBDIR}")
    set(GEANT4_CONFIG_INCLUDE_DIRS "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}/Geant4")
  endif()

  # - Data
  geant4_export_datasets(INSTALL GEANT4_CONFIG_DATASET_DESCRIPTIONS)

  # Configure the install tree script
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4-config.in
    ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config
    @ONLY
    )

  # Install it
  install(FILES ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    COMPONENT Development
    )

  # Win32 helper file geant4-config.cmd
  if(WIN32)
    # No configuration just a copy
    configure_file(
      ${PROJECT_SOURCE_DIR}/cmake/Templates/geant4-config-cmd.in
      ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config.cmd
      @ONLY
      )

    # Install helper
    install(FILES ${PROJECT_BINARY_DIR}/InstallTreeFiles/geant4-config.cmd
      DESTINATION ${CMAKE_INSTALL_BINDIR}
      PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
      COMPONENT Development
      )
endif()
endif()
