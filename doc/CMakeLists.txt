#-----------------------------------------------------------------------
# - Experimental building Doxygen docs for Geant4 toolkit code
#-----------------------------------------------------------------------
# Option done here as we probably won't distribute docs until reasonably
# complete.
option(GEANT4_BUILD_DOCS "Build Doxygen documentation for toolkit code" OFF)
geant4_add_feature(GEANT4_BUILD_DOCS "Building Doxygen documenation for toolkit")

if(NOT GEANT4_BUILD_DOCS)
  return()
endif()

find_package(Doxygen 1.10 REQUIRED)

#-----------------------------------------------------------------------
# Should be done at a higher level (we currently have to maintain versions
# outside of project() due to nested calls, see top level CMake script)
set(PROJECT_VERSION ${${PROJECT_NAME}_VERSION})

#-----------------------------------------------------------------------
# Doxygen settings
#-----------------------------------------------------------------------
# - Organised as per https://www.doxygen.nl/manual/config.html
# Project related options
set(DOXYGEN_FULL_PATH_NAMES NO)
#set(DOXYGEN_PROJECT_LOGO "the current geant4 logo in suitable form")

# Allow autobrief for JavaDoc (/**) and Qt (/*!, //!) comments to
# reduce need to use special commands like @brief and help readability
set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
set(DOXYGEN_QT_AUTOBRIEF YES)

# Build related configuration options
set(DOXYGEN_SHOW_USED_FILES NO)

# Configuration options related to the input files
set(DOXYGEN_FILE_PATTERNS
  "*.hh"
  "*.icc"
  "*.md"
  "*.dox"
)
set(DOXYGEN_EXTENSION_MAPPING "icc=C++")
set(DOXYGEN_EXCLUDE_PATTERNS "*/test/*")
#set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")

# Configuration options related to the HTML output
# - DOXYGEN_HTML_... options
# - Can use for things like MathJax, custom stylesheet, search engines etc

# Turn off all other LaTeX generator, all others are NO by default
# - May want to investigate XML generator for Breathe use
set(DOXYGEN_GENERATE_LATEX NO)

# Only extract what is documented - to be reviewed
set(DOXYGEN_EXTRACT_ALL NO)
set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
set(DOXYGEN_SOURCE_BROWSER YES)

set(DOXYGEN_ALIASES
  [[check="@xrefitem check \"CHECK\" \"To be checked\""]]
)
set(DOXYGEN_VERBATIM_VARS DOXYGEN_ALIASES)

# Using a custom Mainpage layout
set(DOXYGEN_LAYOUT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/CustomLayout.xml)

#-----------------------------------------------------------------------
# Call CMake function from FindDoxygen to generate in build dir
#-----------------------------------------------------------------------
doxygen_add_docs(doc
  "${CMAKE_CURRENT_SOURCE_DIR}/groups.dox"
  "${PROJECT_SOURCE_DIR}/source/global/management"
  "${PROJECT_SOURCE_DIR}/source/geometry/management"
  "${PROJECT_SOURCE_DIR}/source/geometry/biasing"
  "${PROJECT_SOURCE_DIR}/source/geometry/divisions"
  "${PROJECT_SOURCE_DIR}/source/geometry/magneticfield"
  "${PROJECT_SOURCE_DIR}/source/geometry/navigation"
  "${PROJECT_SOURCE_DIR}/source/geometry/solids/boolean"
  "${PROJECT_SOURCE_DIR}/source/geometry/solids/csg"
  "${PROJECT_SOURCE_DIR}/source/geometry/solids/specific"
  "${PROJECT_SOURCE_DIR}/source/geometry/volumes"
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Processing C++ documentation with Doxygen"
)
