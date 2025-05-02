# Category redefinition file for use in testing
# - Builds each module into its own library 
# - This checks that each module is self contained and is not picking up dependencies
#   by virtue of grouping with other modules
# - Cannot check for internal transitive dependencies, i.e. module has a direct
#   dependence on Foo, but picks up usage requirements for Foo via Bar, which
#   exposes Foo as a public dependency
message(STATUS "Recreated Geant4 categories with one library per module")
geant4_get_modules(__allmods)
foreach(__mod ${__allmods})
  geant4_add_category(${__mod} MODULES ${__mod})
endforeach()