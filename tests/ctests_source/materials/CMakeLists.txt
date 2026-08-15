#-----------------------------------------------------------------------------#
# - Setup
#-----------------------------------------------------------------------------#
geant4_setup_unit_tests(LINK_LIBRARIES G4materials)

#-----------------------------------------------------------------------------#
# - Tests
#-----------------------------------------------------------------------------#
geant4_add_unit_test(G4AtomicShells.test.cc)
geant4_add_unit_test(G4ExtDEDXTable.test.cc)
geant4_add_unit_test(G4Material.test.cc)
geant4_add_unit_test(G4MaterialComponents.test.cc)
geant4_add_unit_test(G4MaterialPropertiesTable.test.cc)
geant4_add_unit_test(G4NistBuilders.test.cc)
geant4_add_unit_test(G4NistManager.test.cc)
geant4_add_unit_test(G4NistManagerBuilder.test.cc)
geant4_add_unit_test(G4SandiaTable.test.cc)

# Tests that use obsolete classes, but might still be relevant
#geant4_add_unit_test(obsolete_testDEDXTables.cc)
#geant4_add_unit_test(obsolete_testMaterialPropertiesTable.cc)
#geant4_add_unit_test(obsolete_testMaterialPropertyVector.cc)

