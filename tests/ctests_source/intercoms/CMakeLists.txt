#-----------------------------------------------------------------------------#
# - Setup
#-----------------------------------------------------------------------------#
geant4_setup_unit_tests(LINK_LIBRARIES G4intercoms)

#-----------------------------------------------------------------------------#
# - Tests
#-----------------------------------------------------------------------------#
geant4_add_unit_test(G4UIparameter.test.cc)
geant4_add_unit_test(G4UIcommand.test.cc)
geant4_add_unit_test(G4GenericMessenger.test.cc)
