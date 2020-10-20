#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emlowenergy
# Package: Geant4.src.G4processes.G4electromagnetic.G4emlowenergy
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4emdna-moltypes
    HEADERS
        G4Electron_aq.hh
	G4FakeMolecule.hh
        G4H2.hh
        G4H2O2.hh
        G4H2O.hh
        G4H3O.hh
	G4HO2.hh
        G4Hydrogen.hh
	G4O2.hh
	G4O3.hh
        G4OH.hh
	G4Oxygen.hh
        G4DNAMolecule.hh
    SOURCES
        G4Electron_aq.cc
	G4FakeMolecule.cc
        G4H2.cc
        G4H2O2.cc
        G4H2O.cc
        G4H3O.cc
	G4HO2.cc
        G4Hydrogen.cc
	G4O2.cc
	G4O3.cc
        G4OH.cc
	G4Oxygen.cc
        G4DNAMolecule.cc
    GRANULAR_DEPENDENCIES
        G4geometrymng
        G4globman
        G4heprandom
        G4materials
        G4partman
        G4track
        G4emdna-man
        G4emdna-molman
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

