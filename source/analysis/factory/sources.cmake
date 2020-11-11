#------------------------------------------------------------------------------
# Module : G4analysisfac
# Package: Geant4.src.G4analysis.G4analysisfac
#------------------------------------------------------------------------------

# Optional additional links
if(GEANT4_USE_HDF5)
  set(G4analysisfac_G4hdf5 G4hdf5)
endif()

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4analysisfac
  HEADERS
    g4analysis.hh
    g4analysis_defs.hh
    G4GenericAnalysisManager.hh
    G4GenericAnalysisManager.icc
    G4GenericFileManager.hh
    G4GenericFileManager.icc
  SOURCES
    g4analysis.cc
    G4GenericAnalysisManager.cc
    G4GenericFileManager.cc
  GRANULAR_DEPENDENCIES
    G4globman
    G4intercoms
    G4analysismng
    G4hntools
    G4csv
    G4root
    G4xml
    ${G4analysisfac_G4hdf5}
  GLOBAL_DEPENDENCIES
    G4global
    G4intercoms
)

# List any source specific properties here
