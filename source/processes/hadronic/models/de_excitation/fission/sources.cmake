#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hadronic_deex_fission
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_fission
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
GEANT4_DEFINE_MODULE(NAME G4hadronic_deex_fission
    HEADERS
        G4CompetitiveFission.hh
        G4EvaporationLevelDensityParameter.hh
        G4FissionBarrier.hh
        G4FissionLevelDensityParameter.hh
        G4FissionLevelDensityParameterINCLXX.hh
        G4FissionParameters.hh
        G4FissionProbability.hh
        G4ParaFissionModel.hh
        G4VFissionBarrier.hh
    SOURCES
        G4CompetitiveFission.cc
        G4EvaporationLevelDensityParameter.cc
        G4FissionBarrier.cc
        G4FissionLevelDensityParameter.cc
        G4FissionLevelDensityParameterINCLXX.cc
        G4FissionParameters.cc
        G4FissionProbability.cc
        G4VFissionBarrier.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4globman
        G4had_mod_util
        G4hadronic_deex_management
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_util
        G4hepnumerics
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4track
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

