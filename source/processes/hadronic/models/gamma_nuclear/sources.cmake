#------------------------------------------------------------------------------
# sources.cmake
# Module : G4had_gamm_nuclear
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4lepto_nuclear
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
GEANT4_DEFINE_MODULE(NAME G4had_gamma_nuclear
    HEADERS
	G4LENDorBERTModel.hh
    SOURCES
	G4LENDorBERTModel.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4geometrymng
        G4globman
	G4had_lend
        G4had_mod_man
        G4had_mod_util
        G4had_preequ_exciton
        G4had_string_diff
        G4had_string_frag
        G4had_string_man
        G4had_theo_max
        G4hadronic_bert_cascade
        G4hadronic_deex_evaporation
        G4hadronic_deex_fermi_breakup
        G4hadronic_deex_handler
        G4hadronic_deex_management
        G4hadronic_deex_multifragmentation
        G4hadronic_deex_photon_evaporation
        G4hadronic_deex_util
        G4hadronic_mgt
        G4hadronic_proc
        G4hadronic_util
        G4hadronic_xsect
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
        G4track
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4materials
        G4particles
        G4track
    LINK_LIBRARIES
)

# List any source specific properties here

