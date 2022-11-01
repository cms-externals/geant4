# - G4had_fission module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_fission
  PUBLIC_HEADERS
    G4FissLib.hh
    G4FissionLibrary.hh
    G4LFission.hh
    G4LLNLFission.hh
    G4fissionEvent.hh
  SOURCES
    G4FissLib.cc
    G4FissionLibrary.cc
    G4LFission.cc
    G4LLNLFission.cc
    G4SmpGEng.cc
    G4SmpIsoDir.cc
    G4SmpNEngCf252.cc
    G4SmpNVel.cc
    G4SmpNuDistDataPu239.cc
    G4SmpNuDistDataPu239_241.cc
    G4SmpNuDistDataPu239_241_MC.cc
    G4SmpNuDistDataU232_234_236_238.cc
    G4SmpNuDistDataU232_234_236_238_MC.cc
    G4SmpNuDistDataU233_235.cc
    G4SmpNuDistDataU233_235_MC.cc
    G4SmpNuDistDataU235.cc
    G4SmpNuDistDataU238.cc
    G4SmpNugDist.cc
    G4SmpPVel.cc
    G4SmpSpNuDistData.cc
    G4SmpSpNubarData.cc
    G4SmpSpNugDistData.cc
    G4SmpTerrell.cc
    G4SmpWatt.cc
    G4fissionEvent.cc
    G4fissionerr.cc
    G4rngc.cc)

geant4_module_link_libraries(G4had_fission
  PUBLIC
    G4baryons
    G4bosons
    G4globman
    G4had_par_hp
    G4hadronic_mgt
    G4hadronic_util
    G4heprandom
    G4leptons
    G4track)
