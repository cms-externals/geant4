# G4biasing_mgt module build definition

# Define the Geant4 Module.
geant4_add_module(G4biasing_mgt
  PUBLIC_HEADERS
    G4VProcessPlacer.hh
    G4ProcessPlacer.hh
    G4BiasingAppliedCase.hh
    G4BiasingOperationManager.hh
    G4VBiasingInteractionLaw.hh
    G4VBiasingOperation.hh
    G4VBiasingOperator.hh
    G4MaxTimeCuts.hh
    G4MinEkineCuts.hh
    G4SpecialCuts.hh
    G4UserSpecialCuts.hh
  SOURCES
    G4VProcessPlacer.cc
    G4ProcessPlacer.cc
    G4BiasingOperationManager.cc
    G4VBiasingOperation.cc
    G4VBiasingOperator.cc
    G4MaxTimeCuts.cc
    G4MinEkineCuts.cc
    G4SpecialCuts.cc
    G4UserSpecialCuts.cc)

geant4_module_link_libraries(G4biasing_mgt
  PUBLIC
    G4globman
    G4procman
    G4track
  PRIVATE
    G4emutils
    G4partman
    G4transportation)
