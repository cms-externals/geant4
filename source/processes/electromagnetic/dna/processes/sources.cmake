# - G4emdna-processes module build definition

# Define the Geant4 Module.
geant4_add_module(G4emdna-processes
  PUBLIC_HEADERS
    G4DNAAttachment.hh
    G4DNABrownianTransportation.hh
    G4DNAChargeDecrease.hh
    G4DNAChargeIncrease.hh
    G4DNADissociation.hh
    G4DNAElastic.hh
    G4DNAElectronSolvatation.hh
    G4DNAElectronSolvation.hh
    G4DNAElectronHoleRecombination.hh
    G4DNAExcitation.hh
    G4DNAPlasmonExcitation.hh
    G4DNAIonisation.hh
    G4DNAWaterDissociationDisplacer.hh
    G4DNAMolecularDissociation.hh
    G4DNAPositronium.hh
    G4DNARotExcitation.hh
    G4DNASecondOrderReaction.hh
    G4DNAVibExcitation.hh
    G4DNAScavengerProcess.hh
  SOURCES
    G4DNAAttachment.cc
    G4DNABrownianTransportation.cc
    G4DNAChargeDecrease.cc
    G4DNAChargeIncrease.cc
    G4DNADissociation.cc
    G4DNAElastic.cc
    G4DNAElectronSolvation.cc
    G4DNAElectronHoleRecombination.cc
    G4DNAExcitation.cc
    G4DNAPlasmonExcitation.cc
    G4DNAIonisation.cc
    G4DNAMolecularDissociation.cc
    G4DNAPositronium.cc
    G4DNARotExcitation.cc
    G4DNAWaterDissociationDisplacer.cc
    G4DNASecondOrderReaction.cc
    G4DNAVibExcitation.cc
    G4DNAScavengerProcess.cc)

geant4_module_link_libraries(G4emdna-processes
  PUBLIC
    G4baryons
    G4emdna-man
    G4emdna-models
    G4emdna-molman
    G4emdna-utils
    G4emutils
    G4leptons
  PRIVATE
    G4emdna-moltypes
    G4globman
    G4heprandom
    G4ions
    G4materials
    G4navigation
    G4partman
    G4track)
