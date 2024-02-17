# - G4partman module build definition

# Define the Geant4 Module.
geant4_add_module(G4partman
  PUBLIC_HEADERS
    G4DalitzDecayChannel.hh
    G4DecayProducts.hh
    G4DecayTable.hh
    G4DecayTableMessenger.hh
    G4DynamicParticle.hh
    G4DynamicParticle.icc
    G4ElectronOccupancy.hh
    G4HyperNucleiProperties.hh
    G4IonTable.hh
    G4Ions.hh
    G4IsotopeProperty.hh
    G4KL3DecayChannel.hh
    G4MuonicAtom.hh
    G4MuonicAtomHelper.hh
    G4MuonDecayChannel.hh
    G4MuonDecayChannelWithSpin.hh
    G4MuonRadiativeDecayChannelWithSpin.hh
    G4NeutronBetaDecayChannel.hh
    G4NucleiProperties.hh
    G4NucleiPropertiesTableAME12.hh
    G4NucleiPropertiesTheoreticalTable.hh
    G4NuclideTable.hh
    G4NuclideTableMessenger.hh
    G4PDGCodeChecker.hh
    G4PDefManager.hh
    G4ParticleDefinition.hh
    G4ParticleDefinition.icc
    G4ParticleMessenger.hh
    G4ParticleMomentum.hh
    G4ParticlePropertyData.hh
    G4ParticlePropertyData.icc
    G4ParticlePropertyMessenger.hh
    G4ParticlePropertyTable.hh
    G4ParticleTable.hh
    G4ParticleTable.icc
    G4ParticleTableIterator.hh
    G4ParticleWithCuts.hh
    G4ParticlesWorkspace.hh
    G4PhaseSpaceDecayChannel.hh
    G4PionRadiativeDecayChannel.hh
    G4PrimaryParticle.hh
    G4PrimaryVertex.hh
    G4TauLeptonicDecayChannel.hh
    G4VDecayChannel.hh
    G4VIsotopeTable.hh
    G4VUserPrimaryParticleInformation.hh
    G4VUserPrimaryVertexInformation.hh
    pwdefs.hh
  SOURCES
    G4DalitzDecayChannel.cc
    G4DecayProducts.cc
    G4DecayTable.cc
    G4DecayTableMessenger.cc
    G4DynamicParticle.cc
    G4ElectronOccupancy.cc
    G4HyperNucleiProperties.cc
    G4IonTable.cc
    G4Ions.cc
    G4IsotopeProperty.cc
    G4KL3DecayChannel.cc
    G4MuonicAtom.cc
    G4MuonicAtomHelper.cc
    G4MuonDecayChannel.cc
    G4MuonDecayChannelWithSpin.cc
    G4MuonRadiativeDecayChannelWithSpin.cc
    G4NeutronBetaDecayChannel.cc
    G4NucleiProperties.cc
    G4NucleiPropertiesTableAME12.cc
    G4NucleiPropertiesTheoreticalTableA.cc
    G4NucleiPropertiesTheoreticalTableB.cc
    G4NuclideTable.cc
    G4NuclideTableMessenger.cc
    G4PDGCodeChecker.cc
    G4PDefManager.cc
    G4ParticleDefinition.cc
    G4ParticleMessenger.cc
    G4ParticlePropertyData.cc
    G4ParticlePropertyMessenger.cc
    G4ParticlePropertyTable.cc
    G4ParticleTable.cc
    G4ParticlesWorkspace.cc
    G4PhaseSpaceDecayChannel.cc
    G4PionRadiativeDecayChannel.cc
    G4PrimaryParticle.cc
    G4PrimaryVertex.cc
    G4TauLeptonicDecayChannel.cc
    G4VDecayChannel.cc
    G4VIsotopeTable.cc)

geant4_module_link_libraries(G4partman PUBLIC G4globman G4hepgeometry G4heprandom G4intercoms)
