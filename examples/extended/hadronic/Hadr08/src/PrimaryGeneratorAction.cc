/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), 
                                                   fParticleGun( nullptr ) {
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun( n_particle );
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle = particleTable->FindParticle( particleName="pi+" );
  fParticleGun->SetParticleDefinition( particle );
  fParticleGun->SetParticleMomentumDirection( G4ThreeVector( 0.0, 0.0, 1.0 ) );
  fParticleGun->SetParticleEnergy( 5.0*GeV );
  fParticleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, -50.0*cm ) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {
  fParticleGun->GeneratePrimaryVertex( anEvent );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ParticleGun* PrimaryGeneratorAction::GetParticleGun() {
  return fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
