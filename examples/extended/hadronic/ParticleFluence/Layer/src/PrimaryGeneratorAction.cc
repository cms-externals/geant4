//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::
PrimaryGeneratorAction( const DetectorConstruction* pDetector ) :
  fPointerDetectorConstruction( pDetector )
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun( n_particle );
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  //***LOOKHERE*** Default particle and energy
  fParticleGun->SetParticleDefinition( particleTable->FindParticle( "geantino" ) );
  fParticleGun->SetParticleEnergy( 10.0*GeV );
  SetGunPosition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::SetGunPosition() const {
  // Shoot the particle in the middle between the world and the target layer
  G4double targetThickness =
    ( fPointerDetectorConstruction ? fPointerDetectorConstruction->GetThickness() : 0.0 );
  G4double gunPosition = -0.55*targetThickness;  //***LOOKHERE*** default gun position
                                                 //               along the z-axis
  G4cout << G4endl << "PrimaryGenerationAction::SetGunPosition() : gun position along z = " 
         << gunPosition << " mm " << G4endl << G4endl;  
  fParticleGun->SetParticlePosition( G4ThreeVector( 0.0, 0.0, gunPosition ) );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries( G4Event* anEvent ) {
  G4ThreeVector v( 0.0, 0.0, 1.0 );  //***LOOKHERE*** default shoot along the z-axis
  fParticleGun->SetParticleMomentumDirection( v );
  fParticleGun->GeneratePrimaryVertex( anEvent );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
