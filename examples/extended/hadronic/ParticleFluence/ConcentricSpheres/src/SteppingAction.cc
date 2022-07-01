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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4IonTable.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4VSolid.hh"
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"
#include "Run.hh"

const std::array< G4String, SteppingAction::numberScoringShells >
  SteppingAction::arrayScoringShellNames = { "tracker", "emCalo", "hadCalo" };

const std::array< G4String, SteppingAction::numberKinematicRegions >
  SteppingAction::arrayKinematicRegionNames = { "", "below 20 MeV", "above 20 MeV" };

const std::array< G4String, SteppingAction::numberScoringPositions >
  SteppingAction::arrayScoringPositionNames = { "forward", "backward" };

const std::array< G4String, SteppingAction::numberParticleTypes >
SteppingAction::arrayParticleTypeNames = { "all", "electron", "gamma", "muon", "neutrino",
                                           "pion", "neutron", "proton", "ion", "otherMeson",
                                           "otherBaryon" };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int SteppingAction::getIndex( const G4int iScoringShell, const G4int iKinematicRegion,
                                const G4int iScoringPosition, const G4int iParticleType ) {
  G4int index = -1;
  if ( iScoringShell >= 0     &&  iScoringShell < numberScoringShells        &&
       iKinematicRegion >= 0  &&  iKinematicRegion < numberKinematicRegions  &&
       iScoringPosition >= 0  &&  iScoringPosition < numberScoringPositions  &&
       iParticleType >= 0     &&  iParticleType < numberParticleTypes           ) {
    index = iScoringShell * numberKinematicRegions * numberScoringPositions * numberParticleTypes +
                                  iKinematicRegion * numberScoringPositions * numberParticleTypes +
                                                           iScoringPosition * numberParticleTypes +
                                                                              iParticleType;
  }
  return index;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction() :G4UserSteppingAction() {
  initialize();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::initialize() {
  // Initialization needed at the beginning of each Run
  fPrimaryParticleId = 0;
  fPrimaryParticleEnergy = 0.0;
  fPrimaryParticleDirection = G4ThreeVector( 0.0, 0.0, 1.0 );
  fTrackerMaterialName = fEmCaloMaterialName = fHadCaloMaterialName = "";
  fIsFirstStepOfTheEvent = true;
  fIsFirstStepInTracker = fIsFirstStepInEmCalo = fIsFirstStepInHadCalo = true;
  fIsFirstStepInScoringTrackerShell = fIsFirstStepInScoringEmCaloShell =
    fIsFirstStepInScoringHadCaloShell = true;  
  fCubicVolumeScoringTrackerShell = fCubicVolumeScoringEmCaloShell =
    fCubicVolumeScoringHadCaloShell = 1.0;
  for ( G4int i = 0; i < numberCombinations; ++i ) {
    fArraySumStepLengths[i] = 0.0;
  }
  /*
  for ( G4int i = 0; i < numberCombinations; ++i ) fArraySumStepLengths[i] = 999.9;
  G4cout << " numberCombinations=" << numberCombinations << G4endl;
  for ( G4int i = 0; i < numberScoringShells; ++i ) {
    for ( G4int j = 0; j < numberKinematicRegions; ++j ) {
      for ( G4int k = 0; k < numberScoringPositions; ++k ) {
        for ( G4int ll = 0; ll < numberParticleTypes; ++ll ) {
          G4int index = getIndex( i, j, k, ll );
          G4cout << "(i, j, k, ll)=(" << i << ", " << j << ", " << k << ", "
                 << ll << ")  ->" << index;
          if ( fArraySumStepLengths[ index ] < 1.0 ) G4cout << " <=== REPEATED!";
          else                                       fArraySumStepLengths[ index ] = 0.0;
          G4cout << G4endl;
        }
      }
    }
  }
  for ( G4int i = 0; i < numberCombinations; ++i ) {
    if ( fArraySumStepLengths[i] > 999.0 ) G4cout << " i=" << i << " NOT COVERED !" << G4endl;
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction( const G4Step* theStep ) {
  // Get information on the primary particle
  if ( fIsFirstStepOfTheEvent ) {
    if ( theStep->GetTrack()->GetParentID() == 0 ) {
      fPrimaryParticleId = theStep->GetTrack()->GetDefinition()->GetPDGEncoding();
      fPrimaryParticleEnergy = theStep->GetPreStepPoint()->GetKineticEnergy();
      fPrimaryParticleDirection = theStep->GetPreStepPoint()->GetMomentumDirection();
      if ( fRunPtr ) {
        fRunPtr->setPrimaryParticleId( fPrimaryParticleId );
        fRunPtr->setPrimaryParticleEnergy( fPrimaryParticleEnergy );
        fRunPtr->setPrimaryParticleDirection( fPrimaryParticleDirection );
      }
      fIsFirstStepOfTheEvent = false;
    }
  }
  // Get information on the materials
  if ( fIsFirstStepInTracker  &&
       theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physiTrackerShell" ) {
    fTrackerMaterialName = theStep->GetPreStepPoint()->GetMaterial()->GetName();
    if ( fRunPtr ) fRunPtr->setTrackerMaterialName( fTrackerMaterialName );
    fIsFirstStepInTracker = false;
  }
  if ( fIsFirstStepInEmCalo  &&
       theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physiEmCaloShell" ) {
    fEmCaloMaterialName = theStep->GetPreStepPoint()->GetMaterial()->GetName();
    if ( fRunPtr ) fRunPtr->setEmCaloMaterialName( fEmCaloMaterialName );
    fIsFirstStepInEmCalo = false;
  }
  if ( fIsFirstStepInHadCalo  &&
       theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() == "physiHadCaloShell" ) {
    fHadCaloMaterialName = theStep->GetPreStepPoint()->GetMaterial()->GetName();
    if ( fRunPtr ) fRunPtr->setHadCaloMaterialName( fHadCaloMaterialName );
    fIsFirstStepInHadCalo = false;
  }
  // Get information on step lengths in the scoring shells
  G4int iScoringShell = -1;
  if ( theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
       "physiScoringTrackerShell" ) {   
    iScoringShell = 0;
    if (  fIsFirstStepInScoringTrackerShell ) {
      fCubicVolumeScoringTrackerShell =
        theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      if ( fRunPtr ) fRunPtr->setCubicVolumeScoringTrackerShell( fCubicVolumeScoringTrackerShell );
      fIsFirstStepInScoringTrackerShell = false;
    }
  } else if ( theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
              "physiScoringEmCaloShell" ) {
    iScoringShell = 1;
    if (  fIsFirstStepInScoringEmCaloShell ) {
      fCubicVolumeScoringEmCaloShell =
        theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      if ( fRunPtr ) fRunPtr->setCubicVolumeScoringEmCaloShell( fCubicVolumeScoringEmCaloShell );
      fIsFirstStepInScoringEmCaloShell = false;
    }
  } else if ( theStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() ==
              "physiScoringHadCaloShell" ) {
    iScoringShell = 2;
    if (  fIsFirstStepInScoringHadCaloShell ) {
      fCubicVolumeScoringHadCaloShell =
        theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetSolid()->GetCubicVolume();
      if ( fRunPtr ) fRunPtr->setCubicVolumeScoringHadCaloShell( fCubicVolumeScoringHadCaloShell );
      fIsFirstStepInScoringHadCaloShell = false;
    }
  }
  if ( iScoringShell >= 0 ) {
    G4double stepLength = theStep->GetTrack()->GetStepLength() * theStep->GetTrack()->GetWeight();
    G4int absPdg = theStep->GetTrack()->GetDefinition() == nullptr ? 0 :
      std::abs( theStep->GetTrack()->GetDefinition()->GetPDGEncoding() );
    /*
    G4cout << theStep->GetTrack()->GetDefinition()->GetParticleName() << "  absPdg=" << absPdg
           << "  Ekin[MeV]=" << theStep->GetPreStepPoint()->GetKineticEnergy()
           << "  r[mm]=" << theStep->GetTrack()->GetPosition().mag()
           << "  z[mm]=" << theStep->GetTrack()->GetPosition().z()
           << "  " << theStep->GetTrack()->GetVolume()->GetName()
           << "  " << theStep->GetTrack()->GetMaterial()->GetName()
           << "  L[mm]=" << stepLength << "  " 
           << ( fPrimaryParticleDirection.dot( theStep->GetTrack()->GetPosition().unit() ) > 0.0
                ? "forward" : "backward" ) << G4endl;
    */
    // Three kinematical regions:  [0] : any value ;  [1] : below 20 MeV ;  [2] : above 20 MeV
    G4int iKinematicRegion = theStep->GetPreStepPoint()->GetKineticEnergy() < 20.0 ? 1 : 2;
    // Two scoring positions:  [0] : forward hemisphere ;  [1] : backward hemisphere
    // (with respect to the primary particle initial direction)
    G4int iScoringPosition =
      fPrimaryParticleDirection.dot( theStep->GetTrack()->GetPosition().unit() ) > 0.0 ? 0 : 1;
    G4int iParticleType = -1;
    if      ( absPdg == 11 ) iParticleType = 1;  // electron (and positron)
    else if ( absPdg == 22 ) iParticleType = 2;  // gamma
    else if ( absPdg == 13 ) iParticleType = 3;  // muons (mu- and mu+)
    else if ( absPdg == 12 || absPdg == 14 || absPdg == 16 ) iParticleType = 4;  // neutrinos
                                                           // (and anti-neutrinos), all flavors
    else if ( absPdg == 111 || absPdg == 211 ) iParticleType = 5;  // (charged) pions
    else if ( absPdg == 2112 ) iParticleType = 6;  // neutron (and anti-neutron)
    else if ( absPdg == 2212 ) iParticleType = 7;  // proton  (and anti-proton)
    else if ( G4IonTable::IsIon( theStep->GetTrack()->GetDefinition() ) || // ions (and anti-ions)
              G4IonTable::IsAntiIon( theStep->GetTrack()->GetDefinition() ) ) iParticleType = 8;
    else if ( absPdg < 1000 ) iParticleType = 9;   // other mesons (e.g. kaons) (Note: this works
                                                   // in most cases, but not always!)
    else if ( absPdg > 1000 ) iParticleType = 10;  // other baryons (e.g. hyperons, anti-hyperons,
                                                   // etc.)
    // Consider the specific case : scoring shell, kinematic region, scoring position, and
    // particle type
    G4int index = getIndex( iScoringShell, iKinematicRegion, iScoringPosition, iParticleType );
    fArraySumStepLengths[index] += stepLength;
    // Consider the "all" particle case, with the same scoring shell, kinematic region and
    // scoring position
    index = getIndex( iScoringShell, iKinematicRegion, iScoringPosition, 0 );
    fArraySumStepLengths[index] += stepLength;
    // Consider the "any" kinematic region case, with the same scoring shell, scoring position
    // and particle type    
    index = getIndex( iScoringShell, 0, iScoringPosition, iParticleType );
    fArraySumStepLengths[index] += stepLength;
    // Consider the "any" kinematic region and "all" particle, with the same scoring shell and
    // scoring position    
    index = getIndex( iScoringShell, 0, iScoringPosition, 0 );
    fArraySumStepLengths[index] += stepLength;
    if ( fRunPtr ) fRunPtr->setArray( fArraySumStepLengths );
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
