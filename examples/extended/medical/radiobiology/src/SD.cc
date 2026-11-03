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
//
/// \file radiobiology/src/SD.cc
/// \brief Implementation of the RadioBio::SD class

#include "SD.hh"

#include "G4AccumulableManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"

#include "Hit.hh"
#include "Manager.hh"
#include "VRadiobiologicalAccumulable.hh"
#include "VoxelizedSensitiveDetector.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SD::SD(const G4String& name, const G4String& hitsCollectionName) : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new RadioBioHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection to the hit collection of the current event
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  hce->AddHitsCollection(hcID, fHitsCollection);

  // Track IDs restart for every event, so the anti-double-counting cache
  // must also restart for every event.
  fCountedVoxelEntries.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Create a new hit
  RadioBio::Hit* newHit = new RadioBio::Hit();

  // Get the pre-step kinetic energy
  G4double eKinPre = aStep->GetPreStepPoint()->GetKineticEnergy();

  // Get the post-step kinetic energy
  G4double eKinPost = aStep->GetPostStepPoint()->GetKineticEnergy();

  // Get the average kinetic energy during the step
  G4double eKinMean = (eKinPre + eKinPost) * 0.5;

  // Get all secondary particles generated during the current step
  const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();

  const std::size_t secondarySize = secondary->size();

  G4double energySecondaryElectrons = 0.;

  if (secondarySize > 0) {
    for (std::size_t numsec = 0; numsec < secondarySize; ++numsec) {
      const G4Track* secondaryTrack = (*secondary)[numsec];

      const G4ParticleDefinition* secondaryDefinition = secondaryTrack->GetParticleDefinition();

      // Get the PDG code of the secondary particle
      const G4int pdgSecondary = secondaryDefinition->GetPDGEncoding();

      // Accumulate the kinetic energy of generated secondary electrons
      if (pdgSecondary == 11) {
        energySecondaryElectrons += secondaryTrack->GetKineticEnergy();
      }
    }
  }

  // Update the hit with the necessary quantities
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());

  newHit->SetParentID(aStep->GetTrack()->GetParentID());

  newHit->SetPartType(aStep->GetTrack()->GetParticleDefinition());

  newHit->SetEkinMean(eKinMean);

  newHit->SetDeltaE(aStep->GetTotalEnergyDeposit());

  newHit->SetEinit(aStep->GetPreStepPoint()->GetKineticEnergy());

  newHit->SetTrackLength(aStep->GetStepLength());

  newHit->SetElectronEnergy(energySecondaryElectrons);

  newHit->SetPhysicalVolume(aStep->GetPreStepPoint()->GetPhysicalVolume());

  newHit->SetVoxelIndexes(aStep->GetPreStepPoint()->GetTouchableHandle());

  // Count a particle when it enters a voxel through a geometrical boundary
  // or when it is created inside that voxel. Each track is counted at most
  // once in a given voxel during the event.
  G4bool voxelEntry = false;
  const G4Track* track = aStep->GetTrack();
  const G4StepPoint* preStepPoint = aStep->GetPreStepPoint();

  const auto* beamFluence = Manager::GetInstance()->GetQuantity("BeamFluence");
  if (beamFluence != nullptr && beamFluence->IsCalculationEnabled()
      && (preStepPoint->GetStepStatus() == fGeomBoundary || track->GetCurrentStepNumber() == 1))
  {
    const auto entryKey = std::make_tuple(track->GetTrackID(), newHit->GetXindex(),
                                          newHit->GetYindex(), newHit->GetZindex());

    voxelEntry = fCountedVoxelEntries.insert(entryKey).second;
  }

  newHit->SetVoxelEntry(voxelEntry);

  // Insert the hit into the hits collection
  fHitsCollection->insert(newHit);

  // Accumulables are accumulated only if calculation is enabled
  for (G4int i = 0; i < G4AccumulableManager::Instance()->GetNofAccumulables(); ++i) {
    // Get the accumulable from the accumulable manager
    G4VAccumulable* generalAccumulable = G4AccumulableManager::Instance()->GetAccumulable(i);

    VRadiobiologicalAccumulable* radiobiologicalAccumulable =
      dynamic_cast<VRadiobiologicalAccumulable*>(generalAccumulable);

    // Ignore accumulables that are not radiobiological accumulables
    if (radiobiologicalAccumulable == nullptr) {
      continue;
    }

    // Get the radiobiological quantity only after confirming the accumulable type.
    auto* quantity = Manager::GetInstance()->GetQuantity(generalAccumulable->GetName());
    if (quantity == nullptr) {
      continue;
    }

    // Accumulate only if the corresponding calculation is enabled
    if (quantity->IsCalculationEnabled()) {
      radiobiologicalAccumulable->Accumulate(newHit);
    }
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SD::EndOfEvent(G4HCofThisEvent*)
{
  if (verboseLevel > 1) {
    const G4int numberOfHits = fHitsCollection->entries();

    G4cout << G4endl << "--------> Hits Collection: in this event there are " << numberOfHits
           << " hits in the detector slices:" << G4endl;

    for (G4int i = 0; i < numberOfHits; ++i) {
      (*fHitsCollection)[i]->Print();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
