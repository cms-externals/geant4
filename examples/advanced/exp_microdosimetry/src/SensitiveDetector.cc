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
// Authors: Susanna Guatelli, Francesco Romano and Giuliana Milluzzo
// susanna@uow.edu.au, francesco.romano@ct.infn.it, giuliana.milluzzo@ct.infn.it
//
//  code based on the basic example B2
//
#include "SensitiveDetector.hh"

#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"

SensitiveDetector::SensitiveDetector(const G4String& name, const G4String& hitsCollectionName,
                                     G4bool isThisAMicrodosimeter,
                                     AnalysisManager* analysis_manager,
				     G4bool isArray)
  : G4VSensitiveDetector(name), fHitsCollection(NULL),isSiliconArray(isArray) 
    

{
  collectionName.insert(hitsCollectionName);
  analysis = analysis_manager;
  isMicrod = isThisAMicrodosimeter;  // only false for the E-stage of the telescope
}

SensitiveDetector::~SensitiveDetector() {}

void SensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
  // Clear array scores ARRAY 
  fArrayScores.clear();

  // Create hits collection
  fHitsCollection = new SensitiveDetectorHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

G4bool SensitiveDetector::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // ----- Score SiliconArray hits separately for each sensitive volume ----- // 
  if (isSiliconArray)
    {
      if (!aStep) return false;

  // energy deposit
      G4double edep = aStep->GetTotalEnergyDeposit();
      if (edep <= 0.) return false;

      const G4StepPoint* preStep = aStep->GetPreStepPoint();
      if (!preStep) return false;

    // Each SiliconArray sensor has its own copy number.
      const G4int copyNo =
        preStep->GetTouchableHandle()->GetCopyNumber();

      G4double len = aStep->GetStepLength();

    // Apply the same path-length selection as in the existing detector code.
      const G4ParticleDefinition* particle =
        aStep->GetTrack()->GetParticleDefinition();

      const G4String particleName = particle->GetParticleName();
      const G4int particleZ = particle->GetAtomicNumber();

      if (particleZ < 1 && particleName != "neutron")
	{
	  len = 0.;
	}

    // Accumulate energy and path length for this sensor.
      ArrayScore& score = fArrayScores[copyNo];
      score.edep += edep;
      score.path += len;

      return true;
    }

   // ----- // ----- // ----- // ----- // -----
  
 // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  if (edep == 0.) return false;

  // Uncomment the following lines to only score specific physical volumes
  /*
  G4String volumeName = aStep -> GetPreStepPoint() -> GetPhysicalVolume()-> GetName();
  if(volumeName != "SV_phys1" && volumeName != "sen_bridge" && volumeName !=
  "physSensitiveBridgeVolume" ) return false;
  */

  SensitiveDetectorHit* newHit = new SensitiveDetectorHit();

  newHit->SetEdep(edep);

  // step length
  G4double len = aStep->GetStepLength();

  // only consider ions, protons, and neutrons for path length
  const G4ParticleDefinition* thisParticle = aStep->GetTrack()->GetParticleDefinition();
  G4String particleName = thisParticle->GetParticleName();
  G4int particleZ = thisParticle->GetAtomicNumber();

  if ((particleZ < 1) && (particleName != "neutron")) len = 0.;

  newHit->SetPath(len);

  fHitsCollection->insert(newHit);

  return true;
}

void SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
#ifdef ANALYSIS_USE
  // Initialisation of total energy deposition per event to zero
  G4double totalEdepInOneEvent = 0;
  G4double totalPathLengthInOneEvent = 0;

  // Store SiliconArray energy and path length separately for each SV
  if (isSiliconArray)
  {
    const G4Event* event = G4RunManager::GetRunManager()->GetCurrentEvent();

    if (!event) return;

    const G4int eventID = event->GetEventID();

    for (const auto& entry : fArrayScores)
      {
	const G4int copyNo = entry.first;
	const ArrayScore& score = entry.second;

	if (score.edep > 0.)
	  {
	    analysis->StoreEnergyDeposition(score.edep / keV, score.path / um, copyNo, eventID);
	  }
      }
    return;
  }

  G4int NbHits = fHitsCollection->entries();
  // G4cout << "number of hits " <<NbHits << G4endl;

  G4double edep, len;

  for (G4int i = 0; i < NbHits; i++)
  {
    edep = (*fHitsCollection)[i]->GetEdep();
    totalEdepInOneEvent = totalEdepInOneEvent + edep;

    len = (*fHitsCollection)[i]->GetPath();
    totalPathLengthInOneEvent = totalPathLengthInOneEvent + len;
  }

  G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  if (totalEdepInOneEvent != 0)
  {
    if (isMicrod == true)
      //analysis->StoreEnergyDeposition(totalEdepInOneEvent / keV, totalPathLengthInOneEvent / um,
      //                              eventID);

      analysis->StoreEnergyDeposition(totalEdepInOneEvent / keV, totalPathLengthInOneEvent / um,0,
                                      eventID);
    else
      analysis->StoreSecondStageEnergyDeposition(totalEdepInOneEvent / keV, eventID);
  }
#endif
}
