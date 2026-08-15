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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include "Analysis.hh"
#include "ApplicationParameters.hh"
#include "CalorHit.hh"
#include "CalorimeterSD.hh"

#include <iomanip>

using namespace ApplicationParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction() : G4UserEventAction(), fZAbsHCID(-1), fEdepVector(), fTrackLengthVector()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CalorHitsCollection* EventAction::GetHitsCollection(G4int hcID, const G4Event* event) const
{
  CalorHitsCollection* hitsCollection =
    static_cast<CalorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));

  if (!hitsCollection)
  {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("EventAction::GetHitsCollection()", "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::PrintEventStatistics(G4double absoEdep, G4double absoTrackLength) const
{
  // print event statistics
  G4cout << "   Absorber: total energy: " << std::setw(7) << G4BestUnit(absoEdep, "Energy")
         << "       total track length: " << std::setw(7) << G4BestUnit(absoTrackLength, "Length")
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* event)
{
  // Get hits collections IDs (only once)
  if (fZAbsHCID == -1)
  {
    fZAbsHCID = G4SDManager::GetSDMpointer()->GetCollectionID("zHitsCollection");
  }

  if (!fEdepVector.size())
  {
    // Get hits collections
    CalorHitsCollection* absoHC = GetHitsCollection(fZAbsHCID, event);

    // Add an element to both vectors for each hit
    for (std::size_t i = 0; i < absoHC->entries(); ++i)
    {
      fEdepVector.push_back(0);
      fTrackLengthVector.push_back(0);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
  // Get hits collections
  CalorHitsCollection* absoHC = GetHitsCollection(fZAbsHCID, event);

  // Get hit with total values
  CalorHit* absoHit = (*absoHC)[absoHC->entries() - 1];

  // Fill test vector with absoHit Edep
  for (std::size_t i = 0; i < absoHC->entries(); ++i)
  {
    CalorHit* calorHit = (*absoHC)[i];
    fEdepVector[i] = calorHit->GetEdep();
    fTrackLengthVector[i] = calorHit->GetTrackLength();
  }

  // Print per event (modulo n)
  //
  G4int eventID = event->GetEventID();
  G4int printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ((printModulo > 0) && (eventID % printModulo == 0))
  {
    G4cout << "---> End of event: " << eventID << G4endl;
    PrintEventStatistics(absoHit->GetEdep(), absoHit->GetTrackLength());
  }

  // Generate label for this event
  std::ostringstream label;
  label << "event" << eventID;

  // Fill histograms, ntuple
  //

  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  for (G4int i = 0; i < NofH1; ++i)
  {
    analysisManager->FillH1(i, absoHit->GetEdep());
  }

  // fill ntuple
  if (NofNtuple > 0)
  {
    analysisManager->FillNtupleDColumn(0, 0, absoHit->GetEdep());
    analysisManager->FillNtupleDColumn(0, 1, absoHit->GetTrackLength());
    // analysisManager->FillNtupleSColumn(0, 2, label.str());
    analysisManager->AddNtupleRow();
  }

  for (G4int i = 0; i < 4; ++i)
  {
    if (NofNtuple > i)
    {
      analysisManager->FillNtupleDColumn(i + 1, 0, absoHit->GetTrackLength());
      analysisManager->AddNtupleRow(i + 1);
    }
  }

  // print info
  if ((NofNtuple > 0) && (eventID < 10))
  {
    G4cout << "   EdepVector: ";
    for (G4int i = 0; i < G4int(fEdepVector.size()); ++i)
    {
      G4cout << fEdepVector[i] << "  ";
    }
    G4cout << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
