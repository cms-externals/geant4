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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "EventAction.hh"
#include "DetectorParameters.hh"
#include "ApplicationParameters.hh"
#include "Analysis.hh"

#include "G4GenericMessenger.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

using namespace ApplicationParameters;
using namespace DetectorParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(EventAction* eventAction)
 : G4UserRunAction(),
   fEventAction(eventAction)
{
  // set printing run numbers only
  G4RunManager::GetRunManager()->SetPrintProgress(0);

  // Test writing
  // Create analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // Create directories
  // analysisManager->SetHistoDirectoryName("histograms");
  // analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(VerboseLevel);
  if ( MergeNtuple) {
    analysisManager->SetNtupleMerging(true);
  }
  // The filename will be set in a test macro
  // analysisManager->SetFileName(Output);

  // Book histograms, ntuple

  // H1
  for (G4int i = 0; i < NofH1; ++i) {
    G4String name = "Edep-";
    name += std::to_string(i);
    analysisManager->CreateH1(name, "Edep in absorber",  50, 0., 200*MeV);
  }

  // H2
  for (G4int i = 0; i < NofH2; ++i) {
    G4String name = "PX-PY-";
    name += std::to_string(i);
    analysisManager->CreateH2(name, "PX PY 0 in absorber layer #2",
                              50, -10.*MeV, 10*MeV, 50, -10.*MeV, 10*MeV);
  }

  // H3
  for (G4int i = 0; i < NofH3; ++i) {
    G4String name = "PX-PY-PZ-";
    name += std::to_string(i);
    analysisManager->CreateH3(name, "PX PY PZ in absorber layer #2",
                              50, -10.*MeV, 10*MeV, 50, -10.*MeV, 10*MeV,
                              50, -10.*MeV, 10*MeV);
  }

  // P1
  for (G4int i = 0; i < NofP1; ++i) {
    G4String name = "Profile-Edep-vs-Z-";
    name += std::to_string(i);
    analysisManager->CreateP1(name, "Longit Edep pseudo-profile",
                              NofZLayers, 0, NofZLayers*LayerZSize, 0., 10.*GeV);
  }

  // P2
  for (G4int i = 0; i < NofP2; ++i) {
    G4String name = "Profile-Edep-vs-XY-";
    name += std::to_string(i);
    analysisManager->CreateP2(name, "Longit Edep pseudo-profile",
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              NofXYLayers, - NofXYLayers*LayerXYSize/2, NofXYLayers*LayerXYSize/2,
                              0, 10*GeV);
  }

  // Ntuples
  if ( NofNtuple > 0 ) {
    analysisManager->CreateNtuple("EDep", "Edep and Track Length");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->CreateNtupleSColumn("Label");
    analysisManager
     ->CreateNtupleDColumn("EdepVector", fEventAction->GetEdepVector());
    analysisManager
     ->CreateNtupleDColumn("TrackLengthVector", fEventAction->GetTrackLengthVector());
    analysisManager->FinishNtuple();
  }

  for (G4int i=0; i<6; ++i ) {
    if (NofNtuple > i) {
      G4String name = "TrackL-";
      name += std::to_string(i+1);
      analysisManager->CreateNtuple(name, "Track Length");
      analysisManager->CreateNtupleDColumn("Labs");
      analysisManager->FinishNtuple();
    }
  }

  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintStatistics()
{
  // print histogram statistics

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if(isMaster) {
    // H1
    for (G4int i = 0; i < NofH1; ++i) {
      auto h1 = analysisManager->GetH1(i);
      if ( h1 != nullptr ) {
        G4cout << "   H1: mean: " << h1->mean() << G4endl;
      }
    }

    // H2
    for (G4int i = 0; i < NofH2; ++i) {
      auto h2 =  analysisManager->GetH2(i);
      if ( h2 != nullptr ) {
        G4cout << "   H2: mean_x: " << h2->mean_x()
          << "   mean_y: " << h2->mean_y() << G4endl;
      }
    }

    // H3
    for (G4int i = 0; i < NofH3; ++i) {
      auto h3 =  analysisManager->GetH3(i);
      if ( h3 != nullptr ) {
        G4cout << "   H3: mean_x: " << h3->mean_x()
          << "   mean_y: " << h3->mean_y()
          << "   mean_z: " << h3->mean_z() << G4endl;
      }
    }

    for (G4int i = 0; i < NofP1; ++i) {
      auto p1 =  analysisManager->GetP1(i);
      if ( p1 != nullptr ) {
        G4cout << "   P1: mean: " << p1->mean() << G4endl;
      }
    }

    for (G4int i = 0; i < NofP2; ++i) {
      auto p2 =  analysisManager->GetP2(i);
      if ( p2 != nullptr  ) {
        G4cout << "   P2: mean_x: " << p2->mean_x()
          << "   mean_y: " << p2->mean_y() << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::DefineCommands()
{
  // Define /runAction command directory using generic messenger class
  fMessenger
    = new G4GenericMessenger(this,
                             "/runAction/",
                             "Run action commands");

  // printStatistic command
  auto& printStatisticCmd
    = fMessenger->DeclareMethod("printStatistics",
                                &RunAction::PrintStatistics,
                                "Print statistic at the end of Run.");
  printStatisticCmd.SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
