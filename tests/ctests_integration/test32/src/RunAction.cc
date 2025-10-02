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
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

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
    analysisManager->CreateH1(name, "Edep in absorber",  10, 0., 100*MeV);
  }

  // H2
  for (G4int i = 0; i < NofH2; ++i) {
    G4String name = "PX-PY-";
    name += std::to_string(i);
    analysisManager->CreateH2(name, "PX PY 0 in absorber layer #2",
                              10, -10.*MeV, 10*MeV, 10, -10.*MeV, 10*MeV);
  }

  // H3
  for (G4int i = 0; i < NofH3; ++i) {
    G4String name = "PX-PY-PZ-";
    name += std::to_string(i);
    analysisManager->CreateH3(name, "PX PY PZ in absorber layer #2",
                              10, -10.*MeV, 10*MeV, 10, -10.*MeV, 10*MeV,
                              10, -10.*MeV, 10*MeV);
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
    // analysisManager->CreateNtupleSColumn("Label");
    analysisManager
     ->CreateNtupleDColumn("EdepVector", fEventAction->GetEdepVector());
    analysisManager
     ->CreateNtupleDColumn("TrackLengthVector", fEventAction->GetTrackLengthVector());
    analysisManager->FinishNtuple();
  }

  for (G4int i=0; i<4; ++i ) {
    if (NofNtuple > i) {
      G4String name = "TrackL-";
      name += std::to_string(i+1);
      analysisManager->CreateNtuple(name, "Track Length");
      analysisManager->CreateNtupleDColumn("Labs");
      analysisManager->FinishNtuple();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TestWriting() const
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintStatistics() const
{
  // print histogram statistics

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  if(isMaster) {
    // Check only first histogram for all
    G4int id = 0;

    // H1
    auto h1 = analysisManager->GetH1(id);
    if ( h1 ) {
      G4cout << "   H1: "
             << "   mean: " << h1->mean() << " rms: " << h1->rms() << G4endl
             << "       "
             << "       "
             << "   id: " << id << " id by name:: "
             << analysisManager->GetH1Id(analysisManager->GetH1Name(id)) << G4endl;
    }

    // H2
    auto h2 =  analysisManager->GetH2(id);
    if ( h2 ) {
       G4cout << "   H2: "
              << "   mean_x: " << h2->mean_x() << " rms_x: " << h2->rms_x()
              << G4endl << "       "
              << "   mean_y: " << h2->mean_y() << " rms_y: " << h2->rms_y()
              << G4endl
              << "       "
              << "   id: " << id << " id by name:: "
              << analysisManager->GetH2Id(analysisManager->GetH2Name(id)) << G4endl;
    }

    auto h3 =  analysisManager->GetH3(id);
    if ( h3 ) {
      G4cout << "   H3: "
             << "   mean_x: " << h3->mean_x() << " rms_x: " << h3->rms_x()
             << G4endl << "       "
             << "   mean_y: " << h3->mean_y() << " rms_y: " << h3->rms_y()
             << G4endl << "       "
             << "   mean_z: " << h3->mean_z() << " rms_z: " << h3->rms_z()
             << G4endl
             << "       "
             << "   id: " << id << " id by name:: "
             << analysisManager->GetH3Id(analysisManager->GetH3Name(id)) << G4endl;
    }

    auto p1 =  analysisManager->GetP1(id);
    if ( p1 ) {
      G4cout << "   P1: "
             << "   mean: " << p1->mean() << " rms: " << p1->rms()
             << G4endl
             << "       "
             << "   id: " << id << " id by name:: "
             << analysisManager->GetP1Id(analysisManager->GetP1Name(id)) << G4endl;
    }

    auto p2 =  analysisManager->GetP2(id);
    if ( p2 ) {
      G4cout << "   P2: "
             << "   mean_x: " << p2->mean_x() << " rms_x: " << p2->rms_x()
             << G4endl << "       "
             << "   mean_y: " << p2->mean_y() << " rms_y: " << p2->rms_y()
             << G4endl
             << "       "
             << "   id: " << id << " id by name:: "
               << analysisManager->GetP2Id(analysisManager->GetP2Name(id)) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  if ((! MultipleWrite) ||
      (MultipleWrite && run->GetRunID() == 0 )) {
    // Default file name is set in macro
    // If multiple write is activated open file only at start of the first run
    G4AnalysisManager::Instance()->OpenFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
  // print some information
  PrintStatistics();

  G4AnalysisManager::Instance()->Write();

  // reset data (for the next write)
  if (MultipleWrite) {
    G4AnalysisManager::Instance()->Reset();
  }

  // extra write
  if (isMaster && ! MultipleWrite) {
    G4AnalysisManager::Instance()->WriteH1(8, "Edep-8-extra.csv");
    // G4AnalysisManager::Instance()->WriteH1(9, "Edep-9-extra.hdf5");
    G4AnalysisManager::Instance()->WriteH1(10, "Edep-10-extra.root");
    G4AnalysisManager::Instance()->WriteH1(11, "Edep-11-extra.xml");
  }

  if ((! MultipleWrite) ||
      (MultipleWrite && run->GetRunID() == 1 )) {
    G4AnalysisManager::Instance()->CloseFile();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
