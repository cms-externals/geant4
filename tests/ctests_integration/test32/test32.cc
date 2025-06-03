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
/// \file testAnalysis.cc
/// \brief Main program of the  testAnalysis

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "ApplicationParameters.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT_HP.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

using namespace ApplicationParameters;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " testAnalysis "
           << "    [-m macro ]" << G4endl
           << "    [-nh1  value]" << G4endl
           << "    [-nh2  value]" << G4endl
           << "    [-nh3  value]" << G4endl
           << "    [-np1  value]" << G4endl
           << "    [-np2  value]" << G4endl
           << "    [-nnt  value]" << G4endl
           << "    [-merge on|off]" << G4endl
           << "    [-mWrite on|off]" << G4endl
           << "    [-vl   value]" << G4endl
           << G4endl;
  }

  G4int GetValue(const G4String& value) {
    return std::stoi(value);
  }

  G4bool GetOption(const G4String& option) {
    if      ( option == "on" )  return true;
    else if ( option == "off" ) return false;
    else  {
      PrintUsage();
      return 1;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  G4String macro;
  for ( G4int i=1; i<argc; i=i+2 ) {
    G4cout << "... testing " << G4String(argv[i]) << G4endl;
    if      ( G4String(argv[i]) == "-m" )   macro = argv[i+1];
    else if ( G4String(argv[i]) == "-nh1" ) NofH1 = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-nh2" ) NofH2 = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-nh3" ) NofH3 = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-np1" ) NofP1 = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-np2" ) NofP2 = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-nn1" ) NofNtuple = GetValue(argv[i+1]);
    else if ( G4String(argv[i]) == "-merge" ) MergeNtuple = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-mWrite" ) MultipleWrite = GetOption(argv[i+1]);
    else if ( G4String(argv[i]) == "-vl" )    VerboseLevel = GetValue(argv[i+1]);
    else {
      PrintUsage();
      return 1;
    }
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
  auto* runManager = G4RunManagerFactory::CreateRunManager();

  // Set mandatory initialization classes
  //
  DetectorConstruction* detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);

  ActionInitialization* actionInitialization
     = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {
    // interactive mode : define UI session
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !

  delete visManager;
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
