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
//
// -------------------------------------------------------------
//      GEANT4 test45
//
//  Application demonstrating Geant4 hadronic physics:
//  beam interaction with a thick target
//
//  Authors: A.Ivanchenko & V.Ivanchenko
//
//  Modified: Sept 2008
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv) 
{
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = (argc == 1) ? new G4UIExecutive(argc, argv) : nullptr;

  //Construct a serial run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  //set mandatory initialization classes
  runManager->SetUserInitialization(new DetectorConstruction());

  G4String physName = "";

  //Physics List name defined via 3nd argument
  if (argc==3) { physName = argv[2]; }

  // Physics List name defined via environment variable
  if ("" == physName) {
    char* path = std::getenv("PHYSLIST");
    if (path != nullptr) { physName = G4String(path); }
  }
  if ("" == physName) { physName = "FTFP_BERT"; }

  G4PhysListFactory factory;
  auto* phys = factory.GetReferencePhysList(physName);
  if (phys == nullptr) { exit(1); }
  runManager->SetUserInitialization(phys);

  runManager->SetUserAction(new PrimaryGeneratorAction());

  //set user action classes
  runManager->SetUserAction(new RunAction());
  runManager->SetUserAction(new EventAction());
  runManager->SetUserAction(new StackingAction());

  //get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (ui != nullptr) {
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    ui->SessionStart();
    delete ui;
    delete visManager;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
