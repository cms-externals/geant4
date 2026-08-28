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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();

  G4int nThreads = 1;
  if (argc == 3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
  runManager->SetNumberOfThreads(nThreads);
  G4cout << "===== test55 is started with " << nThreads << " threads =====" << G4endl;

  // set mandatory initialization classes
  DetectorConstruction* det = new DetectorConstruction();
  PhysicsList* phys = new PhysicsList(det);
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(phys);

  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization(det));

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc == 1)  // Define UI terminal for interactive mode
  {
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    G4VisManager* visManager = new G4VisExecutive();
    visManager->Initialize();
    ui->SessionStart();
    delete ui;
    delete visManager;
  }
  else if (argc > 1)  // Batch mode with 1 or more files
  {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command + fileName);
  }

  // job termination
  //
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
