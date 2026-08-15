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
//

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4ios.hh"
#include "globals.hh"

#include "PhysicsList.hh"
#include "Tst50DetectorConstruction.hh"
#include "Tst50EventAction.hh"
#include "Tst50PrimaryGeneratorAction.hh"
#include "Tst50RunAction.hh"
#include "Tst50SteppingAction.hh"

int main(int argc, char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = argc == 1 ? new G4UIExecutive(argc, argv) : nullptr;

  // Construct a serial run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  Tst50DetectorConstruction* tst50Detector = new Tst50DetectorConstruction();
  runManager->SetUserInitialization(tst50Detector);

  PhysicsList* tst50Physics = new PhysicsList();
  runManager->SetUserInitialization(tst50Physics);

  Tst50PrimaryGeneratorAction* tst50PrimaryParticle = new Tst50PrimaryGeneratorAction();
  runManager->SetUserAction(tst50PrimaryParticle);

  Tst50RunAction* tst50Run = new Tst50RunAction();
  runManager->SetUserAction(tst50Run);

  Tst50EventAction* tst50EventAction = new Tst50EventAction();

  runManager->SetUserAction(tst50EventAction);

  Tst50SteppingAction* tst50SteppingAction =
    new Tst50SteppingAction(tst50PrimaryParticle, tst50Run, tst50Detector);
  runManager->SetUserAction(tst50SteppingAction);

  // get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if (ui == nullptr)
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UI->ApplyCommand(command + fileName);
  }
  else
  {
    // interactive mode
    G4VisExecutive* visManager = new G4VisExecutive;
    visManager->Initialize();
    ui->SessionStart();
    delete ui;
    delete visManager;
  }

  delete runManager;
  return 0;
}
