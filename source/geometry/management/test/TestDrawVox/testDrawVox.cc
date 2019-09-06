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
// --------------------------------------------------------------
//      GEANT4 - test01 
//
// --------------------------------------------------------------
// Comments
//
// 
// --------------------------------------------------------------

#include "TstDrawVox01DetectorConstruction.hh"
#include "TstDrawVox01RunAction.hh"
#include "TstDrawVox01PrimaryGeneratorAction.hh"
#include "TstDrawVox01SteppingAction.hh"
#include "TstDrawVox01PhysicsList.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

#include "G4ios.hh"

int main(int argc,char** argv) {

  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Set the default random engine to RanecuEngine
  CLHEP::RanecuEngine defaultEngine;
  CLHEP::HepRandom::setTheEngine(&defaultEngine);

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  runManager->SetUserInitialization(new TstDrawVox01DetectorConstruction);
  runManager->SetUserInitialization(new TstDrawVox01PhysicsList);

  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // UserAction classes
  runManager->SetUserAction(new TstDrawVox01RunAction);
  runManager->SetUserAction(new TstDrawVox01PrimaryGeneratorAction);
  runManager->SetUserAction(new TstDrawVox01SteppingAction);

  // User interactions
  // Define (G)UI for interactive mode
 
  // User interactions
  G4UImanager * UIman = G4UImanager::GetUIpointer();
  
  if(ui)
  { 
    // G4UIterminal is a (dumb) terminal.
    UIman->ApplyCommand("/control/execute vis.mac"); 
    ui->SessionStart();
    delete ui;
  }
  else  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UIman->ApplyCommand(command+fileName);
  }
  
  delete visManager;
  delete runManager;

  return 0;
}

