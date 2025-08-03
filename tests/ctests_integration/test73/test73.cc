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

#include "BeamTestDetectorConstruction.hh"
// User defined event action
#include "BeamTestEventAction.hh"
#include "BeamTestPhysicsList.hh"
#include "BeamTestPrimaryGeneratorAction.hh"
#include "BeamTestStackingAction.hh"
// User defined run action
#include "BeamTestRunAction.hh"
#include "BeamTestFileReader.hh"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "Randomize.hh"
#include "G4PhysListFactory.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main(int argc,char** argv) 
{
  //Construct a serial run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);

  // Mandatory initialization classes
  BeamTestDetectorConstruction* detector = new BeamTestDetectorConstruction();
  runManager->SetUserInitialization(detector);

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = nullptr;
  G4String physName = "";

  //Physics List name defined via 3nd argument
  if (argc==3) { physName = argv[2]; }

  // Physics List name defined via environment variable
  if("" == physName) {
    char* path = std::getenv("PHYSLIST");
    if (path) { physName = G4String(path); }
  }

  //reference PhysicsList via its name
  if ("" != physName && factory.IsReferencePhysList(physName)) {
    phys = factory.GetReferencePhysList(physName);
  }

  //local Physics List
  if(nullptr == phys) { phys = new BeamTestPhysicsList(); }

  // define physics
  runManager->SetUserInitialization(phys);

  // User action classes
  BeamTestRunAction* run_action = new BeamTestRunAction();  
  runManager->SetUserAction(run_action);
	
  BeamTestEventAction* event_action = new BeamTestEventAction(run_action);
  runManager->SetUserAction(event_action);
	
  runManager->SetUserAction(new BeamTestPrimaryGeneratorAction());
    
  //Stacking Action
  runManager->SetUserAction(new BeamTestStackingAction());

  // Initialize G4 kernel
  runManager->Initialize();

  //get the pointer to the User Interface manager
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (argc==1)   // Define UI terminal for interactive mode
    {
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      G4VisManager* visManager = new G4VisExecutive();
      visManager->Initialize();
      ui->SessionStart();
      delete ui;
      delete visManager;
    }
  else if (argc>1) // Batch mode with 1 or more files
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }

  delete runManager;
  return 0;
}

