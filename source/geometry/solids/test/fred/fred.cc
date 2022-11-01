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
// fred: GEANT4 test program
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "FredDetectorConstruction.hh"
#include "FredPhysicsList.hh"
#include "FredPrimaryGeneratorAction.hh"
#include "FredEventAction.hh"
#include "FredMessenger.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

int main(int argc,char *argv[])
{
	// Detect interactive mode (if no arguments) and define UI session
	//
	G4UIExecutive* ui = 0;
	if ( argc == 1 ) {
	  ui = new G4UIExecutive(argc, argv);
	}

	// Construct run manager
	G4RunManager *runManager = new G4RunManager;
	
	// Build our master control messenger
	FredMessenger *messenger = new FredMessenger;
	
	// Build our detector
	runManager->SetUserInitialization( new FredDetectorConstruction(messenger) );
	runManager->SetUserInitialization( new FredPhysicsList );
	

	// Initialize visualization manager
	G4VisManager *visManager = new G4VisExecutive;
	visManager->Initialize();

	// Define our generator
	runManager->SetUserAction( new FredPrimaryGeneratorAction(messenger) );
	
	// Do something interesting at the end of each event
	runManager->SetUserAction( new FredEventAction(messenger) );

	//get the pointer to the User Interface manager 
	G4UImanager * UIman = G4UImanager::GetUIpointer();  

	if (!ui)  // Batch mode
	{
	  for (int i=1;i<argc;i++)
	  {
	    UIman->ApplyCommand("/control/execute "+G4String(argv[i]));
	  }
	}
	else      // Interactive mode
	{
          UIman->ApplyCommand("/control/execute vis.mac");    
          ui->SessionStart();
          delete ui;
	}
	
	// All finished...
	delete visManager;
	delete runManager;

	return 0;
}
