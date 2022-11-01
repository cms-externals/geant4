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
// SBT (Solids Batch Test) main program
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "G4InteractiveSolid.hh"
#include "SBTMessenger.hh"
#include "SBTvoxelMessenger.hh"
#include "SBTVisManager.hh"

/*
MEDERNACH Emmanuel
Aug 2000

You could now run SBT with an argument script
and exit SBT.
*/

int main(int argc,char *argv[])
{
        // Initialize visualization manager
        SBTVisManager *visManager = new SBTVisManager;
        visManager->Initialize();

	// Build our "interactive" volume,
	// a volume that is also a messenger
	G4InteractiveSolid *interactiveSolid = new G4InteractiveSolid( "/solid/" );
	
	// Build our batch test messenger.
	// A new batch test will be created by it.
	// The test messenger gets the target solid
	// from a solid query class G4QuerySolid, as specified in the
	// second argument
	SBTMessenger *testMessenger = new SBTMessenger( "/test/", (G4SolidQuery *)interactiveSolid, visManager );

	//
	// Build our voxel test messenger
	SBTvoxelMessenger *voxelMessenger = new SBTvoxelMessenger( "/voxel/", (G4SolidQuery *)interactiveSolid, visManager );
		
	// Give control to interactive terminal
  G4UIsession *session;
  
#ifdef G4UI_USE_TCSH
  session = new G4UIterminal(new G4UItcsh);      
#else
	session = new G4UIterminal;
#endif
	   
	if (argc > 1)  // when run with an argument, run each scripts and exit
	  {
	    G4UImanager * UI = G4UImanager::GetUIpointer();

	    for (int i=1;i<argc;i++)
	      {
		UI->ApplyCommand("/control/execute "+G4String(argv[i]));
	      }
	  }
	else
	  {
	    session->SessionStart();
	  }
	
	// All finished...
	delete session;
	delete testMessenger;
	delete voxelMessenger;
	delete interactiveSolid;
	return 0;
}
