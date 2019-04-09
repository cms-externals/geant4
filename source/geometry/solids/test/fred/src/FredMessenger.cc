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
// FredMessenger.cc
//
// Implementation of fred's (primary) options
//

#include "FredMessenger.hh"
#include "FredTest3Messenger.hh"
#include "FredVoxelTestMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

//
// Constructor
//
FredMessenger::FredMessenger( )
{
	testVolume = 0;

	volumeNames[NATALIA] = "NATALIA";       // Yeah, yeah. Really cheezy, I know
	volumeNames[VOXEL]   = "VOXEL";

	volumeNames[BOX]     = "BOX";
	volumeNames[CONE]    = "CONE";
	volumeNames[CONE2]   = "CONE2";
	volumeNames[ORB]     = "ORB";
 	volumeNames[PARA]    = "PARA";
	volumeNames[SPHERE]  = "SPHERE";
 	volumeNames[TORUS1]  = "TORUS1";
 	volumeNames[TORUS2]  = "TORUS2";
	volumeNames[TRAP]    = "TRAP";
 	volumeNames[TRD]     = "TRD";
	volumeNames[TUBS]    = "TUBS";

 	volumeNames[ELLIPS]  = "ELLIPS";
 	volumeNames[ELCONE]  = "ELCONE";
 	volumeNames[ELTUBE]  = "ELTUBE";
 	volumeNames[EXTRU1]  = "EXTRU1";
 	volumeNames[EXTRU2]  = "EXTRU2";
 	volumeNames[EXTRU3]  = "EXTRU3";
 	volumeNames[EXTRU4]  = "EXTRU4";
 	volumeNames[HYPE]    = "HYPE";
 	volumeNames[PCON]    = "PCON";
	volumeNames[PCON2]   = "PCON2";
	volumeNames[PCON3]   = "PCON3";
	volumeNames[PCON4]   = "PCON4";
	volumeNames[PGON]    = "PGON";
	volumeNames[PGON2]   = "PGON2";
	volumeNames[PGON3]   = "PGON3";
 	volumeNames[TESSEL1] = "TESSEL1";
 	volumeNames[TESSEL2] = "TESSEL2";
 	volumeNames[TESSEL3] = "TESSEL3";
 	volumeNames[TESSEL4] = "TESSEL4";
 	volumeNames[TET]     = "TET";
 	volumeNames[TWBOX]   = "TWBOX";
 	volumeNames[TWTRAP1] = "TWTRAP1";
 	volumeNames[TWTRAP2] = "TWTRAP2";
 	volumeNames[TWTRD]   = "TWTRD";
 	volumeNames[TWTUBS]  = "TWTUBS";

	volumeNames[BOOL1]  =  "BOOL1";

	gunNames[SPRAY] = "SPRAY";
	gunNames[GRID]  = "GRID";
	gunNames[G4]    = "G4";

	drawNames[NORMAL] = "NORMAL";
	drawNames[SHADOW] = "SHADOW";

	//
	// Defaults
	//
	testVolumeType = BOX;
	gunType = SPRAY;
	drawType = NORMAL;
	
	startPhi = 0;
	deltaPhi = 360;
	
	numSide = 8;
	
	//
	// Declare directory
	//
	fredDirectory = new G4UIdirectory( "/fred/" );
	fredDirectory->SetGuidance( "Fred's options" );
	
	//
	// Volume command
	//
	volumeTypeNameCmd = new G4UIcmdWithAString( "/fred/volume", this );
	volumeTypeNameCmd->SetGuidance( "Test Volume Type" );
	volumeTypeNameCmd->SetParameterName( "VolumeType", true );
	
        G4String choices = volumeNames[NATALIA] + " "
			 + volumeNames[VOXEL] + " "
		         + volumeNames[BOX]  + " "
			 + volumeNames[CONE] + " "
			 + volumeNames[CONE2] + " "
			 + volumeNames[ORB] + " "
			 + volumeNames[PARA] + " "
			 + volumeNames[SPHERE] + " "
			 + volumeNames[TORUS1] + " "
			 + volumeNames[TORUS2] + " "
			 + volumeNames[TRAP] + " "
			 + volumeNames[TRD] + " "
                         + volumeNames[TUBS] + " "
			 + volumeNames[ELLIPS] + " " 
			 + volumeNames[ELCONE] + " " 
			 + volumeNames[ELTUBE] + " " 
			 + volumeNames[EXTRU1] + " " 
			 + volumeNames[EXTRU2] + " " 
			 + volumeNames[EXTRU3] + " " 
			 + volumeNames[EXTRU4] + " " 
			 + volumeNames[HYPE] + " " 
			 + volumeNames[PCON] + " "
			 + volumeNames[PCON2] + " "
			 + volumeNames[PCON3] + " "
			 + volumeNames[PCON4] + " "
			 + volumeNames[PGON] + " "
			 + volumeNames[PGON2] + " "
			 + volumeNames[PGON3] + " "
			 + volumeNames[TESSEL1] + " " 
			 + volumeNames[TESSEL2] + " " 
			 + volumeNames[TESSEL3] + " " 
			 + volumeNames[TESSEL4] + " " 
			 + volumeNames[TET] + " " 
			 + volumeNames[TWBOX] + " " 
			 + volumeNames[TWTRAP1] + " " 
			 + volumeNames[TWTRAP2] + " " 
			 + volumeNames[TWTRD] + " " 
			 + volumeNames[TWTUBS] + " " 
			 + volumeNames[BOOL1];
	volumeTypeNameCmd->SetCandidates( choices );
	
	volumeTypeNameCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
	
	//
	// gun command
	//
	gunTypeNameCmd = new G4UIcmdWithAString( "/fred/gun", this );
	gunTypeNameCmd->SetGuidance( "Type of particle gun to use" );
	gunTypeNameCmd->SetParameterName( "GunType", true );
	
        G4String choices2 = gunNames[SPRAY] + " "
			  + gunNames[GRID] + " "
			  + gunNames[G4];
	gunTypeNameCmd->SetCandidates( choices2 );
	
	//
	// draw command
	//
	drawTypeNameCmd = new G4UIcmdWithAString( "/fred/draw", this );
	drawTypeNameCmd->SetGuidance( "Type of drawing to make" );
	drawTypeNameCmd->SetParameterName( "DrawType", true );
	
        G4String choices3 = drawNames[NORMAL] + " "
			  + drawNames[SHADOW];
	drawTypeNameCmd->SetCandidates( choices3 );
	
	//
	// Start Phi Command
	//
	startPhiCmd = new G4UIcmdWithADouble( "/fred/startPhi", this );
	startPhiCmd->SetGuidance( "Starting phi value (degrees) for test volume" );
	startPhiCmd->SetParameterName( "StartPhi", true );
	startPhiCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
	
	//
	// Delta Phi Command
	//
	deltaPhiCmd = new G4UIcmdWithADouble( "/fred/deltaPhi", this );
	deltaPhiCmd->SetGuidance( "Delta phi value (degrees) for test volume" );
	deltaPhiCmd->SetParameterName( "DeltaPhi", true );
	deltaPhiCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
	
	//
	// Number sides command
	//
	numSideCmd = new G4UIcmdWithAnInteger( "/fred/numSide", this );
	numSideCmd->SetGuidance( "Number phi segments for test volume" );
	numSideCmd->SetParameterName( "NumberSide", true );
	numSideCmd->AvailableForStates( G4State_PreInit, G4State_Idle );

	//
	// Pause command
	//
	pauseCmd = new G4UIcmdWithoutParameter( "/fred/pause", this );
	pauseCmd->SetGuidance( "Prompts for return" );
	pauseCmd->AvailableForStates( G4State_PreInit, G4State_Idle );
	
	//
	// Declare test3 messenger
	//
	test3Messenger = new FredTest3Messenger();
	
	//
	// Declare voxel test messenger
	//
	voxelTestMessenger = new FredVoxelTestMessenger();
}

//
// Destructor
//
FredMessenger::~FredMessenger() 
{
	delete test3Messenger;
	delete voxelTestMessenger;
	delete volumeTypeNameCmd;
	delete fredDirectory;
}


//
// SetTestVolume
//
// Set the target volume for test 3. We need to funnel this down
// to our test 3 messenger
//
void FredMessenger::SetTestVolume( const G4VSolid *theTestVolume )
{
	testVolume = theTestVolume;
	test3Messenger->SetTestVolume( testVolume );
	voxelTestMessenger->SetTestVolume( testVolume );
}


//
// SelectedVolume
//
// Return selected test volume type
//
VolumeType FredMessenger::SelectedVolume()
{
	return(testVolumeType);
}

//
// SelectedGun
//
// Return selected gun
//
GunType FredMessenger::SelectedGun()
{
	return(gunType);
}

//
// SelectedDrawing
//
// Return selected drawing type
//
DrawType FredMessenger::SelectedDrawing()
{
	return( drawType );
}


//
// PauseInput
//
// This is to make up for a deficiency in the basic user-inteface:
// Wait for input from the user
//
void FredMessenger::PauseInput()
{
	G4cout << "Press <return> to continue: " << std::flush;
	char c;
	G4cin.get(c);
}


//
// SetNewValue
//
void FredMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == volumeTypeNameCmd) {
		for( G4int vol = 0; vol < FRED_VOLUMETYPE_NUM; vol++ ) {
		     	if (volumeNames[vol] == newValues) {
				testVolumeType = (VolumeType)vol;
				break;
			}
		}
	}
	else if (command == gunTypeNameCmd) {
		for( G4int gun = 0; gun < FRED_GUNTYPE_NUM; gun++ ) {
		     	if (gunNames[gun] == newValues) {
				gunType = (GunType)gun;
				break;
			}
		}
	}
	else if (command == drawTypeNameCmd) {
		for( G4int draw = 0; draw < FRED_DRAWTYPE_NUM; draw++ ) {
		     	if (drawNames[draw] == newValues) {
				drawType = (DrawType)draw;
				break;
			}
		}
	}
	else if (command == startPhiCmd) {
		startPhi = startPhiCmd->GetNewDoubleValue( newValues );
	}
	else if (command == deltaPhiCmd) {
		deltaPhi = deltaPhiCmd->GetNewDoubleValue( newValues );
	}
	else if (command == numSideCmd) {
		numSide = numSideCmd->GetNewIntValue( newValues );
	}
	else if (command == pauseCmd) {
		PauseInput();
	}
}

//
// GetCurrentValue
//
G4String FredMessenger::GetCurrentValue( G4UIcommand *command )
{
	if (command == volumeTypeNameCmd) {
		return volumeNames[testVolumeType];
	}
	else if (command == gunTypeNameCmd) {
		return gunNames[gunType];
	}
	else if (command == drawTypeNameCmd) {
		return drawNames[drawType];
	}
	else if (command == startPhiCmd) {
		return startPhiCmd->ConvertToString( startPhi );
	}
	else if (command == deltaPhiCmd) {
		return deltaPhiCmd->ConvertToString( deltaPhi );
	}
	else if (command == numSideCmd) {
		return numSideCmd->ConvertToString( numSide );
	}
	return "baloney";
}
