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
// GridParticleGunMessenger.cc
//
// Implemenation of GridParticleGun's messenger
//

#include "GridParticleGunMessenger.hh"
#include "GridParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4ios.hh"

//
// Constructor
// 
GridParticleGunMessenger::GridParticleGunMessenger( GridParticleGun *ourGun )
{
	gun = ourGun;
	
	gunDirectory = new G4UIdirectory( "/gridGun/" );
	gunDirectory->SetGuidance( "Grid Particle Gun control commands." );
	
	directionCmd = new G4UIcmdWith3Vector( "/gridGun/direction", this );
	directionCmd->SetGuidance( "Set direction of particles." );
	directionCmd->SetParameterName( "X", "Y", "Z", true, true );
	
	originCmd = new G4UIcmdWith3VectorAndUnit( "/gridGun/origin", this );
	originCmd->SetGuidance( "Set origin of grid." );
	originCmd->SetParameterName( "X", "Y", "Z", true, true );
	originCmd->SetDefaultUnit( "m" );	

	grid1Cmd = new G4UIcmdWith3VectorAndUnit( "/gridGun/grid1", this );
	grid1Cmd->SetGuidance( "Set first axis of grid." );
	grid1Cmd->SetParameterName( "X", "Y", "Z", true, true );
	grid1Cmd->SetDefaultUnit( "m" );

	grid2Cmd = new G4UIcmdWith3VectorAndUnit( "/gridGun/grid2", this );
	grid2Cmd->SetGuidance( "Set second axis of grid." );
	grid2Cmd->SetParameterName( "X", "Y", "Z", true, true );
	grid2Cmd->SetDefaultUnit( "m" );

	n1Cmd = new G4UIcmdWithAnInteger( "/gridGun/n1", this );
	n1Cmd->SetGuidance( "Set number of grid points along first axis." );
	n1Cmd->SetParameterName( "n1", true, true );
	
	n2Cmd = new G4UIcmdWithAnInteger( "/gridGun/n2", this );
	n2Cmd->SetGuidance( "Set number of grid points along second axis." );
	n2Cmd->SetParameterName( "n2", true, true );
}

//
// Destructor
//
GridParticleGunMessenger::~GridParticleGunMessenger()
{
	delete directionCmd;
	delete originCmd;
	delete grid1Cmd;
	delete grid2Cmd;
	delete n1Cmd;
	delete n2Cmd;
	delete gunDirectory;
}

//
// SetNewValue
//
// Called by the UI when user requests a change
//
void GridParticleGunMessenger::SetNewValue( G4UIcommand *command, G4String newValues )
{
	if (command == directionCmd) {
		gun->SetDirection( directionCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == originCmd) {
		gun->SetOrigin( originCmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == grid1Cmd) {
		gun->SetGrid1( grid1Cmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == grid2Cmd) {
		gun->SetGrid2( grid2Cmd->GetNew3VectorValue( newValues ) );
	}
	else if (command == n1Cmd) {
		gun->SetN1( n1Cmd->GetNewIntValue( newValues ) );
	}
	else if (command == n2Cmd) {
		gun->SetN2( n2Cmd->GetNewIntValue( newValues ) );
	}
}

//
// GetCurrentValue
//
G4String GridParticleGunMessenger::GetCurrentValue( G4UIcommand *command )
{
	if (command == directionCmd) {
		return directionCmd->ConvertToString( gun->GetDirection() );
	}
	else if (command == originCmd) {
		return originCmd->ConvertToString( gun->GetOrigin(), "m" );
	}
	else if (command == grid1Cmd) {
		return grid1Cmd->ConvertToString( gun->GetGrid1(), "m" );
	}
	else if (command == grid2Cmd) {
		return grid2Cmd->ConvertToString( gun->GetGrid2(), "m" );
	}
	else if (command == n1Cmd) {
		return n1Cmd->ConvertToString( gun->GetN1() );
	}
	else if (command == n2Cmd) {
		return n2Cmd->ConvertToString( gun->GetN2() );
	}
	
	return "baloney";
}

	
	
	
	
