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
// --------------------------------------------------------------

#include "globals.hh"
#include "G4ios.hh"

#include "TstVADetectorMessenger.hh"

#include "TstVADetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

TstVADetectorMessenger::TstVADetectorMessenger(TstVADetectorConstruction * myDC)
:myDetector(myDC)
{
  G4String defParam;

  mydetDir = new G4UIdirectory("/geom/");
  mydetDir->SetGuidance("Geometry setup commands.");

  selDetCmd = new G4UIcmdWithAString("/geom/select",this);
  selDetCmd->SetGuidance("Select the way detector geometry is built.");
  selDetCmd->SetGuidance("  classic:  use simple placements ");
  selDetCmd->SetGuidance("  assembly: use assembly volume "  );
  selDetCmd->SetParameterName("choice",true);
  selDetCmd->SetDefaultValue("classic");
  selDetCmd->SetCandidates("classic assembly");
  selDetCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  switchCmd = new G4UIcmdWithAString("/geom/switch",this);
  switchCmd->SetGuidance("Assign the selected geometry to G4RunManager.");
  switchCmd->SetGuidance("In case the choice is present to this command,");
  switchCmd->SetGuidance("\"/geom/select\" will be invoked and then switched.");
  switchCmd->SetParameterName("choice",true);
  switchCmd->SetDefaultValue(" ");
  switchCmd->SetCandidates("classic assembly \" \"");
  switchCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  selMatCmd = new G4UIcmdWithAString("/geom/material",this);
  selMatCmd->SetGuidance("UNUSED IN THIS VERSION...\n\n");
  selMatCmd->SetGuidance("Select Material of the SimpleBox.");
  selMatCmd->SetGuidance("  Choice : Air, Al, Pb (default)");
  selMatCmd->SetParameterName("choice",true);
  selMatCmd->SetDefaultValue("Pb");
  selMatCmd->SetCandidates("Air Al Pb");
  selMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  myDetector->SelectDetector(defParam="classic");
  myDetector->SelectMaterial(defParam="Pb");
}

void TstVADetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command == selDetCmd )
  {
    myDetector->SelectDetector(newValues);
  }
  if( command == switchCmd )
  {
    if(newValues=="classic" || newValues=="assembly")
    { myDetector->SelectDetector(newValues); }
    myDetector->SwitchDetector();
  }
  if( command == selMatCmd )
  {
    myDetector->SelectMaterial(newValues);
    //G4cout << "UNUSED IN THIS VERSION...\a" << G4endl;
  }
  return;
}

