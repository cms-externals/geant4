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

#include "MyDetectorMessenger.hh"

#include "MyDetectorConstruction.hh"
#include "G4UIparameter.hh"
#include "globals.hh"

MyDetectorMessenger::MyDetectorMessenger(MyDetectorConstruction * myDet)
:myDetector(myDet)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/mydet/",this);
  command->SetGuidance("My detector control.");
  fpMyDetCommandDirectory = command;

  command = new G4UIcommand("/mydet/calMaterial",this);
  command->SetGuidance("Define material of the calorimeter.");
  command->SetGuidance(" Available materials :");
  command->SetGuidance("   Air (defaust), Al, Fe, Pb");
  param = new G4UIparameter("material",'c',true);
  param->SetDefaultValue("Air");
  command->SetParameter(param);
  fpCalMaterialCommand = command;
}

MyDetectorMessenger::~MyDetectorMessenger () {
  delete fpCalMaterialCommand;
  delete fpMyDetCommandDirectory;
}

void MyDetectorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if (command == fpCalMaterialCommand)
  {
    myDetector->SetCalMaterial(newValues);
  }
}

