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
// * use.  Please see the license in the file  LICENSE and URL above  *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file radiobiology/src/DetectorMessenger.cc
/// \brief Implementation of the RadioBio::DetectorMessenger class

#include "DetectorMessenger.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"

#include "DetectorConstruction.hh"

namespace RadioBio
{

namespace
{
void ConfigureLengthCommand(G4UIcmdWithADoubleAndUnit* command, const G4String& guidance)
{
  command->SetGuidance(guidance);
  command->SetParameterName("Size", false);
  command->SetRange("Size>0.");
  command->SetUnitCategory("Length");
  command->AvailableForStates(G4State_PreInit, G4State_Idle);
  command->SetToBeBroadcasted(false);
}

void ConfigureVectorCommand(G4UIcmdWith3VectorAndUnit* command, const G4String& guidance)
{
  command->SetGuidance(guidance);
  command->SetParameterName("SizeAlongX", "SizeAlongY", "SizeAlongZ", false);
  command->SetRange("SizeAlongX>0. && SizeAlongY>0. && SizeAlongZ>0.");
  command->SetUnitCategory("Length");
  command->AvailableForStates(G4State_PreInit, G4State_Idle);
  command->SetToBeBroadcasted(false);
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector) : fDetector(detector)
{
  fGeometryDir = new G4UIdirectory("/detectorGeom/");
  fGeometryDir->SetGuidance("Commands controlling the air world and voxelized container.");
  fGeometryDir->SetGuidance("In Idle state, apply /voxels/update after changing dimensions.");

  fMaterCmd = new G4UIcmdWithAString("/detectorGeom/setMat", this);
  fMaterCmd->SetGuidance("Set the material of the voxelized container.");
  fMaterCmd->SetParameterName("choice", false);
  fMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterCmd->SetToBeBroadcasted(false);

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSize", this);
  ConfigureLengthCommand(fSizeCmd, "Set all three dimensions of the voxelized container.");

  fSizeVectorCmd = new G4UIcmdWith3VectorAndUnit("/detectorGeom/setBoxSizes", this);
  ConfigureVectorCommand(fSizeVectorCmd, "Set X, Y, and Z dimensions of the voxelized container.");

  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeX", this);
  ConfigureLengthCommand(fSizeXCmd, "Set the X dimension of the voxelized container.");

  fSizeYCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeY", this);
  ConfigureLengthCommand(fSizeYCmd, "Set the Y dimension of the voxelized container.");

  fSizeZCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeZ", this);
  ConfigureLengthCommand(fSizeZCmd, "Set the Z dimension of the voxelized container.");

  fWorldSizeCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setWorldSize", this);
  ConfigureLengthCommand(fWorldSizeCmd, "Set all three dimensions of the air world.");

  fWorldSizeVectorCmd = new G4UIcmdWith3VectorAndUnit("/detectorGeom/setWorldBoxSizes", this);
  ConfigureVectorCommand(fWorldSizeVectorCmd, "Set X, Y, and Z dimensions of the air world.");

  fWorldSizeXCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setWorldSizeX", this);
  ConfigureLengthCommand(fWorldSizeXCmd, "Set the X dimension of the air world.");

  fWorldSizeYCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setWorldSizeY", this);
  ConfigureLengthCommand(fWorldSizeYCmd, "Set the Y dimension of the air world.");

  fWorldSizeZCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setWorldSizeZ", this);
  ConfigureLengthCommand(fWorldSizeZCmd, "Set the Z dimension of the air world.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fWorldSizeZCmd;
  delete fWorldSizeYCmd;
  delete fWorldSizeXCmd;
  delete fWorldSizeVectorCmd;
  delete fWorldSizeCmd;

  delete fSizeZCmd;
  delete fSizeYCmd;
  delete fSizeXCmd;
  delete fSizeVectorCmd;
  delete fSizeCmd;
  delete fMaterCmd;

  delete fGeometryDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMaterCmd) {
    fDetector->SetMaterial(newValue);
  }
  else if (command == fSizeCmd) {
    fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSizeVectorCmd) {
    fDetector->SetSize(fSizeVectorCmd->GetNew3VectorValue(newValue));
  }
  else if (command == fSizeXCmd) {
    fDetector->SetSizeX(fSizeXCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSizeYCmd) {
    fDetector->SetSizeY(fSizeYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSizeZCmd) {
    fDetector->SetSizeZ(fSizeZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeCmd) {
    fDetector->SetWorldSize(fWorldSizeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeVectorCmd) {
    fDetector->SetWorldSize(fWorldSizeVectorCmd->GetNew3VectorValue(newValue));
  }
  else if (command == fWorldSizeXCmd) {
    fDetector->SetWorldSizeX(fWorldSizeXCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeYCmd) {
    fDetector->SetWorldSizeY(fWorldSizeYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fWorldSizeZCmd) {
    fDetector->SetWorldSizeZ(fWorldSizeZCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
