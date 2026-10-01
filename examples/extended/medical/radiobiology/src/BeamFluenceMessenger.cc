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
/// \file radiobiology/src/BeamFluenceMessenger.cc
/// \brief Implementation of the RadioBio::BeamFluenceMessenger class

#include "BeamFluenceMessenger.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

#include "BeamFluence.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamFluenceMessenger::BeamFluenceMessenger(BeamFluence* beamFluence)
  : G4UImessenger(), fBeamFluence(beamFluence)
{
  fDirectory = new G4UIdirectory("/beamFluence/");
  fDirectory->SetGuidance("Commands to configure voxel-resolved particle-entry counting.");

  fCalculationCmd = new G4UIcmdWithABool("/beamFluence/calculate", this);
  fCalculationCmd->SetGuidance("Enable or disable voxel-resolved particle-entry counting.");
  fCalculationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCalculationCmd->SetToBeBroadcasted(false);

  fVerbosityCmd = new G4UIcmdWithAnInteger("/beamFluence/verbose", this);
  fVerbosityCmd->SetGuidance("Set BeamFluence verbosity: 0 quiet, 1 normal, 2 debug.");
  fVerbosityCmd->SetToBeBroadcasted(false);

  fResetCmd = new G4UIcmdWithoutParameter("/beamFluence/reset", this);
  fResetCmd->SetGuidance("Reset accumulated voxel particle-entry counts.");
  fResetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fResetCmd->SetToBeBroadcasted(false);

  fPathCmd = new G4UIcmdWithAString("/beamFluence/fileName", this);
  fPathCmd->SetGuidance("Set the output filename for particle-entry counts.");
  fPathCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPathCmd->SetToBeBroadcasted(false);

  fPrintCmd = new G4UIcmdWithoutParameter("/beamFluence/print", this);
  fPrintCmd->SetGuidance("Print BeamFluence parameters.");
  fPrintCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPrintCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

BeamFluenceMessenger::~BeamFluenceMessenger()
{
  delete fDirectory;
  delete fCalculationCmd;
  delete fVerbosityCmd;
  delete fResetCmd;
  delete fPathCmd;
  delete fPrintCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void BeamFluenceMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fCalculationCmd) {
    fBeamFluence->SetCalculationEnabled(fCalculationCmd->GetNewBoolValue(newValue));
  }
  else if (command == fVerbosityCmd) {
    fBeamFluence->SetVerboseLevel(fVerbosityCmd->GetNewIntValue(newValue));
  }
  else if (command == fResetCmd) {
    fBeamFluence->Reset();
  }
  else if (command == fPathCmd) {
    fBeamFluence->SetPath(newValue);
  }
  else if (command == fPrintCmd) {
    fBeamFluence->PrintParameters();
  }
}

}  // namespace RadioBio
