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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G24PhysicsList.hh"

#include "G4ChargedGeantino.hh"
#include "G4Electron.hh"
#include "G4Geantino.hh"
#include "G4PhysicsListHelper.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4Transportation.hh"
#include "G4UserSpecialCuts.hh"
#include "G4VUserPhysicsList.hh"

class G4Electron;
class G4Geantino;

G24PhysicsList::G24PhysicsList()
{
  G4StepLimiterPhysics* stepLimitPhys = new G4StepLimiterPhysics();
  RegisterPhysics(stepLimitPhys);
}

G24PhysicsList::~G24PhysicsList() {}

void G24PhysicsList::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
  G4Electron::ElectronDefinition();

  // G4Geantino::GeantinoDefinition();
  // G4ChargedGeantino::ChargedGeantinoDefinition();
}

void G24PhysicsList::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  G4VModularPhysicsList::ConstructProcess();
  AddTransportation();
  ph->RegisterProcess(new G4UserSpecialCuts(), G4Electron::ElectronDefinition());
}
