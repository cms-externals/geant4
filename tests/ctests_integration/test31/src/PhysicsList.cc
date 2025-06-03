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
//---------------------------------------------------------------------------
//
// ClassName:   PhysicsList
//
// Author:      V.Ivanchenko 03.05.2004
//
// Modified:
// 16.11.06 Use components from physics_lists subdirectory (V.Ivanchenko)
// 24.10.12 Migrate to the new stopping and ion physics (A.Ribon)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4SystemOfUnits.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"

#include "StepLimiterBuilder.hh"
#include "G4DecayPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList() 
{
  defaultCutValue = 1.*mm;
  pMessenger = new PhysicsListMessenger(this);

  // Add Physics builders
  RegisterPhysics(new G4EmStandardPhysics(verbose));
  RegisterPhysics(new G4DecayPhysics(verbose));
  emName = "emstandard_opt0";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() 
{
  delete pMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(name == emName) { return; }
  if(verbose > 0) {
    G4cout << "### PhysicsList Add Physics <" << name 
           << "> " << G4endl;
  }
  if (name == "emstandard_opt0") {
    ReplacePhysics(new G4EmStandardPhysics());
    emName = name;

  } else if (name == "emstandard_opt1") {
    ReplacePhysics(new G4EmStandardPhysics_option1());
    emName = name;

  } else if (name == "emstandard_opt2") {
    ReplacePhysics(new G4EmStandardPhysics_option2());
    emName = name;

  } else if (name == "emstandard_opt3") {
    ReplacePhysics(new G4EmStandardPhysics_option3());
    emName = name;

  } else if (name == "emstandard_opt4") {
    ReplacePhysics(new G4EmStandardPhysics_option4());
    emName = name;

  } else if (name == "emlivermore") {
    ReplacePhysics(new G4EmLivermorePhysics());
    emName = name;

  } else if (name == "empenelope") {
    ReplacePhysics(new G4EmPenelopePhysics());
    emName = name;

  } else if (name == "emstandardSS") {
    ReplacePhysics(new G4EmStandardPhysicsSS());
    emName = name;

  } else if (name == "emstandardGS") {
    ReplacePhysics(new G4EmStandardPhysicsGS());
    emName = name;

  } else if (name == "emstandardWVI") {
    ReplacePhysics(new G4EmStandardPhysicsWVI());
    emName = name;

  } else if (name == "step_limit" && !stepLimiterIsRegisted) {
    RegisterPhysics(new StepLimiterBuilder());
    stepLimiterIsRegisted = true;

  } else if (name == "elastic" && !helIsRegisted) {
    RegisterPhysics(new G4HadronElasticPhysics());
    helIsRegisted = true;
    
  } else if (name == "binary" && !bicIsRegisted) {
    RegisterPhysics(new G4HadronInelasticQBBC());
    bicIsRegisted = true;
    
  } else if (name == "binary_ion" && !ionIsRegisted) {
    RegisterPhysics(new G4IonPhysics());
    ionIsRegisted = true;

  } else if (name == "gamma_nuc" && !gnucIsRegisted) {
    RegisterPhysics(new G4EmExtraPhysics());
    gnucIsRegisted = true;

  } else if (name == "stopping" && !stopIsRegisted) {
    RegisterPhysics(new G4StoppingPhysics());
    gnucIsRegisted = true;
    
  } else {
    G4cout << "PhysicsList::AddPhysicsList <" << name << ">" 
           << " fail - module is already regitered or is unknown " << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
