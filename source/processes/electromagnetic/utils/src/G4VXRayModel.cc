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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4VXRayModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 28.04.2025
//
// -------------------------------------------------------------------
//

#include "G4VXRayModel.hh"
#include "G4LogicalVolume.hh"
#include "G4LossTableManager.hh"
#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VXRayModel::G4VXRayModel(const G4String& nam):
  pName(nam)
{
  pEmManager = G4LossTableManager::Instance();
  pEmManager->Register(this);
  isMaster = pEmManager->IsMaster();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VXRayModel::~G4VXRayModel()
{
  pEmManager->DeRegister(this);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::Initialise(std::vector<const G4LogicalVolume*>* ptr)
{
  pLogicalVolumes = ptr;
  auto param = G4EmParameters::Instance();
  verbose = isMaster ? param->Verbose() : param->WorkerVerbose();
  InitialiseModel();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::InitialiseModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VXRayModel::IsApplicable(const G4ParticleDefinition*,
				  const G4LogicalVolume*)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4VXRayModel::ModelDescription(std::ostream& outFile) const
{
  outFile << "The description for this model has not been written yet.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
