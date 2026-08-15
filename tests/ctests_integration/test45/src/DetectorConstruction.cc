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
/////////////////////////////////////////////////////////////////////////
//
// DetectorConstruction
//
// Created: 31.01.2003 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

#include "CheckVolumeSD.hh"
#include "DetectorMessenger.hh"
#include "HistoManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  logicTarget = nullptr;
  logicCheck = nullptr;
  logicWorld = nullptr;
  world = nullptr;
  detectorMessenger = new DetectorMessenger(this);

  width = 2. * cm;

  targetMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
  worldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

  // Prepare sensitive detectors
  checkSD = new CheckVolumeSD("checkSD");
  (G4SDManager::GetSDMpointer())->AddNewDetector(checkSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete detectorMessenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if (world != nullptr)
  {
    return world;
  }
  // Cleanup old geometry

  // Sizes
  G4double checkXY = (width + 1.0 * mm) * 0.5;
  G4double worldXY = (width + 1.0 * cm) * 0.5;
  G4double targetZ = HistoManager::GetPointer()->Length() * 0.5;
  G4double checkZ = targetZ + 1.0 * mm;
  G4double worldZ = targetZ + 1.0 * cm;

  //
  // World
  //
  G4Box* solidW = new G4Box("World", worldXY, worldXY, worldZ);
  logicWorld = new G4LogicalVolume(solidW, worldMaterial, "World");
  world = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "World", 0, false, 0);
  //
  // Check volume
  //
  G4Box* solidC = new G4Box("Check", checkXY, checkXY, checkZ);
  logicCheck = new G4LogicalVolume(solidC, worldMaterial, "Check");
  new G4PVPlacement(0, G4ThreeVector(), logicCheck, "Check", logicWorld, false, 0);
  logicCheck->SetSensitiveDetector(checkSD);

  //
  // Target volume
  //
  G4Box* solidA = new G4Box("Target", width / 2.0, width / 2.0, targetZ);
  logicTarget = new G4LogicalVolume(solidA, targetMaterial, "Target");
  new G4PVPlacement(0, G4ThreeVector(), logicTarget, "Target", logicCheck, false, 0);

  G4cout << "### Target is box made"
         << " of " << targetMaterial->GetName() << " with Width(mm)= " << width / mm
         << "  total Length(mm)= " << 2.0 * targetZ / mm << "  ###" << G4endl;

  // colors
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* regWcolor = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  logicCheck->SetVisAttributes(regWcolor);

  G4VisAttributes* regCcolor = new G4VisAttributes(G4Colour(0., 0.3, 0.7));
  logicTarget->SetVisAttributes(regCcolor);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material != nullptr && material != targetMaterial)
  {
    targetMaterial = material;
    if (logicTarget != nullptr) logicTarget->SetMaterial(targetMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& mat)
{
  // search the material by its name
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  if (material != nullptr && material != worldMaterial)
  {
    worldMaterial = material;
    if (logicWorld != nullptr) logicWorld->SetMaterial(worldMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetWidth(G4double val)
{
  if (val > 0.0)
  {
    width = val;
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
