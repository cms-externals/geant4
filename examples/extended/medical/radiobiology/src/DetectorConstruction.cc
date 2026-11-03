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
/// \file radiobiology/src/DetectorConstruction.cc
/// \brief Implementation of the RadioBio::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SolidStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "DetectorMessenger.hh"
#include "VoxelizedSensitiveDetector.hh"

#include <sstream>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  auto* nist = G4NistManager::Instance();

  // The voxelized container is water by default.
  fMaterial = nist->FindOrBuildMaterial("G4_WATER");

  // The outer world is air.
  fWorldMaterial = nist->FindOrBuildMaterial("G4_AIR");

  fDetectorMessenger = new DetectorMessenger(this);

  // Default voxel dimensions refer to the inner container.
  VoxelizedSensitiveDetector::CreateInstance(this, 0.1 * m, 1. * m, 1. * m);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  VoxelizedSensitiveDetector::DestroyInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  VoxelizedSensitiveDetector::GetInstance()->ConstructSD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::CheckGeometryDimensions() const
{
  if (fWorldSizeX <= fBoxSizeX || fWorldSizeY <= fBoxSizeY || fWorldSizeZ <= fBoxSizeZ) {
    std::ostringstream message;
    message << "The air world must be larger than the voxelized container along every axis.\n"
            << "World size:    (" << G4BestUnit(fWorldSizeX, "Length") << ", "
            << G4BestUnit(fWorldSizeY, "Length") << ", " << G4BestUnit(fWorldSizeZ, "Length")
            << ")\n"
            << "Container size: (" << G4BestUnit(fBoxSizeX, "Length") << ", "
            << G4BestUnit(fBoxSizeY, "Length") << ", " << G4BestUnit(fBoxSizeZ, "Length") << ")";

    G4Exception("DetectorConstruction::CheckGeometryDimensions", "WorldTooSmall", FatalException,
                message.str().c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  CheckGeometryDimensions();

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // --------------------------------------------------------------------------
  // Air world
  // --------------------------------------------------------------------------
  auto* worldSolid = new G4Box("World", fWorldSizeX / 2., fWorldSizeY / 2., fWorldSizeZ / 2.);

  fLWorld = new G4LogicalVolume(worldSolid, fWorldMaterial, "WorldLogical");

  fPWorld =
    new G4PVPlacement(nullptr, G4ThreeVector(), fLWorld, "WorldPhysical", nullptr, false, 0, false);

  // --------------------------------------------------------------------------
  // Voxelized material container
  // --------------------------------------------------------------------------
  auto* boxSolid = new G4Box("Container", fBoxSizeX / 2., fBoxSizeY / 2., fBoxSizeZ / 2.);

  fLBox = new G4LogicalVolume(boxSolid, fMaterial, "ContainerLogical");

  fPBox = new G4PVPlacement(nullptr, G4ThreeVector(), fLBox, "ContainerPhysical", fLWorld, false, 0,
                            false);

  // The voxelization must act on the inner material container, not on the world.
  VoxelizedSensitiveDetector::GetInstance()->InitializeWorldPtr(fPBox);
  VoxelizedSensitiveDetector::GetInstance()->Construct();

  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\nAir world dimensions:" << G4endl << "  x: " << G4BestUnit(fWorldSizeX, "Length")
         << G4endl << "  y: " << G4BestUnit(fWorldSizeY, "Length") << G4endl
         << "  z: " << G4BestUnit(fWorldSizeZ, "Length") << G4endl
         << "  material: " << fWorldMaterial->GetName() << G4endl;

  G4cout << "\nVoxelized container dimensions:" << G4endl
         << "  x: " << G4BestUnit(fBoxSizeX, "Length") << G4endl
         << "  y: " << G4BestUnit(fBoxSizeY, "Length") << G4endl
         << "  z: " << G4BestUnit(fBoxSizeZ, "Length") << G4endl
         << "  material: " << fMaterial->GetName() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (material != nullptr) {
    if (fMaterial != material) {
      fMaterial = material;
      if (fLBox != nullptr) {
        VoxelizedSensitiveDetector::GetInstance()->UpdateVoxelizedGeometry();
      }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
  else {
    std::ostringstream message;
    message << "Material " << materialChoice << " does not exist; keeping material "
            << fMaterial->GetName();

    G4Exception("DetectorConstruction::SetMaterial", "NoContainerMaterial", JustWarning,
                message.str().c_str());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  SetSizeX(value);
  SetSizeY(value);
  SetSizeZ(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(const G4ThreeVector& size)
{
  SetSizeX(size.x());
  SetSizeY(size.y());
  SetSizeZ(size.z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeX(G4double sizeX)
{
  fBoxSizeX = sizeX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeY(G4double sizeY)
{
  fBoxSizeY = sizeY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSizeZ(G4double sizeZ)
{
  fBoxSizeZ = sizeZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSize(G4double value)
{
  SetWorldSizeX(value);
  SetWorldSizeY(value);
  SetWorldSizeZ(value);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSize(const G4ThreeVector& size)
{
  SetWorldSizeX(size.x());
  SetWorldSizeY(size.y());
  SetWorldSizeZ(size.z());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeX(G4double sizeX)
{
  fWorldSizeX = sizeX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeY(G4double sizeY)
{
  fWorldSizeY = sizeY;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldSizeZ(G4double sizeZ)
{
  fWorldSizeZ = sizeZ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
