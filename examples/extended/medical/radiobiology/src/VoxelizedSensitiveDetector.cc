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
/// \file radiobiology/src/VoxelizedSensitiveDetector.cc
/// \brief Implementation of the RadioBio::VoxelizedSensitiveDetector class

#include "VoxelizedSensitiveDetector.hh"

#include "G4Box.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"

#include "DetectorConstruction.hh"
#include "SD.hh"
#include "VoxelizedSensitiveDetectorMessenger.hh"

#include <algorithm>
#include <limits>

namespace
{
G4int CalculateVoxelCount(G4double containerWidth, G4double voxelWidth, const G4String& axis)
{
  const G4double count = containerWidth / voxelWidth;
  if (count > std::numeric_limits<G4int>::max()) {
    G4ExceptionDescription description;
    description << "The requested voxel count along the " << axis
                << " axis exceeds the supported integer range.";
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "TooManyVoxels", FatalException,
                description);
    return 1;
  }
  return std::max(1, static_cast<G4int>(count));
}
}  // namespace

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector* VoxelizedSensitiveDetector::CreateInstance(DetectorConstruction* det,
                                                                       G4double xWidth,
                                                                       G4double yWidth,
                                                                       G4double zWidth)
{
  if (fInstance != nullptr) {
    G4Exception("VoxelizedSensitiveDetector::createInstance", "RecreatingVoxelization",
                FatalException, "Creating another, new, instance of VoxelizedSensitiveDetector");
    return fInstance;
  }
  fInstance = new VoxelizedSensitiveDetector(det, xWidth, yWidth, zWidth);
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector* VoxelizedSensitiveDetector::GetInstance()
{
  return fInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::DestroyInstance()
{
  delete fInstance;
  fInstance = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector::VoxelizedSensitiveDetector(DetectorConstruction* det, G4double xWidth,
                                                       G4double yWidth, G4double zWidth)
  : fDetector(det), fVoxelWidthX(xWidth), fVoxelWidthY(yWidth), fVoxelWidthZ(zWidth)
{
  fVoxelizedSensitiveDetectorMessenger = new VoxelizedSensitiveDetectorMessenger(this);
  UpdateVoxelVolume();
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

VoxelizedSensitiveDetector::~VoxelizedSensitiveDetector()
{
  delete fVoxelizedSensitiveDetectorMessenger;
  fInstance = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::UpdateVoxelVolume()
{
  fVoxelVolume = fVoxelWidthX * fVoxelWidthY * fVoxelWidthZ;
  fVoxelDensity = fDetector->GetMaterial()->GetDensity();
  fVoxelMass = fVoxelVolume * fVoxelDensity;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidth(const G4ThreeVector& voxWidth)
{
  fVoxelWidthX = voxWidth.getX();
  fVoxelWidthY = voxWidth.getY();
  fVoxelWidthZ = voxWidth.getZ();
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthX(G4double voxWidthX)
{
  if (fVoxelWidthX == voxWidthX) return;
  fVoxelWidthX = voxWidthX;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthY(G4double voxWidthY)
{
  if (fVoxelWidthY == voxWidthY) return;
  fVoxelWidthY = voxWidthY;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::SetVoxelWidthZ(G4double voxWidthZ)
{
  if (fVoxelWidthZ == voxWidthZ) return;
  fVoxelWidthZ = voxWidthZ;
  CalculateVoxelNumber();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Calculate the number of voxels, approximating to an integer number of voxels.
// Then recalculates voxels size according to approximations

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::CalculateVoxelNumber()
{
  if (fDetector == nullptr || fDetector->GetSizeX() <= 0. || fDetector->GetSizeY() <= 0.
      || fDetector->GetSizeZ() <= 0. || fVoxelWidthX <= 0. || fVoxelWidthY <= 0.
      || fVoxelWidthZ <= 0.)
  {
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "InvalidVoxelDimensions",
                FatalException, "Container and voxel dimensions must be strictly positive.");
    return;
  }

  // A requested width larger than the container means one voxel along that axis.
  // Clamping the counts also prevents a division by zero while individual UI
  // commands update X, Y, and Z sequentially.
  fVoxelNumberAlongX = CalculateVoxelCount(fDetector->GetSizeX(), fVoxelWidthX, "X");
  fVoxelWidthX = fDetector->GetSizeX() / G4double(fVoxelNumberAlongX);

  fVoxelNumberAlongY = CalculateVoxelCount(fDetector->GetSizeY(), fVoxelWidthY, "Y");
  fVoxelWidthY = fDetector->GetSizeY() / G4double(fVoxelNumberAlongY);

  fVoxelNumberAlongZ = CalculateVoxelCount(fDetector->GetSizeZ(), fVoxelWidthZ, "Z");
  fVoxelWidthZ = fDetector->GetSizeZ() / G4double(fVoxelNumberAlongZ);

  if (fVoxelNumberAlongY % 2 == 0) {
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "VoxelNumberYEven", JustWarning,
                "Trying to voxelize with an even number of voxels along the Y axis. "
                "Please select an odd number to prevent from warnings due to tracking");
  }

  if (fVoxelNumberAlongZ % 2 == 0) {
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "VoxelNumberZEven", JustWarning,
                "Trying to voxelize with an even number of voxels along the Z axis. "
                "Please select an odd number to prevent from warnings due to tracking");
  }

  const G4long voxelNumberXY =
    static_cast<G4long>(fVoxelNumberAlongX) * static_cast<G4long>(fVoxelNumberAlongY);
  const G4long maximumVoxelNumber = std::numeric_limits<G4int>::max();
  if (voxelNumberXY > maximumVoxelNumber / fVoxelNumberAlongZ) {
    G4Exception("VoxelizedSensitiveDetector::CalculateVoxelNumber", "TooManyVoxels", FatalException,
                "The total voxel count exceeds the supported integer range.");
    return;
  }
  const G4long totalVoxelNumber = voxelNumberXY * fVoxelNumberAlongZ;
  fTotalVoxelNumber = static_cast<G4int>(totalVoxelNumber);

  UpdateVoxelVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructXDivision()
{
  if (fWorldLogical == nullptr) {
    G4Exception("VoxelizedSensitiveDetector::ConstructXDivision", "WorldNotInit", FatalException,
                "Voxelizing without having a pointer to world logical volume!");
    return;
  }

  if (fDetector == nullptr) {
    G4Exception("VoxelizedSensitiveDetector::ConstructXDivision", "DetConstInit", FatalException,
                "Voxelizing without having a pointer to DetectorConstruction!");
    return;
  }

  fVoxelizedDetectorXDivision = new G4Box("VoxelizedDetectorXDivision", fVoxelWidthX / 2,
                                          fDetector->GetSizeY() / 2, fDetector->GetSizeZ() / 2);

  fVoxelizedDetectorXDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorXDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorXDivisionLog", 0, 0, 0);

  fVoxelizedDetectorXDivisionPhys =
    new G4PVReplica("VoxelizedDetectorXDivisionPhys", fVoxelizedDetectorXDivisionLog, fWorldLogical,
                    kXAxis, fVoxelNumberAlongX, fVoxelWidthX);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructYDivision()
{
  fVoxelizedDetectorYDivision = new G4Box("VoxelizedDetectorYDivision", fVoxelWidthX / 2,
                                          fVoxelWidthY / 2, fDetector->GetSizeZ() / 2);

  fVoxelizedDetectorYDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorYDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorYDivisionLog", 0, 0, 0);

  fVoxelizedDetectorYDivisionPhys =
    new G4PVReplica("VoxelizedDetectorYDivisionPhys", fVoxelizedDetectorYDivisionLog,
                    fVoxelizedDetectorXDivisionLog, kYAxis, fVoxelNumberAlongY, fVoxelWidthY);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructZDivision()
{
  fVoxelizedDetectorZDivision =
    new G4Box("VoxelizedDetectorZDivision", fVoxelWidthX / 2, fVoxelWidthY / 2, fVoxelWidthZ / 2);

  fVoxelizedDetectorZDivisionLog =
    new G4LogicalVolume(fVoxelizedDetectorZDivision, fWorldLogical->GetMaterial(),
                        "VoxelizedDetectorZDivisionLog", 0, 0, 0);

  fVoxelizedDetectorZDivisionPhys =
    new G4PVReplica("VoxelizedDetectorZDivisionPhys", fVoxelizedDetectorZDivisionLog,
                    fVoxelizedDetectorYDivisionPhys, kZAxis, fVoxelNumberAlongZ, fVoxelWidthZ);

  fSensitiveLogicalVolume = fVoxelizedDetectorZDivisionLog;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// First voxelize along X, then Y, then Z

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool VoxelizedSensitiveDetector::ConstructVoxelizedDetector()
{
  // Creating X division
  ConstructXDivision();

  // Creating Y division
  ConstructYDivision();

  // Creating Z division
  ConstructZDivision();

  // Set last, smallest volumes as sensitive
  fSensitiveLogicalVolume = fVoxelizedDetectorZDivisionLog;
  fIsBuilt = true;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::UpdateVoxelizedGeometry()
{
  CalculateVoxelNumber();

  // Nothing else happens if the geometry is not built. The updated parameters
  // will be used during the first construction.
  if (!fIsBuilt) {
    return;
  }

  // Rebuild the complete detector. DetectorConstruction cleans the Geant4
  // volume stores and recreates the solids, logical volumes, and replicas with
  // mutually consistent dimensions.
  fIsBuilt = false;
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::ConstructSD()
{
  const G4String sensitiveDetectorName = "VoxelizedDetector";
  const G4String hitsCollectionName = "LETdata";

  auto* sensitiveDetectorManager = G4SDManager::GetSDMpointer();
  auto* detectorSD = sensitiveDetectorManager->FindSensitiveDetector(sensitiveDetectorName, false);
  if (detectorSD == nullptr) {
    detectorSD = new SD(sensitiveDetectorName, hitsCollectionName);
    sensitiveDetectorManager->AddNewDetector(detectorSD);
  }
  fSensitiveLogicalVolume->SetSensitiveDetector(detectorSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::Construct()
{
  CalculateVoxelNumber();
  ConstructVoxelizedDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void VoxelizedSensitiveDetector::InitializeWorldPtr(G4VPhysicalVolume* pWorld)
{
  if (pWorld == nullptr) {
    G4Exception("VoxelizedSensitiveDetector::InitializeWorldPtr", "WorldinitNull", FatalException,
                "Initializing Voxelization Class with a Null Pointer to World!");
    return;
  }
  fWorldLogical = pWorld->GetLogicalVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Returns absolute voxel index given matrix indexes

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int VoxelizedSensitiveDetector::GetThisVoxelNumber(G4int x, G4int y, G4int z) const
{
  G4int nz = GetVoxelNumberAlongZ();
  G4int ny = GetVoxelNumberAlongY();

  return z + nz * (y + ny * (x));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio
