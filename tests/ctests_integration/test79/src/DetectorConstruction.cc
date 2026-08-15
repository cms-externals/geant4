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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4Orb.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4VisAttributes.hh"

#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
{
  /*parameters of the following experiment:
    Only channeling: A. Mazzolari et al. Phys. Rev. Lett. 112, 135503 (2014)
    Radition:        L. Bandiera et al. Phys. Rev. Lett. 115, 025504 (2015)
    Published experimental validation of G4ChannelingFastSimModel (only channeling):
    A. Sytov et al. https://arxiv.org/abs/2303.04385 (accepted for publication in
    Journal of Korean Physical Society)
  */

  // Crystal material
  CrystalMaterialStr = "G4_Si";

  // Crystal size
  CrystalSize.setX(20 * mm);
  CrystalSize.setY(20 * mm);
  CrystalSize.setZ(0.0305 * mm);

  // Crystal bending angle
  BendingAngle = 0.905 * mrad;

  // Crystal planes or axes considered
  Lattice = "(111)";

  // Crystal rotation angle (also the angle of crystal planes vs the beam)
  AngleX = 0. * 1e-6;  // rad

  // Boolean variable to activate radiation
  ActivateRadiationModel = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Check overlap option
  G4bool checkOverlaps = true;

  // Materials
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");

  // World
  G4Box* solidWorld = new G4Box("World", 1. * m, 1. * m, 20. * m);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");
  logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VPhysicalVolume* physWorld = new G4PVPlacement(0,  // no rotation
                                                   G4ThreeVector(),  // centre position
                                                   logicWorld,  // its logical volume
                                                   "World",  // its name
                                                   0,  // its mother volume
                                                   false,  // no boolean operation
                                                   0,  // copy number
                                                   checkOverlaps);  // overlaps checking

  // --------------- Crystal ------------------------------------
  // Select crystal material

  crystalMaterial = nist->FindOrBuildMaterial(CrystalMaterialStr);

  // Setting crystal rotation angle (also the angle of crystal planes vs the beam)
  G4RotationMatrix* crystalRotationMatrix = new G4RotationMatrix;
  crystalRotationMatrix->rotateY(-AngleX);

  // Setting crystal position
  G4ThreeVector posCrystal = G4ThreeVector(0. * mm, 0. * mm, CrystalSize.z() / 2.);

  // crystal volume
  G4Box* crystalSolid =
    new G4Box("Crystal", CrystalSize.x() / 2, CrystalSize.y() / 2, CrystalSize.z() / 2.);

  crystalLogic = new G4LogicalVolume(crystalSolid, crystalMaterial, "Crystal");

  // visualization attributes
  G4VisAttributes* CrystalVisAttribute = new G4VisAttributes(G4Colour(0., 0., 1.));
  CrystalVisAttribute->SetForceSolid(true);
  crystalLogic->SetVisAttributes(CrystalVisAttribute);

  new G4PVPlacement(crystalRotationMatrix, posCrystal, crystalLogic, "Crystal", logicWorld, false,
                    0, checkOverlaps);

  // crystal region (necessary for the FastSim model)
  fRegion = new G4Region("Crystal");
  fRegion->AddRootLogicalVolume(crystalLogic);

  // print crystal info
  G4cout << "Crystal size: " << CrystalSize.x() / mm << "x" << CrystalSize.y() / mm << "x"
         << CrystalSize.z() / mm << " mm3" << G4endl;
  G4cout << "Crystal bending angle: " << BendingAngle << " rad" << G4endl;
  G4cout << "Crystal AngleX: " << AngleX << " rad" << G4endl;
  G4cout << "ActivateRadiationModel: " << ActivateRadiationModel << G4endl;

  // --------------- Detector -----------------------------------
  // Setting detector position
  G4ThreeVector posDetector = G4ThreeVector(0, 0, 5973 * mm);

  // particle detector volume
  G4Box* Detector = new G4Box("Detector", 20 * cm / 2, 20 * cm / 2, 0.3 * mm / 2);

  G4LogicalVolume* fLogicDetector = new G4LogicalVolume(Detector, Silicon, "Detector");
  new G4PVPlacement(0, posDetector, fLogicDetector, "Detector", logicWorld, false, 0,
                    checkOverlaps);

  // always return the physical World
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  // --------------- fast simulation ----------------------------
  // extract the region of the crystal from the store
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  G4Region* RegionCh = regionStore->GetRegion("Crystal");

  // create the channeling model for this region
  G4ChannelingFastSimModel* ChannelingModel =
    new G4ChannelingFastSimModel("ChannelingModel", RegionCh);
  // activate the channeling model
  ChannelingModel->Input(crystalMaterial, Lattice);
  // setting bending angle of the crystal planes (default is 0)
  ChannelingModel->GetCrystalData()->SetBendingAngle(BendingAngle, crystalLogic);

  // activate radiation model (do it only when you want to take into account the
  // radiation production in an oriented crystal; it takes a lot of computational power)
  if (ActivateRadiationModel) ChannelingModel->RadiationModelActivate();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
