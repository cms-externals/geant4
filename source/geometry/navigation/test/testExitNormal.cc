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
// Geant4 includes
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Navigator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

// geometrical constants
const double topSize = 10.;
const double coneRadius = 4.;
const double coneHeight = 4.;
const double phi = -45.;
const double theta = 90.;
const double psi = 0.;

G4VPhysicalVolume *g4geometry() {
  // material
  G4Material *waterMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Water");

  // geometry
  const G4bool fCheckOverlaps = true;
  G4Box *solidTOP = new G4Box("TOP", 0.5*topSize*cm, 0.5*topSize*cm, 0.5*topSize*cm);

  G4LogicalVolume *logicTOP =
    new G4LogicalVolume(solidTOP,      // its solid
                        waterMat,      // its material
                        "TOP");        // its name

  G4VPhysicalVolume *physTOP =
    new G4PVPlacement(0,                     // no rotation
                      G4ThreeVector(),       // at (0,0,0)
                      logicTOP,              // its logical volume
                      "TOP",                 // its name
                      0,                     // its mother  volume
                      false,                 // no boolean operation
                      0,                     // copy number
                      fCheckOverlaps);       // checking overlaps

  G4Cons *solidCone = new G4Cons("cone",
                                 0.*cm, coneRadius*cm,     // rMin1, rMax1
                                 0.*cm, 0.*cm,             // rMin2, rMax2
                                 0.5*coneHeight*cm,        // dZ
                                 0., twopi);               // phiMin, deltaPhi

  G4LogicalVolume *logicCone =
    new G4LogicalVolume(solidCone,      // its solid
                        waterMat,       // its material
                        "cone");        // its name

  new G4PVPlacement(new G4RotationMatrix(phi*degree, theta*degree, psi*degree),
                    G4ThreeVector(),       // at (0,0,0)
                    logicCone,            // its logical volume
                    "cone",               // its name
                    logicTOP,              // its mother  volume
                    false,                 // no boolean operation
                    0,                     // copy number
                    fCheckOverlaps);       // checking overlaps
  return physTOP;
}

void g4test(G4double const x, G4double const y, G4double const z,
            G4double const vx, G4double const vy, G4double const vz)
{
  G4Navigator *navigator = new G4Navigator;
  G4VPhysicalVolume *world = g4geometry();
  navigator->SetWorldVolume(world);

  G4ThreeVector point(x, y, z);
  G4ThreeVector omega(vx, vy, vz);

  G4VPhysicalVolume* vol = navigator->LocateGlobalPointAndSetup(point);
  G4double safety = 0.0;
  G4double dist = navigator->ComputeStep(point, omega, kInfinity, safety);
  G4bool valid = false;
  G4ThreeVector new_point = point + dist*omega;
//  G4VPhysicalVolume* new_vol = navigator->LocateGlobalPointAndSetup(new_point);
  G4ThreeVector normal = navigator->GetGlobalExitNormal(new_point, &valid);
  G4cout.precision(16);
  G4cout << "Geant4 TEST\n"
    << " first volume = " << vol->GetName() << '\n'
//    << " last volume = " << new_vol->GetName() << '\n'
    << " point = " << point << '\n'
    << " new point = " << new_point << '\n'
    << " dist = " << dist/cm << '\n'
    << " normal = " << normal << '\n'
    << " valid = " << valid << G4endl;
}

int main()
{
  const G4double x = -1.331514e+00*cm, y=-9.139028e-01*cm, z=-2.876654e-01*cm;
  const G4double vx = 4.215310e-01, vy=8.672989e-01, vz=2.647721e-01;
  g4test(x, y, z, vx, vy, vz);
  return 0;
}
