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
// --------------------------------------------------------------------
#include <assert.h>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4BoundingEnvelope.hh"
#include "G4GeometryTolerance.hh"

/////////////////////////////////////////////////////////////
//
//  Main

int main() 
{ 
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double delta = kCarTolerance;
  G4double EPS   = 1e-12;
  G4cout.precision(16);

  //
  // Create bounding envelope by 4 polygons (3 prisms)
  //
  //                          |
  //                  ........+........
  //                    .     |         .
  //                      ....+....       . 
  //                          |     .       .
  //                          + 20    .       . 
  //                          |         .     .
  //                          + 10        .   .
  //                 -20 -10  |           .   .
  //        --+---+---+---+---o---+---+---+---+--
  //                          |   10  20  .   .
  //                      -10 +           .   .
  //                          |             . .
  //                      -20 +               .
  //                          |
  //

  std::vector<const G4ThreeVectorList*> polygons(4);
  G4ThreeVectorList* pppp;
  pppp = new G4ThreeVectorList(0);
  pppp->push_back(G4ThreeVector(30,-10,-13));
  pppp->push_back(G4ThreeVector(40,-20,-13));
  pppp->push_back(G4ThreeVector(40,-20, 13));
  pppp->push_back(G4ThreeVector(30,-10, 13));
  polygons[0] = pppp;

  pppp = new G4ThreeVectorList(0);
  pppp->push_back(G4ThreeVector(30, 10,-13));
  pppp->push_back(G4ThreeVector(40, 20,-13));
  pppp->push_back(G4ThreeVector(40, 20, 13));
  pppp->push_back(G4ThreeVector(30, 10, 13));
  polygons[1] = pppp;

  pppp = new G4ThreeVectorList(0);
  pppp->push_back(G4ThreeVector(10, 30,-13));
  pppp->push_back(G4ThreeVector(20, 40,-13));
  pppp->push_back(G4ThreeVector(20, 40, 13));
  pppp->push_back(G4ThreeVector(10, 30, 13));
  polygons[2] = pppp;

  pppp = new G4ThreeVectorList(0);
  pppp->push_back(G4ThreeVector(-10, 30,-13));
  pppp->push_back(G4ThreeVector(-20, 40,-13));
  pppp->push_back(G4ThreeVector(-20, 40, 13));
  pppp->push_back(G4ThreeVector(-10, 30, 13));
  polygons[3] = pppp;

  G4BoundingEnvelope envelope(polygons);

  // Various tests
  //
  G4Transform3D transform;
  G4VoxelLimits limits;
  G4double pmin, pmax;
  G4bool reply, exist;

  ///////////////////////////////////////////////////////////
  //
  // Unrotated envelope
  //
  ///////////////////////////////////////////////////////////

  G4cout << "\nUnrotated envelope\n" << G4endl;

  // Unlimited voxel
  limits = G4VoxelLimits();
  reply = envelope.BoundingBoxVsVoxelLimits(kXAxis,limits,transform, pmin, pmax);
  exist = false; if (reply) exist = (pmin < pmax) ? true : false;
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == -20-delta); assert(pmax == 40+delta);

  // Voxel outside the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis,-40,-30);
  reply = envelope.BoundingBoxVsVoxelLimits(kXAxis,limits,transform, pmin, pmax);
  exist = true; if (reply) exist = (pmin < pmax) ? true : false;
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == false); assert(pmin == kInfinity); assert(pmax == -kInfinity);

  // Voxel intersect the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kYAxis,12,15);
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == 12-kCarTolerance); assert(pmax == 15+kCarTolerance);

  // Voxel intersect the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kYAxis, 5,20);
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(std::abs(pmin-(20-2*delta))<EPS); assert(pmax == 40+delta);

  // Voxel inside the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis, 1,5); limits.AddLimit(kYAxis, 31,35);
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(std::abs(pmin-(31-delta))<EPS); assert(pmax == 35+delta);

  // Touch the envelope from outside at Y=30
  limits = G4VoxelLimits(); limits.AddLimit(kYAxis, 5, 30-0.5*delta);
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(std::abs(pmin-(-10-1.5*delta))<EPS); assert(pmax == 40+delta);

  // Touch the envelope from inside at Y=40
  limits = G4VoxelLimits(); limits.AddLimit(kYAxis, 5, 40-0.5*delta);
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == -20-delta); assert(pmax == 40+delta);

  ///////////////////////////////////////////////////////////
  //
  // Rotated envelope
  //
  ///////////////////////////////////////////////////////////

  G4cout << "\nRotated envelope\n" << G4endl;

  transform = G4RotateZ3D(-pi/4);
  G4Point3D p1 = transform*G4Point3D(40,-20,0); 
  G4Point3D p2 = transform*G4Point3D(40, 20,0); 
  G4Point3D p3 = (p1+p2)/2;

  // Unlimited voxel
  limits = G4VoxelLimits();
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true);
  assert(std::abs(pmin-(p1.x()-delta))<EPS);
  assert(std::abs(pmax-(p2.x()+delta))<EPS);

  // Voxel outside the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis,60,70);
  reply = envelope.BoundingBoxVsVoxelLimits(kXAxis,limits,transform, pmin, pmax);
  exist = true; if (reply) exist = (pmin < pmax) ? true : false;
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == false); assert(pmin == kInfinity); assert(pmax == -kInfinity);

  // Voxel inside the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kZAxis,-10,10);
  exist = envelope.CalculateExtent(kZAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == -10-kCarTolerance); assert(pmax == 10+kCarTolerance);

  // Voxel inside the envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis,30,40);
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == 30-kCarTolerance); assert(pmax == 40+kCarTolerance);

  // Limited voxel completely inside the envelope

  limits = G4VoxelLimits();
  limits.AddLimit(kXAxis, 32,38); 
  limits.AddLimit(kYAxis,-16,16);
  limits.AddLimit(kZAxis,-12,11);
  exist = envelope.CalculateExtent(kXAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == 32-kCarTolerance); assert(pmax == 38+kCarTolerance);
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == -16-kCarTolerance); assert(pmax == 16+kCarTolerance);
  exist = envelope.CalculateExtent(kZAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true); assert(pmin == -12-kCarTolerance); assert(pmax == 11+kCarTolerance);

  // Voxel intesect the two branches of the rotated envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis, 0,25);
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true);
  assert(std::abs(pmin-(p1.y()-delta))<EPS);
  assert(std::abs(pmax+(p1.y()-delta))<EPS);

  // Voxel intesect the two branches of the rotated envelope
  limits = G4VoxelLimits(); limits.AddLimit(kXAxis, p3.x()+delta, p3.x()+1);
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true);
  assert(std::abs(pmin-(p3.y()-delta))<EPS);
  assert(std::abs(pmax+(p3.y()-delta))<EPS);

  // Litle touch from outside
  limits = G4VoxelLimits(); limits.AddLimit(kYAxis, -100, p1.y());
  exist = envelope.CalculateExtent(kYAxis,limits,transform, pmin, pmax);
  G4cout << "Reply = " << (exist?"True":"False") << "\t min, max = " << pmin << ", " << pmax << G4endl;
  assert(exist == true);
  assert(std::abs(pmin-(p1.y()-delta))<EPS);
  assert(std::abs(pmax-(p1.y()+kCarTolerance))<EPS);

  return 0;     
}
