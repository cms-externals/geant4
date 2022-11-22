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
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"

#include "G4ExtrudedSolid.hh"
#include "G4Polyhedra.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

////////////////////////////////////////////////////////////////////////
//
// Declare auxiliary routines

G4ThreeVector RandomDirection(G4double cosTheta);

////////////////////////////////////////////////////////////////////////
//
// Declare auxiliary routines

G4ThreeVector PointOnPrism(G4int nsect,
                           G4double sphi, G4double dphi,
                           G4double rmax, G4double rmin, G4double dz);
G4TessellatedSolid* CreatePrism(G4int nsect,
                                G4double sphi, G4double dphi,
                                G4double rmax, G4double rmin, G4double dz);

////////////////////////////////////////////////////////////////////////
//
// Declare bench routines

void Check_Inside       (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SurfaceNormal(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP);
void Check_SafetyToIn   (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP);
void Check_SafetyToOut  (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP);
void Check_DistanceToIn (G4VSolid* t1, const G4ThreeVector& p,
                         const std::vector<G4ThreeVector>& vv, G4int NLOOP);
void Check_DistanceToOut(G4VSolid* t1,
                         const std::vector<G4ThreeVector>& pp,
                         const std::vector<G4ThreeVector>& vv, G4bool calcNorm, G4int NLOOP);

////////////////////////////////////////////////////////////////////////
//
// Main

int main()
{
  G4cout << "*********************************************************************" << G4endl;
  G4cout << "*                                                                   *" << G4endl;
  G4cout << "* This benchmark measures the performance of the navigation methods *" << G4endl;
  G4cout << "* for concave and convex right prisms defined as G4ExtrudedSolid,   *" << G4endl;
  G4cout << "* G4TesselatedSolid and G4Polyhedra                                 *" << G4endl;
  G4cout << "*                                                                   *" << G4endl;
  G4cout << "*          Convex prism             Non-convex prism                *" << G4endl;
  G4cout << "*        (hexagonal prism)   (half of a hollow hexagonal prism)     *" << G4endl;
  G4cout << "*                                                                   *" << G4endl;
  G4cout << "*               .                              .                    *" << G4endl;
  G4cout << "*            .......                        ...|                    *" << G4endl;
  G4cout << "*         .............                  ......|                    *" << G4endl;
  G4cout << "*       |...............|              |......                      *" << G4endl;
  G4cout << "*       |...............|              |...|                        *" << G4endl;
  G4cout << "*       |...............|              |...|                        *" << G4endl;
  G4cout << "*       |...............|              |......                      *" << G4endl;
  G4cout << "*         .............                   .....|                    *" << G4endl;
  G4cout << "*            .......                        ...|                    *" << G4endl;
  G4cout << "*               .                              .                    *" << G4endl;
  G4cout << "*                                                                   *" << G4endl;
  G4cout << "*********************************************************************" << G4endl;

  // Set number of calls
  G4int NLOOP = 10000000;
  G4int NMIL  = NLOOP/1000000;
  G4cout << "\n       Number of calls : "<< NMIL << " 000 000" << G4endl;

  // Select solids to test: 0 - omit, 1 - test
  G4int ifxtru1 = 1, iftess1 = 1, ifpoly1 = 1; // convex prisms
  G4int ifxtru2 = 1, iftess2 = 1, ifpoly2 = 1; // non-convex prisms

  // Define number of random points
  const G4int NP = 100000;

  // Set auxiliary coeffifients
  G4double ksur = 1e-12, kin = 0.5, kout = 2.;

  // Generate random points for convex (hexagonal) prism
  G4int nsect1;
  G4double sphi1,dphi1,rmax1,rmin1,dz1;
  std::vector<G4ThreeVector> prism1_Surface(NP);
  std::vector<G4ThreeVector> prism1_Inside(NP);
  std::vector<G4ThreeVector> prism1_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnPrism(nsect1=6, sphi1=30*deg, dphi1=360*deg, rmax1=25*cm, rmin1=0, dz1=20*cm);
    prism1_Surface[i] = (1.0 + k*ksur)*p;
    prism1_Inside [i] = kin  * p;
    prism1_Outside[i] = kout * p;
  }
  PointOnPrism(0,0,0,0,0,0);

  // Generate random points for non-convex (half of hollow hexagonal) prism
  G4int nsect2;
  G4double sphi2,dphi2,rmax2,rmin2,dz2;
  std::vector<G4ThreeVector> prism2_Surface(NP);
  std::vector<G4ThreeVector> prism2_Inside(NP);
  std::vector<G4ThreeVector> prism2_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnPrism(nsect2=3, sphi2=90*deg, dphi2=180*deg, rmax2=25*cm, rmin2=12.5*cm, dz2=20*cm);
    prism2_Surface[i] = (1.0 + k*ksur)*p;
    prism2_Inside [i] = PointOnPrism(nsect2, sphi2+5*deg, dphi2-10*deg, rmax2-5*cm, rmin2+4*cm, dz2-2*cm);
    prism2_Outside[i] = PointOnPrism(nsect2, sphi2-5*deg, dphi2+10*deg, rmax2+5*cm, rmin2-4*cm, dz2+2*cm);
  }
  PointOnPrism(0,0,0,0,0,0);

  // ========= Create convex (hexagonal) prisms using different solids
  //
  std::vector<G4TwoVector> polygon1(nsect1);
  G4double ang1 = dphi1/nsect1;
  for (G4int i=0; i<nsect1; ++i)
  {
    G4double phi = i*ang1 + sphi1;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    polygon1[i].set(rmax1*cosphi,rmax1*sinphi);
  }

  // G4ExtrudedSolid
  G4TwoVector offA(0,0), offB(0,0);
  G4double scaleA = 1, scaleB = 1;
  G4VSolid* xtru1 = new G4ExtrudedSolid("Extru 1", polygon1, dz1, offA, scaleA, offB, scaleB);

  // G4TesselatedSold
  G4VSolid* tess1 = CreatePrism(nsect1,sphi1,dphi1,rmax1,rmin1,dz1);

  // G4Polyhedra
  const G4int nrz1 = 4;
  G4double zz1[nrz1] = {-dz1,-dz1, dz1, dz1}, rr1[nrz1] = { 0, rmax1, rmax1, 0 };
  G4VSolid* poly1 = new G4Polyhedra("Polyhedra 1", sphi1, dphi1, nsect1, nrz1, rr1, zz1);

  // ========= Create non-convex (half of hollow hexagonal) prisms using different solids
  //
  G4int nvert = 2*(nsect2+1);
  std::vector<G4TwoVector> polygon2(nvert);
  G4double ang2 = dphi2/nsect2;
  for (G4int i=0; i<nsect2+1; ++i)
  {
    G4double phi = i*ang2 + sphi2;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    polygon2[i].set(rmax2*cosphi,rmax2*sinphi);
    polygon2[nvert-1-i].set(rmin2*cosphi,rmin2*sinphi);
  }

  // G4ExtrudedSolid
  G4VSolid* xtru2 = new G4ExtrudedSolid("Extru 2", polygon2, dz2, offA, scaleA, offB, scaleB);

  // G4TesselatedSold
  G4VSolid* tess2 = CreatePrism(nsect2,sphi2,dphi2,rmax2,rmin2,dz2);

  // G4Polyhedra
  const G4int nrz2 = 4;
  G4double zz2[nrz2] = {-dz2,-dz2, dz2, dz2}, rr2[nrz2] = { rmin2, rmax2, rmax2, rmin2 };
  G4VSolid* poly2 = new G4Polyhedra("Polyhedra 2", sphi2, dphi2, nsect2, nrz2, rr2, zz2);

  // Set start points and generate random directions
  //
  G4ThreeVector p1(0,0,-62.2*cm);
  G4ThreeVector p2(-5*cm,0,-45.6*cm);

  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i=0; i<NP; i++) {
    vecInCone[i]   = RandomDirection(std::cos(40*deg));
    vecOnSphere[i] = RandomDirection(-1);
  }

  // Mesure time of Inside(p)
  //
  G4cout << "\n   Inside(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::Inside(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(xtru1, prism1_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(xtru1, prism1_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(xtru1, prism1_Outside, NLOOP, kOutside);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::Inside(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tess1, prism1_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(tess1, prism1_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tess1, prism1_Outside, NLOOP, kOutside);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::Inside(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(poly1, prism1_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(poly1, prism1_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(poly1, prism1_Outside, NLOOP, kOutside);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::Inside(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(xtru2, prism2_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(xtru2, prism2_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(xtru2, prism2_Outside, NLOOP, kOutside);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::Inside(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tess2, prism2_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(tess2, prism2_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tess2, prism2_Outside, NLOOP, kOutside);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::Inside(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(poly2, prism2_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(poly2, prism2_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(poly2, prism2_Outside, NLOOP, kOutside);
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::SurfaceNormal(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(xtru1, prism1_Surface, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::SurfaceNormal(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tess1, prism1_Surface, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::SurfaceNormal(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(poly1, prism1_Surface, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::SurfaceNormal(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(xtru2, prism2_Surface, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::SurfaceNormal(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tess2, prism2_Surface, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::SurfaceNormal(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(poly2, prism2_Surface, NLOOP);
  }
 
  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(xtru1, prism1_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(xtru1, prism1_Outside, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tess1, prism1_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tess1, prism1_Outside, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p) - Convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(poly1, prism1_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(poly1, prism1_Outside, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(xtru2, prism2_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(xtru2, prism2_Outside, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tess2, prism2_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tess2, prism2_Outside, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(poly2, prism2_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(poly2, prism2_Outside, NLOOP);
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(xtru1, prism1_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(xtru1, prism1_Surface, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tess1, prism1_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(tess1, prism1_Surface, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(poly1, prism1_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(poly1, prism1_Surface, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(xtru2, prism2_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(xtru2, prism2_Surface, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tess2, prism2_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(tess2, prism2_Surface, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(poly2, prism2_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(poly2, prism2_Surface, NLOOP);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p,v) - Convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(xtru1, p1, vecInCone, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p,v) - Convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tess1, p1, vecInCone, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p,v) - Convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(poly1, p1, vecInCone, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(xtru2, p2, vecInCone, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tess2, p2, vecInCone, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(poly2, p2, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru1, prism1_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru1, prism1_Surface, vecOnSphere, false, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess1, prism1_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess1, prism1_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(poly1, prism1_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(poly1, prism1_Surface, vecOnSphere, false, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru2, prism2_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru2, prism2_Surface, vecOnSphere, false, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess2, prism2_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess2, prism2_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(poly2, prism2_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(poly2, prism2_Surface, vecOnSphere, false, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru1, prism1_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru1, prism1_Surface, vecOnSphere, true, NLOOP);
  }
  if (iftess1) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess1, prism1_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess1, prism1_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifpoly1) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) - Convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(poly1, prism1_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(poly1, prism1_Surface, vecOnSphere, true, NLOOP);
  }
  G4cout << G4endl;
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru2, prism2_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru2, prism2_Surface, vecOnSphere, true, NLOOP);
  }
  if (iftess2) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess2, prism2_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess2, prism2_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifpoly2) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) - Non-convex prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(poly2, prism2_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(poly2, prism2_Surface, vecOnSphere, true, NLOOP);
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////
// 
//  Check Inside

void Check_Inside(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in)
{
  G4Timer time;
  clock_t t;
  G4int NP = pp.size();

  // Check points
  //
  t = clock(); 
  //time.Start();
  G4int n = 0, k = 0;
  for(int i=0; i<(NLOOP/20); ++i) {
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (k >= NP) k = 0;

    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (t1->Inside(pp[k++]) != in) { ++n; }
    if (k >= NP) k = 0;
  }
  t = clock() - t;

  G4cout <<" "<< n <<" inconsistencies";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
// 
//  Check SurfaceNormal(p)

void Check_SurfaceNormal(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP)
{
  G4Timer time;
  clock_t t;
  G4int NP = pp.size();

  // Check points
  //
  G4ThreeVector norm[20];
  t = clock(); 
  //time.Start();

  G4int k = 0;
  for(int i=0; i<(NLOOP/20); i++) {
    norm[0] = t1->SurfaceNormal(pp[k++]);
    norm[1] = t1->SurfaceNormal(pp[k++]);
    norm[2] = t1->SurfaceNormal(pp[k++]);
    norm[3] = t1->SurfaceNormal(pp[k++]);
    norm[4] = t1->SurfaceNormal(pp[k++]);
    norm[5] = t1->SurfaceNormal(pp[k++]);
    norm[6] = t1->SurfaceNormal(pp[k++]);
    norm[7] = t1->SurfaceNormal(pp[k++]);
    norm[8] = t1->SurfaceNormal(pp[k++]);
    norm[9] = t1->SurfaceNormal(pp[k++]);
    if (k >= NP) k = 0;

    norm[10] = t1->SurfaceNormal(pp[k++]);
    norm[11] = t1->SurfaceNormal(pp[k++]);
    norm[12] = t1->SurfaceNormal(pp[k++]);
    norm[13] = t1->SurfaceNormal(pp[k++]);
    norm[14] = t1->SurfaceNormal(pp[k++]);
    norm[15] = t1->SurfaceNormal(pp[k++]);
    norm[16] = t1->SurfaceNormal(pp[k++]);
    norm[17] = t1->SurfaceNormal(pp[k++]);
    norm[18] = t1->SurfaceNormal(pp[k++]);
    norm[19] = t1->SurfaceNormal(pp[k++]);
    if (k >= NP) k = 0;
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<"                  ";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
// 
//  Check DistanceToIn(p)

void Check_SafetyToIn(G4VSolid* t1, const std::vector<G4ThreeVector>& pp,G4int NLOOP)
{
  G4Timer time;
  clock_t t;
  G4int NP = pp.size();

  // Check points
  //
  t = clock(); 
  //time.Start();
  G4int n = 0, k = 0;
  for(int i=0; i<(NLOOP/20); ++i) {
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (k >= NP) k = 0;

    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToIn(pp[k++]) < 0) { ++n; }
    if (k >= NP) k = 0;
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" negative distances";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
// 
//  Check SafetyToOut(p)

void Check_SafetyToOut(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP)
{
  G4Timer time;
  clock_t t;
  G4int NP = pp.size();

  // Check points
  //
  t = clock(); 
  //time.Start();
  G4int n = 0, k = 0;
  for(int i=0; i<(NLOOP/20); i++) {
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (k >= NP) k = 0;

    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (t1->DistanceToOut(pp[k++]) < 0) { ++n; }
    if (k >= NP) k = 0;
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" negative distances";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
// 
//  Check DistanceToIn(p,v)

void Check_DistanceToIn(G4VSolid* t1,
                        const G4ThreeVector& p,
                        const std::vector<G4ThreeVector>& v, G4int NLOOP)
{
  G4Timer time;
  clock_t t;
  G4int NP = v.size();

  t = clock();
  //time.Start();
  
  G4int n = 0, k = 0;
  for(int i=0; i<(NLOOP/20); i++) {
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (k >= NP) k = 0;
      
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (t1->DistanceToIn(p,v[k++]) == kInfinity) { ++n; }
    if (k >= NP) k = 0;
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" kInfinity (no hit)";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
//
//  Check DistanceToOut(p,v)

void Check_DistanceToOut(G4VSolid* t1,
                         const std::vector<G4ThreeVector>& p,
                         const std::vector<G4ThreeVector>& vv, G4bool calcNorm, G4int NLOOP)
{
  G4Timer time;
  clock_t t;
  G4int NP = p.size();

  std::vector<G4ThreeVector> v(NP);
  for (G4int i=0; i<NP; ++i)
  {
    v[i] = vv[i];
  }

  G4bool trueNorm;
  G4ThreeVector norm;
  G4int n = 0, k = 0, ntrue = 0;

  // Check points
  //
  if (calcNorm)
  {
    t = clock();
    //time.Start();

    for(int i=0; i<(NLOOP/20); i++) {
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (k >= NP) k = 0;
      
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k; if (trueNorm) ++ntrue;
      if (k >= NP) k = 0;
    }
    t = clock() - t; //time.Stop();
  }
  else
  {
    t = clock();
    //time.Start();

    for(int i=0; i<(NLOOP/20); i++) {
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (k >= NP) k = 0;
      
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k]) == 0) { ++n; }; ++k;
      if (k >= NP) k = 0;
    }
    t = clock() - t; //time.Stop();
  }

  G4cout <<" "<< n <<" zero distances" << "; trueNorm = "<< ntrue; 
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random unit vector within cone

G4ThreeVector RandomDirection(G4double cosTheta)
{
  G4double z   = (1. - cosTheta)*G4UniformRand() + cosTheta;
  G4double rho = std::sqrt((1.+z)*(1.-z));
  G4double phi = CLHEP::twopi*G4UniformRand();
  return G4ThreeVector(rho*std::cos(phi), rho*std::sin(phi), z);
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of a right prism

G4ThreeVector PointOnPrism(G4int nsect,
                           G4double sphi, G4double dphi,
                           G4double rmax, G4double rmin, G4double dz)
{
  // Set vertices
  //
  std::vector<G4ThreeVector> baseRmaxA;
  std::vector<G4ThreeVector> baseRminA;
  std::vector<G4ThreeVector> baseRmaxB;
  std::vector<G4ThreeVector> baseRminB;

  G4double ang = dphi/nsect;
  for (G4int i=0; i<nsect+1; ++i)
  {
    G4double phi = sphi + ang*i;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    baseRmaxA.push_back(G4ThreeVector(rmax*cosphi,rmax*sinphi,-dz));
    baseRminA.push_back(G4ThreeVector(rmin*cosphi,rmin*sinphi,-dz));
    baseRmaxB.push_back(G4ThreeVector(rmax*cosphi,rmax*sinphi, dz));
    baseRminB.push_back(G4ThreeVector(rmin*cosphi,rmin*sinphi, dz));
  }

  // Find areas
  //
  G4double hmax = rmax*std::cos(ang*0.5);
  G4double hmin = rmin*std::cos(ang*0.5);
  G4double lmax = rmax*std::sin(ang*0.5);
  G4double lmin = rmin*std::sin(ang*0.5);
  G4double sbase = (lmax+lmin)*(hmax-hmin);
  G4double sext = 4*lmax*dz;
  G4double sint = 4*lmin*dz;
  G4double scut = (dphi < CLHEP::twopi) ? 2*(rmax-rmin)*dz : 0;

  // Set array of surface ares
  //
  G4double ssurf[5] = { sint*nsect, sext*nsect, sbase*nsect, sbase*nsect, scut*2 };
  for (G4int i=1; i<5; ++i) { ssurf[i] += ssurf[i-1]; }

  // Statistics, output if nsect = 0
  //
  static int counts[5] = { 0, 0, 0, 0, 0 }; // int, ext, base, base, cut
  if (counts[0] == 0 && counts[1] == 0 && counts[2] == 0 && counts[3] == 0 && counts[4] == 0)
  {
    G4cout << ((dphi == CLHEP::twopi && rmin == 0) ? "\nConvex prism " : "\nNon-convex prism ")
           << "(Nsect = " << nsect << ", Sphi = " << sphi/deg << ", Dphi = " << dphi/deg
            << ", Rmax = " << rmax  << ", Rmin = " << rmin << ", Dz = " << dz << ")" << G4endl;
    G4cout << "areas   = " << sint*nsect << ", " << sext*nsect
           << ", " << sbase*nsect << ", " << sbase*nsect << ", " << scut*2 << G4endl;
  }
  if (nsect == 0)
  {
    G4cout << "npoints = ";
    for (int i=0; i<5; ++i)
    {
      if (i > 0) G4cout << ", ";
      G4cout << counts[i];
      counts[i] = 0;
    }
    G4cout << G4endl;
    return G4ThreeVector(0,0,0);
  }

  // Select surface
  //
  G4double select = ssurf[4]*G4UniformRand();
  G4int k = 4;
  if (select <= ssurf[3]) k = 3;
  if (select <= ssurf[2]) k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;
  ++counts[k];

  // Generate point on selected surface
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // internal lateral surface
    {
      G4int i = nsect*G4UniformRand();
      G4ThreeVector v1 = baseRminA[i+1] - baseRminA[i];
      G4ThreeVector v2(0,0,2*dz);
      p = baseRminA[i] + v1*G4UniformRand() + v2*G4UniformRand();
      break;
    }
    case 1: // external lateral surface
    {
      G4int i = nsect*G4UniformRand();
      G4ThreeVector v1 = baseRmaxA[i+1] - baseRmaxA[i];
      G4ThreeVector v2(0,0,2*dz);
      p = baseRmaxA[i] + v1*G4UniformRand() + v2*G4UniformRand();
      break;
    }
    case 2: // base at -dz
    {
      G4double s1 = lmin*(hmax-hmin);
      G4double s2 = lmax*(hmax-hmin);
      G4int i = nsect*G4UniformRand();
      G4ThreeVector p0 = baseRminA[i];
      G4ThreeVector p1 = baseRminA[i+1];
      G4ThreeVector p2 = baseRmaxA[i+1];
      G4ThreeVector p3 = baseRmaxA[i];
      if ((s1+s2)*G4UniformRand() > s1) p0 = p2;
      G4double u = G4UniformRand();
      G4double v = G4UniformRand();
      if (u + v > 1.) { u = 1. - u; v = 1. - v; }
      p = (1.-u-v)*p0 + u*p1 + v*p3;
      break;
    }
    case 3: // base at +dz
    {
      G4double s1 = lmin*(hmax-hmin);
      G4double s2 = lmax*(hmax-hmin);
      G4int i = nsect*G4UniformRand();
      G4ThreeVector p0 = baseRminB[i];
      G4ThreeVector p1 = baseRminB[i+1];
      G4ThreeVector p2 = baseRmaxB[i+1];
      G4ThreeVector p3 = baseRmaxB[i];
      if ((s1+s2)*G4UniformRand() > s1) p0 = p2;
      G4double u = G4UniformRand();
      G4double v = G4UniformRand();
      if (u + v > 1.) { u = 1. - u; v = 1. - v; }
      p = (1.-u-v)*p0 + u*p1 + v*p3;
      break;
    }
    case 4: // surface at phi cut
    {
      if (G4UniformRand() < 0.5)
      {
	G4ThreeVector v1 = baseRmaxA[0] - baseRminA[0];
	G4ThreeVector v2(0,0,2*dz);
	p = baseRminA[0] + v1*G4UniformRand() + v2*G4UniformRand();
      }
      else
      {
	G4ThreeVector v1 = baseRmaxA[nsect] - baseRminA[nsect];
	G4ThreeVector v2(0,0,2*dz);
	p = baseRminA[nsect] + v1*G4UniformRand() + v2*G4UniformRand();
      }
      break;
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Create tesselated solid for regular right prism

G4TessellatedSolid* CreatePrism(G4int nsect,
                                G4double sphi, G4double dphi,
                                G4double rmax, G4double rmin, G4double dz)
{
  // Set vertices
  //
  std::vector<G4ThreeVector> baseRmaxA;
  std::vector<G4ThreeVector> baseRminA;
  std::vector<G4ThreeVector> baseRmaxB;
  std::vector<G4ThreeVector> baseRminB;

  G4double ang = dphi/nsect;
  for (G4int i=0; i<nsect+1; ++i)
  {
    G4double phi = sphi + ang*i;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    baseRmaxA.push_back(G4ThreeVector(rmax*cosphi,rmax*sinphi,-dz));
    baseRminA.push_back(G4ThreeVector(rmin*cosphi,rmin*sinphi,-dz));
    baseRmaxB.push_back(G4ThreeVector(rmax*cosphi,rmax*sinphi, dz));
    baseRminB.push_back(G4ThreeVector(rmin*cosphi,rmin*sinphi, dz));
  }

  // Set base at -dz
  //
  std::vector<G4VFacet *> faces;
  for (G4int i=0; i<nsect; ++i)
  {
    faces.push_back(new G4TriangularFacet(baseRminA[i],baseRmaxA[i+1],baseRmaxA[i],ABSOLUTE));
    if (rmin > 0)
    {
      faces.push_back(new G4TriangularFacet(baseRminA[i],baseRminA[i+1],baseRmaxA[i+1],ABSOLUTE));
    }
  }

  // Set base at +dz
  //
  for (G4int i=0; i<nsect; ++i)
  {
    faces.push_back(new G4TriangularFacet(baseRminB[i],baseRmaxB[i],baseRmaxB[i+1],ABSOLUTE));
    if (rmin > 0)
    {
      faces.push_back(new G4TriangularFacet(baseRminB[i],baseRmaxB[i+1],baseRminB[i+1],ABSOLUTE));
    }
  }

  // Set external lateral surface
  //
  for (G4int i=0; i<nsect; ++i)
  {
    faces.push_back(new G4QuadrangularFacet(baseRmaxA[i],baseRmaxA[i+1],
                                            baseRmaxB[i+1],baseRmaxB[i],ABSOLUTE));
  }

  // Set internal lateral surface, if it exists
  //
  if (rmin > 0)
  {
    for (G4int i=0; i<nsect; ++i)
    {
      faces.push_back(new G4QuadrangularFacet(baseRminA[i+1],baseRminA[i],
                                              baseRminB[i],baseRminB[i+1],ABSOLUTE));
    }
  }

  // Set two faces at phi cut, if it exists
  //
  if (dphi < CLHEP::twopi)
  {
   faces.push_back(new G4QuadrangularFacet(baseRminA[0],baseRmaxA[0],
                                           baseRmaxB[0],baseRminB[0],ABSOLUTE));
   faces.push_back(new G4QuadrangularFacet(baseRmaxA[nsect],baseRminA[nsect],
                                           baseRminB[nsect],baseRmaxB[nsect],ABSOLUTE));
  }

  // Create tesselated solid
  //
  G4TessellatedSolid* solid = new G4TessellatedSolid("TessellatedSolid");
  G4int nface = faces.size();
  for (G4int i=0; i<nface; ++i) { solid->AddFacet(faces[i]); }
  solid->SetSolidClosed(true);

  return solid;
}
