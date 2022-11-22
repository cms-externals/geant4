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

#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"
#include "G4GeomTools.hh"

#include "G4Orb.hh"
#include "G4Trd.hh"
#include "G4Tet.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

////////////////////////////////////////////////////////////////////////
//
// Generate random direction within cone

G4ThreeVector RandomDirection(G4double cosTheta)
{
  G4double z   = (1. - cosTheta)*G4UniformRand() + cosTheta;
  G4double rho = std::sqrt((1.+z)*(1.-z));
  G4double phi = CLHEP::twopi*G4UniformRand();
  return G4ThreeVector(rho*std::cos(phi), rho*std::sin(phi), z);
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point in triangle

G4ThreeVector PointInTriangle(const G4ThreeVector& p1,
                              const G4ThreeVector& p2,
                              const G4ThreeVector& p3)
{
  G4double r1 = G4UniformRand();
  G4double r2 = G4UniformRand();
  return (r1 + r2 < 1.) ?
    p1 + (p2 - p1)*r1 + (p3 - p1)*r2 :
    p1 + (p2 - p1)*(1. - r1) + (p3 - p1)*(1. - r2);
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of a tetrahedron

G4ThreeVector PointOnTet(const G4ThreeVector& p1,
                         const G4ThreeVector& p2,
                         const G4ThreeVector& p3,
                         const G4ThreeVector& p4)
{
  G4double s1 = G4GeomTools::TriangleAreaNormal(p1, p2, p3).mag();
  G4double s2 = G4GeomTools::TriangleAreaNormal(p1, p3, p4).mag();
  G4double s3 = G4GeomTools::TriangleAreaNormal(p1, p4, p2).mag();
  G4double s4 = G4GeomTools::TriangleAreaNormal(p2, p3, p4).mag();
  G4double select = (s1 + s2 + s3 + s4)*G4UniformRand();

  G4int k = 0;
  if (select > s1) k = 1;
  if (select > s1 + s2) k = 2;
  if (select > s1 + s2 + s3) k = 3;

  G4ThreeVector p(0., 0., 0.);
  if (k == 0) p = PointInTriangle(p1, p2, p3);
  if (k == 1) p = PointInTriangle(p1, p3, p4);
  if (k == 2) p = PointInTriangle(p1, p4, p2);
  if (k == 3) p = PointInTriangle(p2, p3, p4);
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Declare bench routines

void Check_Inside       (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SurfaceNormal(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SafetyToIn   (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SafetyToOut  (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_DistanceToIn (G4VSolid* t1, const G4ThreeVector& p, const std::vector<G4ThreeVector>& vv, G4int NLOOP);
void Check_DistanceToOut(G4VSolid* t1,
                         const std::vector<G4ThreeVector>& pp,
                         const std::vector<G4ThreeVector>& vv, G4bool calcNorm, G4int NLOOP);

////////////////////////////////////////////////////////////////////////
//
// Main

int main()
{
  // Set number of calls
  G4int NLOOP = 100000000;
  G4int NM = NLOOP/1000000;
  G4cout << "\n       Number of calls : "<< NM << " 000 000" << G4endl;

  // Select shape to test: 0 - omit, 1 - test
  G4int iforb = 1, iftrd = 1, iftet = 1;

  // Generate random points
  const G4int NP = 100000;
  G4double ksur = 1e-12, kin = 0.8, kout = 1.2;

  // Random points for Orb
  //
  G4double rmax = 8*cm;
  G4VSolid* orb = new G4Orb("Orb", rmax);

  std::vector<G4ThreeVector> orb_Surface(NP);
  std::vector<G4ThreeVector> orb_Inside(NP);
  std::vector<G4ThreeVector> orb_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = rmax *  RandomDirection(-1);
    orb_Surface[i] = (1.0 + k*ksur)*p;
    orb_Inside [i] = kin  * p;
    orb_Outside[i] = kout * p;
  }

  // Create Trd
  G4double dx = 6*cm;
  G4double dy = 6*cm;
  G4double dz = 6*cm;
  G4VSolid* trd = new G4Trd("Trd", dx, 0., 0., dy, dz);

  // Random points for Tet
  //
  G4ThreeVector p1(-dx,   0.,-dz);
  G4ThreeVector p2( dx,   0.,-dz);
  G4ThreeVector p3(  0.,-dy,  dz);
  G4ThreeVector p4(  0., dy,  dz);
  G4VSolid* tet = new G4Tet("Tet", p1, p2, p3, p4);

  std::vector<G4ThreeVector> tet_Surface(NP);
  std::vector<G4ThreeVector> tet_Inside(NP);
  std::vector<G4ThreeVector> tet_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTet(p1, p2, p3, p4);
    tet_Surface[i] = (1.0 + k*ksur)*p;
    tet_Inside [i] = kin  * p;
    tet_Outside[i] = kout * p;
  }

  // Set points and generate random directions
  //
  G4ThreeVector p_orb(5*cm, 0.,  15*cm);
  G4ThreeVector p_tet(0.,   0., 9.5*cm);

  G4double cosTheta=std::cos(40*deg);
  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i=0; i<NP; i++) {
    vecInCone[i] = -RandomDirection(cosTheta);
    vecOnSphere[i] = RandomDirection(-1);
  }

  // Mesure time of Inside(p)
  //
  G4cout << "\n   Inside(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(orb, orb_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(orb, orb_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(orb, orb_Outside, NLOOP, kOutside);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(trd, tet_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(trd, tet_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(trd, tet_Outside, NLOOP, kOutside);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tet, tet_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(tet, tet_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tet, tet_Outside, NLOOP, kOutside);
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(orb, orb_Surface, NLOOP, kSurface);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(trd, tet_Surface, NLOOP, kSurface);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tet, tet_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(orb, orb_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(orb, orb_Outside, NLOOP, kOutside);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(trd, tet_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(trd, tet_Outside, NLOOP, kOutside);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tet, tet_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tet, tet_Outside, NLOOP, kOutside);
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(orb, orb_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(orb, orb_Surface, NLOOP, kSurface);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(trd, tet_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(trd, tet_Surface, NLOOP, kSurface);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tet, tet_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(tet, tet_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(orb, p_orb, vecInCone, NLOOP);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(trd, p_tet, vecInCone, NLOOP);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tet, p_tet, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(orb, orb_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(orb, orb_Surface, vecOnSphere, false, NLOOP);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trd, tet_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trd, tet_Surface, vecOnSphere, false, NLOOP);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tet, tet_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tet, tet_Surface, vecOnSphere, false, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(orb, orb_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(orb, orb_Surface, vecOnSphere, true, NLOOP);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trd, tet_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trd, tet_Surface, vecOnSphere, true, NLOOP);
  }
  if (iftet) {
    G4cout << "****** Time test for G4Tet::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tet, tet_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tet, tet_Surface, vecOnSphere, true, NLOOP);
  }

  return 0;
}

//////////////////////////////////////////////////////////////
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
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    k *= (k < NP);

    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    n += (t1->Inside(pp[k++]) != in);
    k *= (k < NP);
  }
  t = clock() - t;

  G4cout <<" "<< n <<" inconsistencies";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

//////////////////////////////////////////////////////////////
//
//  Check SurfaceNormal(p)

void Check_SurfaceNormal(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in)
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
    k *= (k < NP);

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
    k *= (k < NP);
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<"                  ";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

//////////////////////////////////////////////////////////////
//
//  Check DistanceToIn(p)

void Check_SafetyToIn(G4VSolid* t1, const std::vector<G4ThreeVector>& pp,G4int NLOOP, EInside in)
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
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    k *= (k < NP);

    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    n += (t1->DistanceToIn(pp[k++]) < 0);
    k *= (k < NP);
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" negative distances";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

//////////////////////////////////////////////////////////////
//
//  Check SafetyToOut(p)

void Check_SafetyToOut(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in)
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
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    k *= (k < NP);

    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    n += (t1->DistanceToOut(pp[k++]) < 0);
    k *= (k < NP);
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" negative distances";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

//////////////////////////////////////////////////////////////
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
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    k *= (k < NP);

    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    n += (t1->DistanceToIn(p,v[k++]) == kInfinity);
    k *= (k < NP);
  }
  t = clock() - t;
  //time.Stop();

  G4cout <<" "<< n <<" kInfinity (no hit)";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}

//////////////////////////////////////////////////////////////
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
  G4int n = 0, k = 0;

  // Check points
  //
  if (calcNorm)
  {
    t = clock();
    //time.Start();

    for(int i=0; i<(NLOOP/20); i++) {
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      k *= (k < NP);

      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0); ++k;
      k *= (k < NP);
    }
    t = clock() - t; //time.Stop();
  }
  else
  {
    t = clock();
    //time.Start();

    for(int i=0; i<(NLOOP/20); i++) {
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      k *= (k < NP);

      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      n += (t1->DistanceToOut(p[k],v[k]) == 0); ++k;
      k *= (k < NP);
    }
    t = clock() - t; //time.Stop();
  }

  G4cout <<" "<< n <<" zero distances";
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}
