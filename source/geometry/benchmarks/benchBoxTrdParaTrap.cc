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

#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Para.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

G4ThreeVector RandomDirection(G4double cosTheta)
{
  G4double z   = (1. - cosTheta)*G4UniformRand() + cosTheta;
  G4double rho = std::sqrt((1.+z)*(1.-z));
  G4double phi = CLHEP::twopi*G4UniformRand();
  return G4ThreeVector(rho*std::cos(phi), rho*std::sin(phi), z);
}

G4ThreeVector PointOnBox(G4double dx, G4double dy, G4double dz)
{
  G4double sxy = dx*dy, sxz = dx*dz, syz = dy*dz;
  G4double select = (sxy + sxz + syz)*G4UniformRand();

  if (select < sxy)
  {
    return G4ThreeVector((2*G4UniformRand() - 1)*dx,
                         (2*G4UniformRand() - 1)*dy,
                         (select < 0.5*sxy) ? -dz : dz);
  }
  if (select < sxy + sxz)
  {
    return G4ThreeVector((2*G4UniformRand() - 1)*dx,
                         (select < sxy + 0.5*sxz) ? -dy : dy,
                         (2*G4UniformRand() - 1)*dz);
  }
  else
  {
    return G4ThreeVector((select < sxy + sxz + 0.5*syz) ? -dx : dx,
                         (2*G4UniformRand() - 1)*dy,
                         (2*G4UniformRand() - 1)*dz);
  }
}

void Check_Inside       (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SurfaceNormal(G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SafetyToIn   (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_SafetyToOut  (G4VSolid* t1, const std::vector<G4ThreeVector>& pp, G4int NLOOP, EInside in);
void Check_DistanceToIn (G4VSolid* t1, const G4ThreeVector& p, const std::vector<G4ThreeVector>& vv, G4int NLOOP);
void Check_DistanceToOut(G4VSolid* t1,
                         const std::vector<G4ThreeVector>& pp,
                         const std::vector<G4ThreeVector>& vv, G4bool calcNorm, G4int NLOOP);

int main() 
{
  // Set number of calls
  G4int NLOOP = 100000000;
  G4int NM = NLOOP/1000000; 
  G4cout << "\n       Number of calls : "<< NM << " 000 000" << G4endl;

  // Select shape to test: 0 - omit, 1 - test
  G4int ifbox = 1, iftrd = 1, ifpara = 1, iftrap = 1;

  // Generate random points
  const G4int NP = 100000;
  G4double dx = 20*cm, dy = 20*cm, dz = 20*cm;

  std::vector<G4ThreeVector>  pointsSurface(NP);
  std::vector<G4ThreeVector>  pointsInside(NP);
  std::vector<G4ThreeVector>  pointsOutside(NP);

  for (int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1 
    G4ThreeVector p = PointOnBox(dx,dy,dz);
    pointsSurface[i] = (1.0 + (1.E-12)*k)*p;
    pointsInside[i]  = 0.5*p;
    pointsOutside[i] = 2.0*p;
  }

  // Generate random directions
  G4ThreeVector p0(0,0,-60*cm);
  G4double cosTheta=std::cos(40*deg);
  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i=0; i<NP; i++) { 
    vecInCone[i] = RandomDirection(cosTheta);
    vecOnSphere[i] = RandomDirection(-1);
  }

  // -----------------------------------------------------------------------------

  G4double Dz,Theta,Phi, Dy1,Dx1,Dx2,Alpha1, Dy2,Dx3,Dx4,Alpha2;
  G4double del    = 1.E-12; // vertex disturbance
  G4double angdel = 1.E-12; // angle disturbance

  // Define box
  G4VSolid* box = new G4Box("Box", dx, dy, dz);

  // Define G4Trd
  G4VSolid* trd = new G4Trd("Trd",
    Dx1 = dx, Dx2 = dx-del,
    Dy1 = dy, Dy2 = dy,
    Dz = dz);

  // Define G4Para
  G4VSolid* para = new G4Para("Para", dx, dy, dz, Alpha1=0*deg, Theta = angdel*deg, Phi = 45*deg);

  // Define G4Trap
  G4VSolid* trap = new G4Trap("Trap",
    Dz  = dz, Theta = 0, Phi = 0,
    Dy1 = dy, Dx1 = dx, Dx2 = dx, Alpha1 = 0,
    Dy2 = dy-del, Dx3 = dx-del, Dx4 = dx-del, Alpha2 = 0
  );

  G4VSolid* trap1 = new G4Trap("Trap",
    Dz  = dz, Theta = angdel*deg, Phi = 0,
    Dy1 = dy, Dx1 = dx, Dx2 = dx, Alpha1 = 0,
    Dy2 = dy, Dx3 = dx-del, Dx4 = dx-del, Alpha2 = 0
  );

  G4VSolid* trap2 = new G4Trap("Trap",
    Dz  = dz, Theta = 0, Phi = 0,
    Dy1 = dy, Dx1 = dx, Dx2 = dx, Alpha1 = 0,
    Dy2 = dy, Dx3 = dx-del, Dx4 = dx-del, Alpha2 = 0
  );

  // -----------------------------------------------------------------------------

  // Mesure time of Inside(p)
  //
  G4cout << "\n   Inside(p)\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(box, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(box, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(box, pointsOutside, NLOOP, kOutside);  
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(trd, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(trd, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(trd, pointsOutside, NLOOP, kOutside);  
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(para, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(para, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(para, pointsOutside, NLOOP, kOutside);  
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(trap, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(trap, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(trap, pointsOutside, NLOOP, kOutside);  

    G4cout << "****** Time test for G4Trap::Inside(p) Type 1 ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(trap1, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(trap1, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(trap1, pointsOutside, NLOOP, kOutside);  

    G4cout << "****** Time test for G4Trap::Inside(p) Type 2 ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(trap2, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(trap2, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(trap2, pointsOutside, NLOOP, kOutside);  
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(box, pointsSurface, NLOOP, kSurface);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(trd, pointsSurface, NLOOP, kSurface);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(para, pointsSurface, NLOOP, kSurface);
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(trap, pointsSurface, NLOOP, kSurface);

    G4cout << "****** Time test for G4Trap::SurfaceNormal(p) Type 1 ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(trap1, pointsSurface, NLOOP, kSurface);

    G4cout << "****** Time test for G4Trap::SurfaceNormal(p) Type 2 ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(trap2, pointsSurface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(box, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(box, pointsOutside, NLOOP, kOutside);  
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(trd, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(trd, pointsOutside, NLOOP, kOutside);  
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(para, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(para, pointsOutside, NLOOP, kOutside);  
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(trap, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(trap, pointsOutside, NLOOP, kOutside);  

    G4cout << "****** Time test for G4Trap::DistanceToIn(p) Type 1 ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(trap1, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(trap1, pointsOutside, NLOOP, kOutside);  

    G4cout << "****** Time test for G4Trap::DistanceToIn(p) Type 2 ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(trap2, pointsSurface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(trap2, pointsOutside, NLOOP, kOutside);  
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(box, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(box, pointsSurface, NLOOP, kSurface);
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(trd, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(trd, pointsSurface, NLOOP, kSurface);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(para, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(para, pointsSurface, NLOOP, kSurface);
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(trap, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(trap, pointsSurface, NLOOP, kSurface);

    G4cout << "****** Time test for G4Trap::DistanceToOut(p) Type 1 ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(trap1, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(trap1, pointsSurface, NLOOP, kSurface);

    G4cout << "****** Time test for G4Trap::DistanceToOut(p) Type 2 ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(trap2, pointsInside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(trap2, pointsSurface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::DistanceToIn(p,v) ******" <<G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(box, p0, vecInCone, NLOOP);  
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToIn(p,v) ******" <<G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(trd, p0, vecInCone, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::DistanceToIn(p,v) ******" <<G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(para, p0, vecInCone, NLOOP);
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::DistanceToIn(p,v) ******" <<G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(trap, p0, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(box, pointsInside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(box, pointsSurface, vecOnSphere, false, NLOOP);  
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trd, pointsInside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trd, pointsSurface, vecOnSphere, false, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(para, pointsInside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(para, pointsSurface, vecOnSphere, false, NLOOP);
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trap, pointsInside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trap, pointsSurface, vecOnSphere, false, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (ifbox) {
    G4cout << "****** Time test for G4Box::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(box, pointsInside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(box, pointsSurface, vecOnSphere, true, NLOOP);  
  }
  if (iftrd) {
    G4cout << "****** Time test for G4Trd::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trd, pointsInside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trd, pointsSurface, vecOnSphere, true, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Para::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(para, pointsInside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(para, pointsSurface, vecOnSphere, true, NLOOP);
  }
  if (iftrap) {
    G4cout << "****** Time test for G4Trap::DistanceToOut(p,v) ******" <<G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(trap, pointsInside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(trap, pointsSurface, vecOnSphere, true, NLOOP);
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
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (k >= NP) k = 0;
      
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
      if (t1->DistanceToOut(p[k],v[k],true,&trueNorm,&norm) == 0) { ++n; }; ++k;
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

  G4cout <<" "<< n <<" zero distances"; 
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC << " sec" <<  G4endl;
}
