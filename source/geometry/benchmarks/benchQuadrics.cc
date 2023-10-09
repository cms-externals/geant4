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

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4EllipticalCone.hh"
#include "G4Paraboloid.hh"
#include "G4Hype.hh"
#include "G4GenericTrap.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

// ------------------------------------------------------------------------------------------------------
double MinCosine(const std::vector<G4TwoVector>& vv, double dz, int i1, int i2, int i3, int i4)
{
  G4ThreeVector v1(vv[i1].x(), vv[i1].y(),-dz);
  G4ThreeVector v2(vv[i2].x(), vv[i2].y(),-dz);
  G4ThreeVector v3(vv[i3].x(), vv[i3].y(), dz);
  G4ThreeVector v4(vv[i4].x(), vv[i4].y(), dz);

  return std::min(std::abs(((v3 - v1).unit()).z()), std::abs(((v4 - v2).unit()).z()));
}

// ------------------------------------------------------------------------------------------------------
double GetCosine(double t, double s, const std::vector<G4TwoVector>& vv, double dz, int i1, int i2, int i3, int i4)
{
  G4ThreeVector v1(vv[i1].x(), vv[i1].y(),-dz);
  G4ThreeVector v2(vv[i2].x(), vv[i2].y(),-dz);
  G4ThreeVector v3(vv[i3].x(), vv[i3].y(), dz);
  G4ThreeVector v4(vv[i4].x(), vv[i4].y(), dz);


  G4ThreeVector x1 = v1 + (v3 - v1)*t;
  G4ThreeVector x2 = v2 + (v4 - v2)*t;
  G4ThreeVector x12 = x2  - x1;

  G4ThreeVector z1 = v1 + (v2 - v1)*s;
  G4ThreeVector z2 = v3 + (v4 - v3)*s;
  G4ThreeVector z12 = z2  - z1;

  G4ThreeVector normA = G4ThreeVector(x12.y(),-x12.x(), 0).unit();
  G4ThreeVector normB = x12.cross(z12).unit();
  /*
  std::cout << "\n" << v1 << v2 << v3 << v4 << std::endl;
  std::cout << "   x12, z12 = " << x12 << z12 << std::endl;
  std::cout << "   normA,B = " << normA << normB
            << " t,s = " << t << ", " << s
            << " cosAB = " << normA.dot(normB) << std::endl;
  */
  return normA.dot(normB);
}

// ------------------------------------------------------------------------------------------------------
void CheckCosine(const std::vector<G4TwoVector>& vv, double dz, int i1, int i2, int i3, int i4)
{
  double mincos = MinCosine(vv, dz, i1, i2, i3 ,i4);
  constexpr int nparam = 100;
  for (int it = 0; it <= nparam; ++it) {
    double t = double(it)/nparam;
    for (int is = 0; is <= nparam; ++is) {
      double s = double(is)/nparam;
      double cosine = GetCosine(t, s, vv, dz, i1, i2, i3 ,i4);
      if (cosine < mincos) {
        std::cout << "=== COSINE "<< mincos << " in corners is not MIN " << cosine << " !!!"
                  << " t,s = " << t << ", " << s << std::endl;
        mincos = cosine;
      }
    }
  }
  std::cout << std::endl;
  std::cout << "=== Min COSINE = " << mincos << " !!!" << std::endl;
  std::cout << "=== Min COSINE = " << mincos << " !!!" << std::endl;
  std::cout << "=== Min COSINE = " << mincos << " !!!" << std::endl;
}

////////////////////////////////////////////////////////////////////////
//
// Declare auxiliary routines

G4double      comp_ellint_2(G4double e);
G4double      EllipticConeLateralArea(G4double pA, G4double pB, G4double pH);
G4double      RandomRadiusInRing(G4double rmin, G4double rmax);
G4TwoVector   RandomPointInRing(G4double rmin, G4double rmax);
G4TwoVector   RandomPointOnEllipse(G4double a, G4double b);
G4TwoVector   RandomPointInEllipse(G4double a, G4double b);
G4ThreeVector RandomDirection(G4double cosTheta);
G4ThreeVector RandomPointOnEllipsoid(G4double a, G4double b, G4double c);
G4double      GenericTrapFaceArea(const std::vector<G4TwoVector>& vertices,
                                  G4int i1, G4int i2, G4int i3, G4int i4, G4double z12, G4double z34);

////////////////////////////////////////////////////////////////////////
//
// Declare routines for picking random points

G4ThreeVector PointOnEllipsoid(G4double xaxis, G4double yaxis, G4double zaxis, G4double zcutmin, G4double zcutmax);
G4ThreeVector PointOnEllipticTube(G4double dx, G4double dy, G4double dz);
G4ThreeVector PointOnEllipticCone(G4double xSemiAxis, G4double ySemiAxis, G4double zheight, G4double zTopCut);
G4ThreeVector PointOnParaboloid(G4double r1, G4double r2, G4double dz);
G4ThreeVector PointOnHype(G4double innerR,  G4double outerR, G4double innerSt, G4double outerSt, G4double dz);
G4ThreeVector PointOnGenericTrap(G4double dz, const std::vector<G4TwoVector>& vertices);

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
  G4int NLOOP = 1000000;
  G4int NM = NLOOP/1000000;
  G4cout << "\n       Number of calls : "<< NM << " 000 000" << G4endl;

  // Select types of solid to test: 0 - omit, 1 - test
  G4int iforb = 1, ifellip = 1, ifetube = 1, ifecone = 1, ifpara = 1, ifhype = 1, ifgentrap = 1;

  // Generate random points
  const G4int NP = 100000;
  G4double ksur = 1e-12, kin = 0.8, kout = 1.2, del=1*cm;

  // Random points for Orb
  //
  G4double rmax;
  G4VSolid* orb = new G4Orb("Orb", rmax=8*cm);

  std::vector<G4ThreeVector> orb_Surface(NP);
  std::vector<G4ThreeVector> orb_Inside(NP);
  std::vector<G4ThreeVector> orb_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = rmax *  RandomDirection(-1);
    orb_Surface[i] = (1.0 + k*ksur)*p;
    orb_Inside [i] = kin  * p;
    orb_Outside[i] = kout * p;
  }

  // Random points for Ellipsoid
  //
  G4double ax, by, cz, zmin, zmax;
  G4VSolid* ellip = new G4Ellipsoid("Ellipsoid", ax=5*cm, by=6*cm, cz=10*cm, zmin=-5*cm, zmax=6*cm);

  std::vector<G4ThreeVector> ellip_Surface(NP);
  std::vector<G4ThreeVector> ellip_Inside(NP);
  std::vector<G4ThreeVector> ellip_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnEllipsoid(ax,by,cz,zmin,zmax);
    ellip_Surface[i] = (1.0 + k*ksur)*p;
    ellip_Inside [i] = kin  * p;
    ellip_Outside[i] = kout * p;
  }
  PointOnEllipsoid(0,0,0,0,0);

  // Random points for Elliptical tube
  //
  G4double dx, dy, dz;
  G4VSolid* etube = new G4EllipticalTube("ElTube", dx=7*cm, dy=5*cm, dz=6*cm);

  std::vector<G4ThreeVector> etube_Surface(NP);
  std::vector<G4ThreeVector> etube_Inside(NP);
  std::vector<G4ThreeVector> etube_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnEllipticTube(dx,dy,dz);
    etube_Surface[i] = (1.0 + k*ksur)*p;
    etube_Inside [i] = kin  * p;
    etube_Outside[i] = kout * p;
  }
  PointOnEllipticTube(0,0,0);

  // Random points for Elliptical cone
  //
  G4double xaxis, yaxis, zh, zcut;
  G4VSolid* econe = new G4EllipticalCone("ElCone", xaxis=1.5, yaxis=1.0, zh=10*cm, zcut=6*cm);

  std::vector<G4ThreeVector> econe_Surface(NP);
  std::vector<G4ThreeVector> econe_Inside(NP);
  std::vector<G4ThreeVector> econe_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnEllipticCone(xaxis,yaxis,zh,zcut);
    econe_Surface[i] = (1.0 + k*ksur)*p;
    econe_Inside [i] = kin  * p;
    econe_Outside[i] = kout * p;
  }
  PointOnEllipticCone(0,0,0,0);

  // Random points for Paraboloid
  //
  G4double r1, r2;
  G4VSolid* para = new G4Paraboloid("Paraboloid", dz=6*cm, r1=6*cm, r2=7*cm);

  std::vector<G4ThreeVector> para_Surface(NP);
  std::vector<G4ThreeVector> para_Inside(NP);
  std::vector<G4ThreeVector> para_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnParaboloid(r1,r2,dz);
    para_Surface[i] = (1.0 + k*ksur)*p;
    para_Inside [i] = kin  * p;
    para_Outside[i] = kout * p;
  }
  PointOnParaboloid(0,0,0);

  // Random points for Hype
  //
  G4double innerR, outerR, innerSt, outerSt;
  G4VSolid* hype = new G4Hype("Hype", innerR=2*cm, outerR=6*cm, innerSt=15*deg, outerSt=15*deg, dz=6*cm);

  std::vector<G4ThreeVector> hype_Surface(NP);
  std::vector<G4ThreeVector> hype_Inside(NP);
  std::vector<G4ThreeVector> hype_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnHype(innerR,outerR,innerSt,outerSt,dz);
    hype_Surface[i] = (1.0 + k*ksur)*p;
    hype_Inside [i] = PointOnHype(innerR+del,outerR-del,innerSt,outerSt,dz-del);
    hype_Outside[i] = PointOnHype(innerR-del,outerR+del,innerSt,outerSt,dz+del);
  }
  PointOnHype(0,0,0,0,0);

  // Random points for Generic trap
  //
  std::vector<G4TwoVector> vertices(8);
  vertices[0].set(-70*mm,-60*mm);
  vertices[1].set(-60*mm, 60*mm);
  vertices[2].set( 60*mm, 60*mm);
  vertices[3].set( 70*mm,-60*mm);
  vertices[4].set(-60*mm,  0*mm);
  vertices[5].set(  0*mm, 50*mm);
  vertices[6].set( 60*mm,  0*mm);
  vertices[7].set(  0*mm,-60*mm);
  G4VSolid* gentrap = new G4GenericTrap("Generic trap", dz = 60*mm, vertices);

  CheckCosine(vertices, dz, 0, 3, 4, 7);
  CheckCosine(vertices, dz, 1, 0, 5, 4);
  CheckCosine(vertices, dz, 2, 1, 6, 5);
  CheckCosine(vertices, dz, 3, 2, 7, 6);

  std::vector<G4ThreeVector> gentrap_Surface(NP);
  std::vector<G4ThreeVector> gentrap_Inside(NP);
  std::vector<G4ThreeVector> gentrap_Outside(NP);
  for (G4int i = 0; i < NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnGenericTrap(dz, vertices);
    gentrap_Surface[i] = (1.0 + k*ksur)*p;
    gentrap_Inside [i] = kin  * p;
    gentrap_Outside[i] = kout * p;
  }
  std::cout << "=== GetSurfaceArea : " << gentrap->GetSurfaceArea() << std::endl;
  std::cout << "=== G4VSolid::GetSurfaceArea : " << gentrap->G4VSolid::GetSurfaceArea() << std::endl;
  PointOnGenericTrap(0, vertices);

  // Check
  for (G4int i = 0; i < NP; ++i)
  {
    if (gentrap->Inside(gentrap_Surface[i]) != kSurface) std::cout << "=== Point is not on Surface !!!" << std::endl;
    if (gentrap->Inside(gentrap_Inside[i])  != kInside)  std::cout << "=== Point is not Inside !!!" << std::endl;
    if (gentrap->Inside(gentrap_Outside[i]) != kOutside) std::cout << "=== Point is not Outside !!!" << std::endl;
  }

  // Set points and generate random directions
  //
  G4ThreeVector p_orb    (5*cm,  0, 15*cm);
  G4ThreeVector p_ellip  (4.5*cm,0,  7.2*cm);
  G4ThreeVector p_etube  (5*cm,  0, 13.95*cm);
  G4ThreeVector p_econe  (5*cm,  0, 30*cm);
  G4ThreeVector p_para   (5*cm,  0, 16.6*cm);
  G4ThreeVector p_hype   (5*cm,  0, 13.85*cm);
  G4ThreeVector p_gentrap(4*cm,  0, 12.05*cm);

  G4double cosTheta=std::cos(40*deg);
  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i = 0; i < NP; ++i) {
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
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(ellip, ellip_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(ellip, ellip_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(ellip, ellip_Outside, NLOOP, kOutside);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(etube, etube_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(etube, etube_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(etube, etube_Outside, NLOOP, kOutside);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(econe, econe_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(econe, econe_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(econe, econe_Outside, NLOOP, kOutside);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(para, para_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(para, para_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(para, para_Outside, NLOOP, kOutside);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(hype, hype_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(hype, hype_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(hype, hype_Outside, NLOOP, kOutside);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(gentrap, gentrap_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(gentrap, gentrap_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(gentrap, gentrap_Outside, NLOOP, kOutside);
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(orb, orb_Surface, NLOOP, kSurface);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(ellip, ellip_Surface, NLOOP, kSurface);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(etube, etube_Surface, NLOOP, kSurface);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(econe, econe_Surface, NLOOP, kSurface);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(para, para_Surface, NLOOP, kSurface);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(hype, hype_Surface, NLOOP, kSurface);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(gentrap, gentrap_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(orb, orb_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(orb, orb_Outside, NLOOP, kOutside);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(ellip, ellip_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(ellip, ellip_Outside, NLOOP, kOutside);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(etube, etube_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(etube, etube_Outside, NLOOP, kOutside);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(econe, econe_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(econe, econe_Outside, NLOOP, kOutside);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(para, para_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(para, para_Outside, NLOOP, kOutside);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(hype, hype_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(hype, hype_Outside, NLOOP, kOutside);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(gentrap, gentrap_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(gentrap, gentrap_Outside, NLOOP, kOutside);
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(orb, orb_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(orb, orb_Surface, NLOOP, kSurface);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(ellip, ellip_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(ellip, ellip_Surface, NLOOP, kSurface);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(etube, etube_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(etube, etube_Surface, NLOOP, kSurface);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(econe, econe_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(econe, econe_Surface, NLOOP, kSurface);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(para, para_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(para, para_Surface, NLOOP, kSurface);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(hype, hype_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(hype, hype_Surface, NLOOP, kSurface);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(gentrap, gentrap_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_SafetyToOut(gentrap, gentrap_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(orb, p_orb, vecInCone, NLOOP);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(ellip, p_ellip, vecInCone, NLOOP);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(etube, p_etube, vecInCone, NLOOP);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(econe, p_econe, vecInCone, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(para, p_para, vecInCone, NLOOP);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(hype, p_hype, vecInCone, NLOOP);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(gentrap, p_gentrap, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(orb, orb_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(orb, orb_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(ellip, ellip_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(ellip, ellip_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(etube, etube_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(etube, etube_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(econe, econe_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(econe, econe_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(para, para_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(para, para_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(hype, hype_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(hype, hype_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(gentrap, gentrap_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(gentrap, gentrap_Surface, vecOnSphere, false, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (iforb) {
    G4cout << "****** Time test for G4Orb::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(orb, orb_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(orb, orb_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifellip) {
    G4cout << "****** Time test for G4Ellipsoid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(ellip, ellip_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(ellip, ellip_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifetube) {
    G4cout << "****** Time test for G4EllipticalTube::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(etube, etube_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(etube, etube_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifecone) {
    G4cout << "****** Time test for G4EllipticalCone::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(econe, econe_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(econe, econe_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifpara) {
    G4cout << "****** Time test for G4Paraboloid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(para, para_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(para, para_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifhype) {
    G4cout << "****** Time test for G4Hype::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(hype, hype_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(hype, hype_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifgentrap) {
    G4cout << "****** Time test for G4GenericTrap::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(gentrap, gentrap_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(gentrap, gentrap_Surface, vecOnSphere, true, NLOOP);
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
  for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
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
  for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
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
  for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
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
  for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
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
  for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
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
  for (G4int i = 0; i < NP; ++i)
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

    for (G4int i = 0; i < (NLOOP/20); ++i) {
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

    for (G4int i = 0; i < (NLOOP/20); ++i) {
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
  G4cout <<"   Time: "<< (double)t/CLOCKS_PER_SEC*1000 << " ms" <<  G4endl;
}

////////////////////////////////////////////////////////////////////////
//
// Compute Elliptical Integral of the Second Kind

G4double comp_ellint_2(G4double e)
{
  const G4double eps = 1. / 134217728.; // 1 / 2^27

  G4double a = 1.;
  G4double b = std::sqrt((1. - e)*(1. + e));
  if (b == 1.) return CLHEP::halfpi;
  if (b == 0.) return 1.;

  G4double x = 1.;
  G4double y = b;
  G4double S = 0.;
  G4double M = 1.;
  while (x - y > eps*y) {
    G4double tmp = (x + y) * 0.5;
    y = std::sqrt(x*y);
    x = tmp;
    M += M;
    S += M * (x - y)*(x - y);
  }
  return 0.5 * CLHEP::halfpi * ((a + b)*(a + b) - S) / (x + y);
}

////////////////////////////////////////////////////////////////////////
//
// Compute the lateral surface area of an elliptic cone

G4double EllipticConeLateralArea(G4double pA, G4double pB, G4double pH)
{
  G4double x = std::abs(pA);
  G4double y = std::abs(pB);
  G4double h = std::abs(pH);
  G4double a = std::max(x,y);
  G4double b = std::min(x,y);
  G4double e = std::sqrt((1. - b/a)*(1. + b/a)) / std::hypot(1.,b/h);
  return 2. * a * std::hypot(b,h) * comp_ellint_2(e);
}

////////////////////////////////////////////////////////////////////////
//
// Returns random radius in annular ring

G4double RandomRadiusInRing(G4double rmin, G4double rmax)
{
  if (rmin == rmax) return rmin;
  G4double k = G4UniformRand();
  return (rmin <= 0) ? rmax*std::sqrt(k)
                     : std::sqrt(k*rmax*rmax + (1.-k)*rmin*rmin);
}

////////////////////////////////////////////////////////////////////////
//
// Returns random point in annular ring

G4TwoVector RandomPointInRing(G4double rmin, G4double rmax)
{
  G4double rho = RandomRadiusInRing(rmin, rmax);
  G4double phi = CLHEP::twopi*G4UniformRand();
  return G4TwoVector(rho*std::cos(phi), rho*std::sin(phi));
}

////////////////////////////////////////////////////////////////////////
//
// Returns random point in ellipse (rejection sampling)

G4TwoVector RandomPointInEllipse(G4double a, G4double b)
{
  G4double aa = (a*a == 0) ? 0 : 1/(a*a);
  G4double bb = (b*b == 0) ? 0 : 1/(b*b);
  for (G4int i = 0; i < 1000; ++i)
  {
    G4double x = a*(2*G4UniformRand() - 1);
    G4double y = b*(2*G4UniformRand() - 1);
    if (x*x*aa + y*y*bb <= 1) return G4TwoVector(x,y);
  }
  return G4TwoVector(0,0);
}

////////////////////////////////////////////////////////////////////////
//
// Returns random point on ellipse (rejection sampling)

G4TwoVector RandomPointOnEllipse(G4double a, G4double b)
{
  G4double A = std::abs(a);
  G4double B = std::abs(b);
  G4double mu_max = std::max(A,B);

  G4double x,y;
  for (G4int i = 0; i < 1000; ++i)
  {
    G4double phi = CLHEP::twopi*G4UniformRand();
    x = std::cos(phi);
    y = std::sin(phi);
    G4double mu = std::sqrt((B*x)*(B*x) + (A*y)*(A*y));
    if (mu_max*G4UniformRand() <= mu) break;
  }
  return G4TwoVector(A*x,B*y);
}

////////////////////////////////////////////////////////////////////////
//
// Returns random point on ellipsoid (rejection sampling)

G4ThreeVector RandomPointOnEllipsoid(G4double a, G4double b, G4double c)
{
  G4double A = std::abs(a);
  G4double B = std::abs(b);
  G4double C = std::abs(c);
  G4double mu_max = std::max(std::max(A*B,A*C),B*C);

  G4ThreeVector p;
  for (G4int i = 0; i < 1000; ++i)
  {
    p = RandomDirection(-1.);
    G4double xbc = p.x()*B*C;
    G4double yac = p.y()*A*C;
    G4double zab = p.z()*A*B;
    G4double mu = std::sqrt(xbc*xbc + yac*yac + zab*zab);
    if (mu_max*G4UniformRand() <= mu) break;
  }
  return G4ThreeVector(A*p.x(),B*p.y(),C*p.z());
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
// Find area of generic trap face

G4double GenericTrapFaceArea(const std::vector<G4TwoVector>& vertices,
                             G4int i1, G4int i2, G4int i3, G4int i4, G4double z12, G4double z34)
{
  constexpr G4int NSTEP = 200;
  constexpr G4double DSTEP = 1.0/NSTEP;
  G4ThreeVector v1(vertices[i1].x(), vertices[i1].y(), z12);
  G4ThreeVector v2(vertices[i2].x(), vertices[i2].y(), z12);
  G4ThreeVector v3(vertices[i3].x(), vertices[i3].y(), z34);
  G4ThreeVector v4(vertices[i4].x(), vertices[i4].y(), z34);
  if (z12 == z34) return ((v3 - v1).cross(v4 - v2)).mag()*0.5;
  G4double area = 0.;
  for (G4int i = 0; i < NSTEP; ++i)
  {
    G4double t1 = DSTEP*i;
    G4double t2 = t1 + DSTEP;
    G4ThreeVector a1 = v1 + (v2 - v1)*t1;
    G4ThreeVector a2 = v1 + (v2 - v1)*t2;
    G4ThreeVector a3 = v4 + (v3 - v4)*t2;
    G4ThreeVector a4 = v4 + (v3 - v4)*t1;
    for (G4int k = 0; k < NSTEP; ++k)
    {
      G4double s1 = DSTEP*k;
      G4double s2 = s1 + DSTEP;
      G4ThreeVector b1 = a1 + (a4 - a1)*s1;
      G4ThreeVector b2 = a1 + (a4 - a1)*s2;
      G4ThreeVector b3 = a2 + (a3 - a2)*s2;
      G4ThreeVector b4 = a2 + (a3 - a2)*s1;
      area += ((b3 - b1).cross(b4 - b2)).mag();
    }
  }
  return area*0.5;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of ellipsoid

G4ThreeVector PointOnEllipsoid(G4double xaxis, G4double yaxis, G4double zaxis, G4double zcutmin, G4double zcutmax)
{
  // Statistics, output if xaxis = 0
  //
  static G4int counts[3] = { 0, 0, 0 };
  if (counts[0] == 0  && counts[1] == 0  && counts[2] == 0)
  {
    G4cout << "\nEllipsoid "
           << "(a=" << xaxis << ",b=" << yaxis  << ",c=" << zaxis
           << ",zmin=" << zcutmin << ",zmax=" << zcutmax << ")" << G4endl;
  }
  if (xaxis == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << G4endl;
    for (int i = 0; i < 3; ++i) { counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }

  // Generate point on unit sphere, then project to ellipsoid
  //
  G4ThreeVector p = RandomDirection(-1.);
  G4double x = p.x()*xaxis;
  G4double y = p.y()*yaxis;
  G4double z = p.z()*zaxis;
  if (z < zcutmin)
  {
    z = zcutmin;
    ++counts[0];
  }
  else if (z > zcutmax)
  {
    z = zcutmax;
    ++counts[2];
  }
  else
  {
    ++counts[1];
  }
  return G4ThreeVector(x,y,z);
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of elliptic tube

G4ThreeVector PointOnEllipticTube(G4double dx, G4double dy, G4double dz)
{
  G4double a = std::max(dx,dy);
  G4double b = std::min(dx,dy);
  G4double e = std::sqrt((1. - b/a)*(1. + b/a));
  G4double sside =  4. * a * comp_ellint_2(e);
  G4double sbase = CLHEP::pi*dx*dy;

  // Set areas (base at -Z, side surface, base at +Z)
  //
  G4double ssurf[3] = { sbase, sside, sbase };
  for (G4int i = 1; i < 3; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[2]*G4UniformRand();
  G4int k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Statistics, output if dx == 0
  //
  static G4int counts[3] = { 0, 0, 0 };

  G4double tot = 0;
  for (G4int i = 0; i < 3; ++i) { tot += counts[i]; }
  if (tot == 0)
  {
    G4cout << "\nEllipic tube "<< "(dx=" << dx << ",dy=" << dy << ",dz=" << dz << ")" << G4endl;
    G4cout << "areas   = " << sbase << ", " << sside << ", " << sbase << G4endl;
  }
  if (dx == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << G4endl;
    for (G4int i = 0; i < 3; ++i) { counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }
  ++counts[k];

  // Pick random point
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // base at -Z, uniform distribution
    {
      G4TwoVector r = RandomPointInEllipse(dx,dy);
      p.set(r.x(),r.y(),-dz);
      break;
    }
    case 1: // side surface, uniform distribution
    {
      G4double zh = dz*(2.*G4UniformRand() - 1.);
      G4TwoVector rho = RandomPointOnEllipse(dx,dy);
      p.set(rho.x(),rho.y(),zh);
      break;
    }
    case 2: // base at +Z, uniform distribution
    {
      G4TwoVector r = RandomPointInEllipse(dx,dy);
      p.set(r.x(),r.y(),dz);
      break;
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of elliptical cone

G4ThreeVector PointOnEllipticCone(G4double xSemiAxis, G4double ySemiAxis, G4double zheight, G4double zTopCut)
{
  G4double x0 = xSemiAxis*zheight; // x semi axis at z=0
  G4double y0 = ySemiAxis*zheight; // y semi axis at z=0
  G4double s0 = EllipticConeLateralArea(x0,y0,zheight);
  G4double kmin = (zTopCut >= zheight ) ? 0. : (zheight - zTopCut)/zheight;
  G4double kmax = (zTopCut >= zheight ) ? 2. : (zheight + zTopCut)/zheight;

  // Set areas (base at -Z, side surface, base at +Z)
  //
  G4double szmin =  pi*x0*y0*kmax*kmax;
  G4double szmax =  pi*x0*y0*kmin*kmin;
  G4double sside =  s0*(kmax*kmax - kmin*kmin);
  G4double ssurf[3] = { szmin, sside, szmax };
  for (G4int i = 1; i < 3; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[2]*G4UniformRand();
  G4int k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Statistics, output if zTopCut == 0
  //
  static G4int counts[3] = { 0, 0, 0 };

  G4double tot = 0;
  for (G4int i = 0; i < 3; ++i) { tot += counts[i]; }
  if (tot == 0)
  {
    G4cout << "\nEllipic cone "
           << "(a=" << xSemiAxis << ",b=" << ySemiAxis
           << ",h=" << zheight << ",zcut=" << zTopCut << ")" << G4endl;
    G4cout << "areas   = " << szmin << ", " << sside << ", " << szmax << G4endl;
  }
  if (zTopCut == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << G4endl;
    for (G4int i = 0; i < 3; ++i) { counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }
  ++counts[k];

  // Pick random point
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // base at -Z, uniform distribution, rejection sampling
    {
      G4double zh = zheight + zTopCut;
      G4TwoVector rho = RandomPointInEllipse(zh*xSemiAxis,zh*ySemiAxis);
      p.set(rho.x(),rho.y(),-zTopCut);
      break;
    }
    case 1: // side surface, uniform distribution, rejection sampling
    {
      G4double zh = RandomRadiusInRing(zheight-zTopCut, zheight+zTopCut);
      G4double a = x0;
      G4double b = y0;

      G4double hh = zheight*zheight;
      G4double aa = a*a;
      G4double bb = b*b;
      G4double R  = std::max(a,b);
      G4double mu_max = R*std::sqrt(hh + R*R);

      G4double x,y;
      for (G4int i = 0; i < 1000; ++i)
      {
        G4double phi = CLHEP::twopi*G4UniformRand();
        x = std::cos(phi);
        y = std::sin(phi);
        G4double xx = x*x;
        G4double yy = y*y;
        G4double E = hh + aa*xx + bb*yy;
        G4double F = (aa-bb)*x*y;
        G4double G = aa*yy + bb*xx;
        G4double mu = std::sqrt(E*G - F*F);
        if (mu_max*G4UniformRand() <= mu) break;
      }
      p.set(zh*xSemiAxis*x,zh*ySemiAxis*y,zheight-zh);
      break;
    }
    case 2: // base at +Z, uniform distribution, rejection sampling
    {
      G4double zh = zheight - zTopCut;
      G4TwoVector rho = RandomPointInEllipse(zh*xSemiAxis,zh*ySemiAxis);
      p.set(rho.x(),rho.y(),zTopCut);
      break;
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of paraboloid :
// r^2 = k1 * z + k2;

G4ThreeVector PointOnParaboloid(G4double r1, G4double r2, G4double dz)
{
  // Find areas
  //
  G4double rr1 = r1*r1;
  G4double rr2 = r2*r2;
  G4double h = 2*dz;

  G4double k1 = (rr2 - rr1) / h;
  G4double k2 = (rr2 + rr1) / 2;

  G4double h1 = k2/k1 - dz;
  G4double h2 = k2/k1 + dz;
  G4double hh1 = h1*h1;
  G4double hh2 = h2*h2;

  G4double A1 = rr1 + 4*hh1;
  A1 = (h1 == 0) ? 0 : pi*r1 / (6*hh1) * (std::sqrt(A1*A1*A1) - r1*r1*r1);

  G4double A2 = rr2 + 4*hh2;
  A2 = pi*r2 / (6*hh2) * (std::sqrt(A2*A2*A2) - r2*r2*r2);

  G4double szmin = pi*rr1;
  G4double sside = A2 - A1;
  G4double szmax = pi*rr2;

  // Set areas (base at -Z, lateral surface, base at +Z)
  //
  G4double ssurf[3] = { szmin, sside, szmax };
  for (G4int i = 1; i < 3; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[2]*G4UniformRand();
  G4int k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Statistics, output if dz == 0
  //
  static G4int counts[3] = { 0, 0, 0 };

  G4double tot = 0;
  for (G4int i = 0; i < 3; ++i) { tot += counts[i]; }
  if (tot == 0)
  {
    G4cout << "\nParabolid "<< "(rmin=" << r1 << ",rmax=" << r2 << ",dz=" << dz << ")" << G4endl;
    G4cout << "areas   = " << szmin << ", " << sside << ", " << szmax << G4endl;
  }
  if (dz == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << G4endl;
    for (G4int i = 0; i < 3; ++i) { counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }
  ++counts[k];

  // Pick random point
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // base at -Z, uniform distribution
    {
      G4TwoVector r = RandomPointInRing(0,r1);
      p.set(r.x(),r.y(),-dz);
      break;
    }
    case 1: // side surface, quazi-uniform distribution
    {
      G4double zh = dz*(2.*G4UniformRand() - 1.);
      G4double rho = std::sqrt(k1*zh + k2);
      G4double phi = CLHEP::twopi*G4UniformRand();
      p.set(rho*std::cos(phi),rho*std::sin(phi),zh);
      break;
    }
    case 2: // base at +Z, uniform distribution
    {
      G4TwoVector r = RandomPointInRing(0,r2);
      p.set(r.x(),r.y(),dz);
      break;
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of hype :
// r^2 - z^2 * tg^2(st) = a^2

G4ThreeVector PointOnHype(G4double innerR,  G4double outerR,
                          G4double innerSt, G4double outerSt, G4double dz)
{
  G4double r, cosSt, tanSt, alpha, t;
  G4double rrext, rrint, sext, sint, sbase;

  // Find area of outer surface
  //
  r = outerR;
  cosSt = std::cos(outerSt);
  tanSt = std::tan(outerSt);
  rrext  = dz*dz*tanSt*tanSt + r*r;

  if (outerSt == 0)
  {
    sext = CLHEP::twopi * r * (dz+dz);
  }
  else
  {
    alpha = CLHEP::twopi * r*r * cosSt/tanSt;
    t     = dz*tanSt / (r*cosSt);
    t     = std::log( t + std::sqrt(t*t + 1) );
    sext  = std::abs( alpha * (std::sinh(2*t)/2 + t) );
  }

  // Find area of inner surface
  //
  r  = innerR;
  cosSt = std::cos(innerSt);
  tanSt = std::tan(innerSt);
  rrint  = dz*dz*tanSt*tanSt + r*r;

  if (innerSt == 0)
  {
    sint = CLHEP::twopi * r * (dz+dz);
  }
  else
  {
    alpha = CLHEP::twopi * r*r * cosSt/tanSt;
    t     = dz*tanSt / (r*cosSt);
    t     = std::log( t + std::sqrt(t*t + 1) );
    sint  = std::abs( alpha * (std::sinh(2*t)/2 + t) );
  }

  // Find area of bases
  //
  sbase = pi*(rrext - rrint);

  // Set areas (inner, outer, base, base)
  //
  G4double ssurf[4] = { sint, sext, sbase, sbase };
  for (G4int i = 1; i < 4; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[3]*G4UniformRand();
  G4int k = 3;
  if (select <= ssurf[2]) k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Statistics, output if dz == 0
  //
  static G4int counts[4] = { 0, 0, 0, 0 };

  G4double tot = 0;
  for (G4int i = 0; i < 4; ++i) { tot += counts[i]; }
  if (tot == 0)
  {
    G4cout << "\nHype "<< "(innerR=" << innerR << ",outerR=" << outerR
           << ",innerSt=" << innerSt/deg << ",outerSt=" << outerSt/deg << ",dz=" << dz << ")" << G4endl;
    G4cout << "areas   = " << sint << ", " << sext << ", " << sbase << ", " << sbase << G4endl;
  }
  if (dz == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << ", " << counts[3] << G4endl;
    for (G4int i = 0; i < 4; ++i) { counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }
  ++counts[k];

  // Pick random point
  //
  G4ThreeVector p;
  switch(k)
  {
    case 0: // inner surface, quazi-uniform distribution
    {
      G4double zh  = dz*(2.*G4UniformRand() - 1.);
      G4double ts  = std::tan(innerSt);
      G4double rho = std::sqrt(zh*zh*ts*ts + innerR*innerR);
      G4double phi = CLHEP::twopi*G4UniformRand();
      p.set(rho*std::cos(phi),rho*std::sin(phi),zh);
      break;
    }
    case 1: // outer surface, quazi-uniform distribution
    {
      G4double zh  = dz*(2.*G4UniformRand() - 1.);
      G4double ts  = std::tan(outerSt);
      G4double rho = std::sqrt(zh*zh*ts*ts + outerR*outerR);
      G4double phi = CLHEP::twopi*G4UniformRand();
      p.set(rho*std::cos(phi),rho*std::sin(phi),zh);
      break;
    }
    case 2: // base at -Z, uniform distribution
    {
      G4TwoVector rho = RandomPointInRing(std::sqrt(rrint),std::sqrt(rrext));
      p.set(rho.x(),rho.y(),-dz);
      break;
    }
    case 3: // base at +Z, uniform distribution
    {
      G4TwoVector rho = RandomPointInRing(std::sqrt(rrint),std::sqrt(rrext));
      p.set(rho.x(),rho.y(),dz);
      break;
    }
  }
  return p;
}

////////////////////////////////////////////////////////////////////////
//
// Generate random point on the surface of generic trap
//

G4ThreeVector PointOnGenericTrap(G4double dz, const std::vector<G4TwoVector>& vertices)
{
  static G4double ssurf[6] = {-1.};
  static G4int counts[6] = { 0, 0, 0, 0, 0, 0 };
  G4int iface[6][4] = { {0,1,2,3}, {7,6,5,4}, {1,0,4,5}, {2,1,5,6}, {3,2,6,7}, {0,3,7,4} };

  // First call, get area of faces, print areas
  if (ssurf[0] < 0)
  {
    for (G4int i = 0; i < 6; ++i)
    {
      G4double z12 = (i == 1) ?  dz : -dz;
      G4double z34 = (i == 0) ? -dz :  dz;
      G4int i0 = iface[i][0];
      G4int i1 = iface[i][1];
      G4int i2 = iface[i][2];
      G4int i3 = iface[i][3];
      ssurf[i] = GenericTrapFaceArea(vertices, i0, i1, i2, i3, z12, z34);
    }

    G4ThreeVector bmin(DBL_MAX, DBL_MAX,-dz), bmax(-DBL_MAX,-DBL_MAX, dz);
    for (G4int i = 0; i < 8; ++i)
    {
      if (vertices[i].x() < bmin.x()) bmin.setX(vertices[i].x());
      if (vertices[i].y() < bmin.y()) bmin.setY(vertices[i].y());
      if (vertices[i].x() > bmax.x()) bmax.setX(vertices[i].x());
      if (vertices[i].y() > bmax.y()) bmax.setY(vertices[i].y());
    }

    G4cout << G4endl;
    G4cout << "GenericTrap "<< "(dz=" << dz << ",vertices)" << "   extent: "<< bmin << " " << bmax << G4endl;
    G4cout << "areas   = "
           << ssurf[0] << ", " << ssurf[1] << ", " << ssurf[2] << ", "
           << ssurf[3] << ", " << ssurf[4] << ", " << ssurf[5] << G4endl;

    for (G4int i = 1; i < 6; ++i) { ssurf[i] += ssurf[i - 1]; }
    G4cout << "=== Total area  : " << ssurf[5] << std::endl;
  }

  // Last call, print statistics
  if (dz == 0)
  {
    G4cout << "npoints = "
           << counts[0] << ", " << counts[1] << ", " << counts[2] << ", "
           << counts[3] << ", " << counts[4] << ", " << counts[5] << G4endl;
    for (G4int i = 0; i < 6; ++i) { ssurf[i] = -1.;  counts[i] = 0; }
    return G4ThreeVector(0,0,0);
  }

  // Set nodes
  G4ThreeVector pt[8];
  pt[0].set(vertices[0].x(), vertices[0].y(),-dz);
  pt[1].set(vertices[1].x(), vertices[1].y(),-dz);
  pt[2].set(vertices[2].x(), vertices[2].y(),-dz);
  pt[3].set(vertices[3].x(), vertices[3].y(),-dz);
  pt[4].set(vertices[4].x(), vertices[4].y(), dz);
  pt[5].set(vertices[5].x(), vertices[5].y(), dz);
  pt[6].set(vertices[6].x(), vertices[6].y(), dz);
  pt[7].set(vertices[7].x(), vertices[7].y(), dz);

  // Select face
  G4double select = ssurf[5]*G4UniformRand();
  G4int kface = 0;
  kface += (G4int)(select > ssurf[0]);
  kface += (G4int)(select > ssurf[1]);
  kface += (G4int)(select > ssurf[2]);
  kface += (G4int)(select > ssurf[3]);
  kface += (G4int)(select > ssurf[4]);

  // Generate random point
  G4int i0 = iface[kface][0];
  G4int i1 = iface[kface][1];
  G4int i2 = iface[kface][2];
  G4int i3 = iface[kface][3];
  if (kface < 2)
  {
    counts[kface]++;
    G4double s2 = ((pt[1] - pt[2]).cross(pt[3] - pt[2])).mag()*0.5;
    if (select > ssurf[kface] - s2) i0 = i2;
    G4double u = G4UniformRand();
    G4double v = G4UniformRand();
    if (u + v > 1.) { u = 1. - u; v = 1. - v; }
    return (1. - u - v)*pt[i0] + u*pt[i1] + v*pt[i3];
  }
  counts[kface]++;
  G4double llmax = std::max((pt[i1] - pt[i0]).mag2(), (pt[i3] - pt[i2]).mag2());
  G4ThreeVector p;
  for (G4int i = 0; i < 1000; ++i)
  {
    G4double t = G4UniformRand();
    G4ThreeVector a1 = pt[i0] +  (pt[i3] - pt[i0])*t;
    G4ThreeVector a2 = pt[i1] +  (pt[i2] - pt[i1])*t;
    G4double reject = G4UniformRand();
    if (reject*reject > (a2 - a1).mag2()/llmax) continue;
    G4double s = G4UniformRand();
    return a1 + (a2 - a1)*s;
  }
  G4cout << "=== PointOnGenericTrap: all 1000 attempts to generate random points failed !!!" << G4endl;
  return G4ThreeVector();
}
