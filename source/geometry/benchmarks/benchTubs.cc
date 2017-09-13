#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"

#include "G4Tubs.hh"
#include "G4CutTubs.hh"

#include "G4Timer.hh"
#include "Randomize.hh"

////////////////////////////////////////////////////////////////////////
//
// Generate random radius in annular ring according to uniform area

G4double RandomRadiusInRing(G4double rmin, G4double rmax)
{
  if (rmin == rmax) return rmin;
  G4double k = G4UniformRand();
  return (rmin <= 0) ? rmax*std::sqrt(k) : std::sqrt(k*rmax*rmax + (1.-k)*rmin*rmin);
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
// Generate random point on the surface of tube or tube sector 

G4ThreeVector PointOnTube(G4double fRMin, G4double fRMax, G4double fDz, G4double fSPhi, G4double fDPhi)
{
  G4double hz    = 2.*fDz;      // height
  G4double lext  = fDPhi*fRMax; // length of external circular arc
  G4double lint  = fDPhi*fRMin; // length of internal circular arc
  G4double sbase = 0.5*fDPhi*(fRMax*fRMax-fRMin*fRMin);
  G4double scut  = (fDPhi == CLHEP::twopi) ? 0. : hz*(fRMax-fRMin);

  G4double ssurf[6] = { scut, scut, sbase, sbase, hz*lext, hz*lint };
  for (G4int i=1; i<6; ++i) { ssurf[i] += ssurf[i-1]; }

  // Select surface
  //
  G4double select = ssurf[5]*G4UniformRand();
  G4int k = 5;
  if (select <= ssurf[4]) k = 4;
  if (select <= ssurf[3]) k = 3;
  if (select <= ssurf[2]) k = 2;
  if (select <= ssurf[1]) k = 1;
  if (select <= ssurf[0]) k = 0;

  // Statistics, output if Rmax == 0
  //
  static int counts[6] = {0,0,0,0,0,0};

  G4double tot = 0; for (int i=0; i<6; ++i) tot += counts[i];
  if (tot == 0)
  {
    G4cout << "\nareas   = " << scut  << ", " << scut<< ", "
	   << sbase<< ", " << sbase << ", " << hz*lext << ", " << hz*lint << G4endl;
  }
  if (fRMax == 0)
  {
    G4cout << "npoints = " << counts[0];
    for (int i=1; i<6; ++i) G4cout << ", " << counts[i];
    G4cout << G4endl;
    for (int i=0; i<6; ++i) counts[i] = 0;
    return G4ThreeVector(0,0,0);
  }
  ++counts[k];

  G4ThreeVector p;
  switch(k)
  {
    case 0: // start phi
    {
      G4double r = fRMin + (fRMax - fRMin)*G4UniformRand();
      G4double x = r*std::cos(fSPhi);
      G4double y = r*std::sin(fSPhi);
      G4double z = -fDz + 2*fDz*G4UniformRand();
      p.set(x,y,z);
      break;
    }
    case 1: // end phi
    {
      G4double r = fRMin + (fRMax - fRMin)*G4UniformRand();
      G4double x = r*std::cos(fSPhi+fDPhi);
      G4double y = r*std::sin(fSPhi+fDPhi);
      G4double z = -fDz + 2*fDz*G4UniformRand();
      p.set(x,y,z);
      break;
    }
    case 2: // base at -dz
    {
      G4double r = RandomRadiusInRing(fRMin, fRMax);
      G4double phi = fSPhi + fDPhi*G4UniformRand();
      G4double x = r*std::cos(phi);
      G4double y = r*std::sin(phi);
      G4double z = -fDz;
      p.set(x,y,z);
      break;
    }
    case 3: // base at +dz
    {
      G4double r = RandomRadiusInRing(fRMin, fRMax);
      G4double phi = fSPhi + fDPhi*G4UniformRand();
      G4double x = r*std::cos(phi);
      G4double y = r*std::sin(phi);
      G4double z = fDz;
      p.set(x,y,z);
      break;
    }
    case 4: // external surface 
    {
      G4double phi = fSPhi + fDPhi*G4UniformRand();
      G4double x = fRMax*std::cos(phi);
      G4double y = fRMax*std::sin(phi);
      G4double z = -fDz + 2*fDz*G4UniformRand();
      p.set(x,y,z);
      break;
    }
    case 5: // internal surface
    {
      G4double phi = fSPhi + fDPhi*G4UniformRand();
      G4double x = fRMin*std::cos(phi);
      G4double y = fRMin*std::sin(phi);
      G4double z = -fDz + 2*fDz*G4UniformRand();
      p.set(x,y,z);
      break;
    }
  }
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

  // Select types of tubs to test: 0 - omit, 1 - test
  G4int itype1 = 1, itype2 = 0, itype3 = 0, itype4 = 0;
  G4int itype5 = 1, itype6 = 0, itype7 = 0, itype8 = 1;

  // Generate random points
  const G4int NP = 100000;
  G4double rmin,rmax,dz,sphi,dphi, del=2*cm;

  // Type = 1 : rmin = 0, dphi = 2pi
  //
  G4VSolid* tube1 = new G4Tubs("Tube 1", rmin=0, rmax=20*cm, dz=20*cm, sphi=0*deg, dphi=360*deg);

  std::vector<G4ThreeVector> tube1_Surface(NP);
  std::vector<G4ThreeVector> tube1_Inside(NP);
  std::vector<G4ThreeVector> tube1_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube1_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube1_Inside [i] = PointOnTube(rmin,rmax-del,dz-del,sphi,dphi);
    tube1_Outside[i] = PointOnTube(rmin,rmax+del,dz+del,sphi,dphi);
  }
  PointOnTube(0,0,0,0,0);

  // Type = 2 : rmin = 0, dphi > pi
  //
  G4VSolid* tube2 = new G4Tubs("Tube 2", rmin, rmax, dz, sphi=45*deg, dphi=270*deg);

  std::vector<G4ThreeVector> tube2_Surface(NP);
  std::vector<G4ThreeVector> tube2_Inside(NP);
  std::vector<G4ThreeVector> tube2_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube2_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube2_Inside [i] = PointOnTube(rmin,rmax-del,dz-del,sphi,dphi) - G4ThreeVector(1*cm,0,0); 
    tube2_Outside[i] = PointOnTube(rmin,rmax+del,dz+del,sphi,dphi) + G4ThreeVector(1*cm,0,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Type = 3 : rmin = 0, dphi = pi
  //
  G4VSolid* tube3 = new G4Tubs("Tube 3", rmin, rmax, dz, sphi=10*deg, dphi=180*deg);

  std::vector<G4ThreeVector> tube3_Surface(NP);
  std::vector<G4ThreeVector> tube3_Inside(NP);
  std::vector<G4ThreeVector> tube3_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube3_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube3_Inside [i] = PointOnTube(rmin,rmax-del,dz-del,sphi,dphi) + G4ThreeVector(0,1*cm,0); 
    tube3_Outside[i] = PointOnTube(rmin,rmax+del,dz+del,sphi,dphi) - G4ThreeVector(0,1*cm,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Type = 4 : rmin = 0, dphi = pi
  //
  G4VSolid* tube4 = new G4Tubs("Tube 4", rmin, rmax, dz, sphi= -45*deg, dphi=90*deg);

  std::vector<G4ThreeVector> tube4_Surface(NP);
  std::vector<G4ThreeVector> tube4_Inside(NP);
  std::vector<G4ThreeVector> tube4_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube4_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube4_Inside [i] = PointOnTube(rmin,rmax-del,dz-del,sphi,dphi) + G4ThreeVector(1*cm,0,0); 
    tube4_Outside[i] = PointOnTube(rmin,rmax+del,dz+del,sphi,dphi) - G4ThreeVector(1*cm,0,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Type = 5 : rmin > 0, dphi = 2pi
  //
  G4VSolid* tube5 = new G4Tubs("Tube 5", rmin=10*cm, rmax=20*cm, dz=20*cm, sphi=0*deg, dphi=360*deg);

  std::vector<G4ThreeVector> tube5_Surface(NP);
  std::vector<G4ThreeVector> tube5_Inside(NP);
  std::vector<G4ThreeVector> tube5_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube5_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube5_Inside [i] = PointOnTube(rmin+del,rmax-del,dz-del,sphi,dphi);
    tube5_Outside[i] = PointOnTube(rmin-del,rmax+del,dz+del,sphi,dphi);
  }
  PointOnTube(0,0,0,0,0);

  // Type = 6 : rmin > 0, dphi > pi
  //
  G4VSolid* tube6 = new G4Tubs("Tube 6", rmin, rmax, dz, sphi=45*deg, dphi=270*deg);

  std::vector<G4ThreeVector> tube6_Surface(NP);
  std::vector<G4ThreeVector> tube6_Inside(NP);
  std::vector<G4ThreeVector> tube6_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube6_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube6_Inside [i] = PointOnTube(rmin+del,rmax-del,dz-del,sphi,dphi) - G4ThreeVector(1*cm,0,0); 
    tube6_Outside[i] = PointOnTube(rmin-del,rmax+del,dz+del,sphi,dphi) + G4ThreeVector(1*cm,0,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Type = 7 : rmin > 0, dphi = pi
  //
  G4VSolid* tube7 = new G4Tubs("Tube 7", rmin, rmax, dz, sphi=10*deg, dphi=180*deg);

  std::vector<G4ThreeVector> tube7_Surface(NP);
  std::vector<G4ThreeVector> tube7_Inside(NP);
  std::vector<G4ThreeVector> tube7_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube7_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube7_Inside [i] = PointOnTube(rmin+del,rmax-del,dz-del,sphi,dphi) + G4ThreeVector(0,1*cm,0); 
    tube7_Outside[i] = PointOnTube(rmin-del,rmax+del,dz+del,sphi,dphi) - G4ThreeVector(0,1*cm,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Type = 8 : rmin > 0, dphi = pi
  //
  G4VSolid* tube8 = new G4Tubs("Tube 8", rmin, rmax, dz, sphi= -45*deg, dphi=90*deg);

  std::vector<G4ThreeVector> tube8_Surface(NP);
  std::vector<G4ThreeVector> tube8_Inside(NP);
  std::vector<G4ThreeVector> tube8_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnTube(rmin,rmax,dz,sphi,dphi);
    tube8_Surface[i] = (1.0 + (1.E-12)*k)*p;
    tube8_Inside [i] = PointOnTube(rmin+del,rmax-del,dz-del,sphi,dphi) + G4ThreeVector(1*cm,0,0); 
    tube8_Outside[i] = PointOnTube(rmin-del,rmax+del,dz+del,sphi,dphi) - G4ThreeVector(1*cm,0,0);  
  }
  PointOnTube(0,0,0,0,0);

  // Generate random directions
  G4ThreeVector p1(     0,0,-57*cm);
  G4ThreeVector p2(     0,0,-51*cm);
  G4ThreeVector p3(-10*cm,0,-37*cm);
  G4ThreeVector p4( 11*cm,0,-38.3*cm);
  G4ThreeVector p5(     0,0,-55.5*cm);
  G4ThreeVector p6(     0,0,-49.8*cm);
  G4ThreeVector p7(-10*cm,2*cm,-36.4*cm);
  //G4ThreeVector p8( 13*cm,0,-35.3*cm);
  G4ThreeVector p8( 7.6*cm,0,-25*cm);

  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i=0; i<NP; i++) { 
    vecInCone[i]   = RandomDirection(std::cos(40*deg));
    vecOnSphere[i] = RandomDirection(-1);
  }

  // Mesure time of Inside(p)
  //
  G4cout << "\n   Inside(p)\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube1, tube1_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube1, tube1_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube1, tube1_Outside, NLOOP, kOutside);  
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube2, tube2_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube2, tube2_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube2, tube2_Outside, NLOOP, kOutside);  
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube3, tube3_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube3, tube3_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube3, tube3_Outside, NLOOP, kOutside);  
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube4, tube4_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube4, tube4_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube4, tube4_Outside, NLOOP, kOutside);  
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube5, tube5_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube5, tube5_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube5, tube5_Outside, NLOOP, kOutside);  
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube6, tube6_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube6, tube6_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube6, tube6_Outside, NLOOP, kOutside);  
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube7, tube7_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube7, tube7_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube7, tube7_Outside, NLOOP, kOutside);  
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::Inside(p) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tube8, tube8_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_Inside(tube8, tube8_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tube8, tube8_Outside, NLOOP, kOutside);  
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube1, tube1_Surface, NLOOP, kSurface);
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube2, tube2_Surface, NLOOP, kSurface);
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube3, tube3_Surface, NLOOP, kSurface);
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube4, tube4_Surface, NLOOP, kSurface);
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube5, tube5_Surface, NLOOP, kSurface);
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube6, tube6_Surface, NLOOP, kSurface);
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube7, tube7_Surface, NLOOP, kSurface);
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::SurfaceNormal(p) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tube8, tube8_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube1, tube1_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube1, tube1_Outside, NLOOP, kOutside);  
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube2, tube2_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube2, tube2_Outside, NLOOP, kOutside);  
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube3, tube3_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube3, tube3_Outside, NLOOP, kOutside);  
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube4, tube4_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube4, tube4_Outside, NLOOP, kOutside);  
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube5, tube5_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube5, tube5_Outside, NLOOP, kOutside);  
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube6, tube6_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube6, tube6_Outside, NLOOP, kOutside);  
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube7, tube7_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube7, tube7_Outside, NLOOP, kOutside);  
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tube8, tube8_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tube8, tube8_Outside, NLOOP, kOutside);  
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube1, tube1_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube1, tube1_Surface, NLOOP, kSurface);
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube2, tube2_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube2, tube2_Surface, NLOOP, kSurface);
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube3, tube3_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube3, tube3_Surface, NLOOP, kSurface);
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube4, tube4_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube4, tube4_Surface, NLOOP, kSurface);
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube5, tube5_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube5, tube5_Surface, NLOOP, kSurface);
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube6, tube6_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube6, tube6_Surface, NLOOP, kSurface);
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube7, tube7_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube7, tube7_Surface, NLOOP, kSurface);
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tube8, tube8_Inside,  NLOOP, kInside);  
    G4cout << "   Points on Surface "; Check_SafetyToOut(tube8, tube8_Surface, NLOOP, kSurface);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube1, p1, vecInCone, NLOOP);
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube2, p2, vecInCone, NLOOP);
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube3, p3, vecInCone, NLOOP);
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube4, p4, vecInCone, NLOOP);
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube5, p5, vecInCone, NLOOP);
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube6, p6, vecInCone, NLOOP);
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube7, p7, vecInCone, NLOOP);
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::DistanceToIn(p,v) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tube8, p8, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube1, tube1_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube1, tube1_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube2, tube2_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube2, tube2_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube3, tube3_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube3, tube3_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube4, tube4_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube4, tube4_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube5, tube5_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube5, tube5_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube6, tube6_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube6, tube6_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube7, tube7_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube7, tube7_Surface, vecOnSphere, false, NLOOP);  
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube8, tube8_Inside, vecOnSphere, false, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube8, tube8_Surface, vecOnSphere, false, NLOOP);  
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (itype1) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 1: Rmin = 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube1, tube1_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube1, tube1_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype2) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 2: Rmin = 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube2, tube2_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube2, tube2_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype3) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 3: Rmin = 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube3, tube3_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube3, tube3_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype4) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 4: Rmin = 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube4, tube4_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube4, tube4_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype5) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 5: Rmin > 0 Dphi = 2pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube5, tube5_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube5, tube5_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype6) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 6: Rmin > 0 Dphi > pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube6, tube6_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube6, tube6_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype7) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 7: Rmin > 0 Dphi = pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube7, tube7_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube7, tube7_Surface, vecOnSphere, true, NLOOP);  
  }
  if (itype8) {
    G4cout << "****** Time test for G4Tubs::DistanceToOut(p,v) Type 8: Rmin > 0 Dphi < pi ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tube8, tube8_Inside, vecOnSphere, true, NLOOP);  
    G4cout << "   Points on Surface "; Check_DistanceToOut(tube8, tube8_Surface, vecOnSphere, true, NLOOP);  
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
