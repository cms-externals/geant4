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
// Declare euxilary routines

G4ThreeVector PointOnPrism(G4int nside, G4double rho, G4double dz);
G4TessellatedSolid* CreatePrism(G4int nside, G4double rho, G4double dz);

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
  // Set number of calls
  G4int NLOOP = 10000000;
  G4int NMIL  = NLOOP/1000000;
  G4cout << "\n       Number of calls : "<< NMIL << " 000 000" << G4endl;

  // Select types of solid to test: 0 - omit, 1 - test
  G4int ifxtru1 = 1, ifxtru2 = 1, iftess = 1, ifpoly = 1;

  // Generate random points for Hexagonal prism
  const G4int NP = 100000;

  G4int nside = 6;
  G4double rho = 25*cm, dz = 20*cm, sphi = 30*deg, dphi = 360*deg;
  G4double ksur = 1e-12, kin = 0.5, kout = 2.;

  std::vector<G4ThreeVector> prism_Surface(NP);
  std::vector<G4ThreeVector> prism_Inside(NP);
  std::vector<G4ThreeVector> prism_Outside(NP);
  for (G4int i=0; i<NP; ++i)
  {
    G4int k = (i % 3) - 1; // k = -1, 0, +1
    G4ThreeVector p = PointOnPrism(nside,rho,dz);
    prism_Surface[i] = (1.0 + k*ksur)*p;
    prism_Inside [i] = kin  * p;
    prism_Outside[i] = kout * p;
  }
  PointOnPrism(0,0,0);

  // Create solids
  //
  std::vector<G4TwoVector> polygon(nside);
  G4double ang = CLHEP::twopi/nside;
  for (G4int i=0; i<nside; ++i)
  {
    G4double phi = (0.5 - i)*ang;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    polygon[i].setX(rho*cosphi);
    polygon[i].setY(rho*sinphi);
  }

  // 
  G4TwoVector off1(0,0), off2(0,0), off3(1e-10,1e-10);
  G4double scale1, scale2;
  G4VSolid* xtru1 = new G4ExtrudedSolid("Extru 1", polygon, dz, off1, scale1=1.,     off2, scale2=1.    );

  G4double del = 1e-12;
  G4VSolid* xtru2 = new G4ExtrudedSolid("Extru 2", polygon, dz, off1, scale1=1.+del, off2, scale2=1.-del);

  G4VSolid* tess = CreatePrism(nside,rho,dz);

  const G4int nrz = 4;
  G4double zz[nrz] = {-dz,-dz, dz, dz}, rr[nrz] = { 0, rho, rho, 0 };
  G4VSolid* phedr = new G4Polyhedra("Polyhedra 1", sphi, dphi, nside, nrz, rr, zz);

  // Set points and generate random directions
  //
  G4ThreeVector p0(0,0,-62.2*cm);
  G4double cosTheta=std::cos(40*deg);

  std::vector<G4ThreeVector> vecInCone(NP);
  std::vector<G4ThreeVector> vecOnSphere(NP);
  for (int i=0; i<NP; i++) {
    vecInCone[i]   = RandomDirection(cosTheta);
    vecOnSphere[i] = RandomDirection(-1);
  }

  // Mesure time of Inside(p)
  //
  G4cout << "\n   Inside(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::Inside(p) - Right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(xtru1, prism_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(xtru1, prism_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(xtru1, prism_Outside, NLOOP, kOutside);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::Inside(p) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(xtru2, prism_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(xtru2, prism_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(xtru2, prism_Outside, NLOOP, kOutside);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(tess, prism_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(tess, prism_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(tess, prism_Outside, NLOOP, kOutside);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::Inside(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_Inside(phedr, prism_Inside,  NLOOP, kInside);
    G4cout << "   Points on Surface "; Check_Inside(phedr, prism_Surface, NLOOP, kSurface);
    G4cout << "   Points Outside    "; Check_Inside(phedr, prism_Outside, NLOOP, kOutside);
  }

  // Mesure time of SurfaceNormal(p)
  //
  G4cout << "\n   SurfaceNormal(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::SurfaceNormal(p) - Right prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(xtru1, prism_Surface, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::SurfaceNormal(p) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(xtru2, prism_Surface, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(tess, prism_Surface, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::SurfaceNormal(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SurfaceNormal(phedr, prism_Surface, NLOOP);
  }
 
  // Mesure time of DistanceToIn(p)
  //
  G4cout << "\n   SafetyToIn(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p) - Right prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(xtru1, prism_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(xtru1, prism_Outside, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(xtru2, prism_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(xtru2, prism_Outside, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(tess, prism_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(tess, prism_Outside, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p) ******" << G4endl;
    G4cout << "   Points on Surface "; Check_SafetyToIn(phedr, prism_Surface, NLOOP);
    G4cout << "   Points Outside    "; Check_SafetyToIn(phedr, prism_Outside, NLOOP);
  }

  // Mesure time of DistanceToOut(p)
  //
  G4cout << "\n   SafetyToOut(p)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p) - Right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(xtru1, prism_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(xtru1, prism_Surface, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(xtru2, prism_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(xtru2, prism_Surface, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(tess, prism_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(tess, prism_Surface, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_SafetyToOut(phedr, prism_Inside,  NLOOP); 
    G4cout << "   Points on Surface "; Check_SafetyToOut(phedr, prism_Surface, NLOOP);
  }

  // Mesure time of DistanceToIn(p,v)
  //
  G4cout << "\n   DistanceToIn(p,v)\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p,v) - Right prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(xtru1, p0, vecInCone, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToIn(p,v) - Quasi-right prism ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(xtru2, p0, vecInCone, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(tess, p0, vecInCone, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToIn(p,v) ******" << G4endl;
    G4cout << "   DistanceToIn      "; Check_DistanceToIn(phedr, p0, vecInCone, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) without calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) without calculation of Normal\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru1, prism_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru1, prism_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru2, prism_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru2, prism_Surface, vecOnSphere, false, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess, prism_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess, prism_Surface, vecOnSphere, false, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(phedr, prism_Inside, vecOnSphere, false, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(phedr, prism_Surface, vecOnSphere, false, NLOOP);
  }

  // Mesure time of DistanceToOut(p,v) with calculation of Normal
  //
  G4cout << "\n   DistanceToOut(p,v) with calculation of Normal\n" << G4endl;
  if (ifxtru1) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru1, prism_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru1, prism_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifxtru2) {
    G4cout << "****** Time test for G4ExtrudedSolid::DistanceToOut(p,v) - Quasi-right prism ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(xtru2, prism_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(xtru2, prism_Surface, vecOnSphere, true, NLOOP);
  }
  if (iftess) {
    G4cout << "****** Time test for G4TessellatedSolid::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(tess, prism_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(tess, prism_Surface, vecOnSphere, true, NLOOP);
  }
  if (ifpoly) {
    G4cout << "****** Time test for G4Polyhedra::DistanceToOut(p,v) ******" << G4endl;
    G4cout << "   Points Inside     "; Check_DistanceToOut(phedr, prism_Inside, vecOnSphere, true, NLOOP);
    G4cout << "   Points on Surface "; Check_DistanceToOut(phedr, prism_Surface, vecOnSphere, true, NLOOP);
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

G4ThreeVector PointOnPrism(G4int nside, G4double rho, G4double dz)
{
  G4double ang = CLHEP::twopi/nside;
  G4double dx = rho*std::cos(0.5*ang);
  G4double dy = rho*std::sin(0.5*ang);

  // Compute areas
  //
  G4double sbase = dx*dy;
  G4double slat  = 4*dy*dz;

  // Statistics, output if nside = 0
  //
  static int counts[3] = { 0, 0, 0 };
  if (counts[0] == 0  && counts[1] == 0  && counts[2] == 0)
  {
    G4cout << "\nRight prism "
           << "(nside=" << nside << ",rho=" << rho  << ",dz=" << dz << ")" << G4endl;
    G4cout << "areas   = " << sbase*nside << ", " << sbase*nside << ", " << slat*nside << G4endl;
  }
  if (nside == 0)
  {
    G4cout << "npoints = " << counts[0] << ", " << counts[1] << ", " << counts[2] << G4endl;
    for (int i=0; i<3; ++i) counts[i] = 0;
    return G4ThreeVector(0,0,0);
  }

  // Set vertices
  //
  G4ThreeVector a1(dx,-dy,-dz), b1(dx,dy,-dz);
  G4ThreeVector a2(dx,-dy, dz), b2(dx,dy, dz);

  // Generate point
  //
  G4double ka = G4UniformRand();
  G4double kb = G4UniformRand();

  G4ThreeVector p;
  G4double select = (sbase + sbase + slat)*G4UniformRand();
  if (select <= sbase)
  {
    p = a1*ka + b1*kb;
    if (ka + kb > 1.) p = a1 + b1 - p;
    p.setZ(-dz);
    ++counts[0];
  }
  else if (select <= sbase + sbase)
  {
    p = a2*ka + b2*kb;
    if (ka + kb > 1.) p = a2 + b2 - p;
    p.setZ(dz);
    ++counts[1];
  }
  else
  {   
    p = a1 + (a2-a1)*ka + (b1-a1)*kb; 
    ++counts[2];
  }

  // Rotate point
  G4int n = nside*G4UniformRand();
  G4double phi = n*ang;
  G4double cosphi = std::cos(phi); 
  G4double sinphi = std::sin(phi);
  G4double x = p.x()*cosphi - p.y()*sinphi;  
  G4double y = p.x()*sinphi + p.y()*cosphi;  
  G4double z = p.z();
  return G4ThreeVector(x,y,z);
}

////////////////////////////////////////////////////////////////////////
//
// Create tesselated solid for regular right prism

G4TessellatedSolid* CreatePrism(G4int nside, G4double rho, G4double dz)
{
  std::vector<G4ThreeVector> baseA;
  std::vector<G4ThreeVector> baseB;

  G4double ang = CLHEP::twopi/nside;
  for (G4int i=0; i<nside+1; ++i)
  {
    G4double phi = (0.5 - i)*ang;
    G4double cosphi = std::cos(phi);
    G4double sinphi = std::sin(phi);
    baseA.push_back(G4ThreeVector(rho*cosphi,rho*sinphi,-dz));
    baseB.push_back(G4ThreeVector(rho*cosphi,rho*sinphi, dz));
  }

  // Set base at -dz
  //
  std::vector<G4VFacet *> faces;
  G4ThreeVector A0(0,0,-dz);
  for (G4int i=0; i<nside; ++i)
  { 
    faces.push_back(new G4TriangularFacet(A0,baseA[i],baseA[i+1],ABSOLUTE)); 
  }

  // Set base at +dz
  //
  G4ThreeVector B0(0,0,dz);
  for (G4int i=0; i<nside; ++i)
  { 
    faces.push_back(new G4TriangularFacet(B0,baseB[i+1],baseB[i],ABSOLUTE)); 
  }

  // Set lateral surface
  //
  for (G4int i=0; i<nside; ++i)
  { 
    faces.push_back(new G4QuadrangularFacet(baseA[i],baseB[i],baseB[i+1],baseA[i+1],ABSOLUTE)); 
  }

  G4TessellatedSolid* solid = new G4TessellatedSolid("TessellatedSolid");
  for (G4int i=0; i<3*nside; ++i)
  { 
    solid->AddFacet(faces[i]);
  }
  solid->SetSolidClosed(true);  

  return solid;
}
