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
// Unit test for G4Box
//
// Author: E.Tcherniaev (evgueni.tcherniaev@cern.ch)
//

#undef NDEBUG  // ensure asserts are compiled in
#include "G4Box.hh"
#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4Para.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "globals.hh"

#include "geomdefs.hh"
#include <assert.h>

void CheckScratching(G4int k, G4double dist)
{
  if (dist != kInfinity)
  {
    G4cout << " Test for scratching, case #" << k << " failed,   dist = " << dist
           << ", expected kInfinity" << G4endl;
  }
}

G4bool ApproxEqual(G4double a, const G4double b)
{
  G4double tolerance = (std::abs(a) + std::abs(b)) * 0.5e-14;
  G4double diff = std::abs(a - b);
  if (diff > tolerance)
  {
    G4cout << "\na = " << a << " b = " << b << " abs(a-b) = " << diff
           << " tolerance = " << tolerance << G4endl;
  }
  return diff <= tolerance;
}

G4bool ApproxEqual(const G4ThreeVector& a, const G4ThreeVector& b, G4double tolerance)
{
  G4double diff = (a - b).mag();
  if (diff > tolerance)
  {
    G4cout << "\na" << a << " b" << b << " mag(a-b) = " << diff << " tolerance = " << tolerance
           << G4endl;
  }
  return diff <= tolerance;
}

void getVertices(const G4VSolid* solid, std::vector<G4ThreeVector>& vv)
{
  G4double dx = 0, dy = 0, dz = 0;
  if (solid->GetEntityType() == "G4Box")
  {
    dx = ((const G4Box*)solid)->GetXHalfLength();
    dy = ((const G4Box*)solid)->GetYHalfLength();
    dz = ((const G4Box*)solid)->GetZHalfLength();
  }
  else if (solid->GetEntityType() == "G4Para")
  {
    dx = ((const G4Para*)solid)->GetXHalfLength();
    dy = ((const G4Para*)solid)->GetYHalfLength();
    dz = ((const G4Para*)solid)->GetZHalfLength();
  }
  else if (solid->GetEntityType() == "G4Trd")
  {
    dx = ((const G4Trd*)solid)->GetXHalfLength1();
    dy = ((const G4Trd*)solid)->GetYHalfLength1();
    dz = ((const G4Trd*)solid)->GetZHalfLength();
  }
  else if (solid->GetEntityType() == "G4Trap")
  {
    dx = ((const G4Trap*)solid)->GetXHalfLength1();
    dy = ((const G4Trap*)solid)->GetYHalfLength1();
    dz = ((const G4Trap*)solid)->GetZHalfLength();
  }
  vv[0].set(-dx, -dy, -dz);
  vv[1].set(-dx, dy, -dz);
  vv[2].set(dx, dy, -dz);
  vv[3].set(dx, -dy, -dz);
  vv[4].set(-dx, -dy, dz);
  vv[5].set(-dx, dy, dz);
  vv[6].set(dx, dy, dz);
  vv[7].set(dx, -dy, dz);
}

G4ThreeVector getPoint(const G4ThreeVector& p1,  // p2    p3
                       const G4ThreeVector& p2,  //  +----+
                       const G4ThreeVector& p3,  //  |    |
                       const G4ThreeVector& p4,  //  +----+
                       G4double s, G4double t)  // p1    p4
{
  G4ThreeVector a = p1 + (p2 - p1) * s;
  G4ThreeVector b = p4 + (p3 - p4) * s;
  return a + (b - a) * t;
}

G4ThreeVector getNormal(const G4ThreeVector& p1,  // p2    p3
                        const G4ThreeVector& p2,  //  +----+
                        const G4ThreeVector& p3,  //  |    |
                        const G4ThreeVector& p4,  //  +----+
                        G4double s, G4double t)  // p1    p4
{
  G4ThreeVector a = p1 + (p2 - p1) * s;
  G4ThreeVector b = p4 + (p3 - p4) * s;
  G4ThreeVector c = p1 + (p4 - p1) * t;
  G4ThreeVector d = p2 + (p3 - p2) * t;
  G4ThreeVector normal = ((d - c).cross(b - a)).unit();
  if (normal == G4ThreeVector(0, 0, 0))
  {
    normal = (p1 + p2 + p3 + p4).unit();
  }
  return normal;
}

void CheckBoundingLimits(const G4VSolid* solid);
void CheckInside(const G4VSolid* solid);
void CheckSurfaceNormal(const G4VSolid* solid);
void CheckSafety(const G4VSolid* solid);
void CheckDistanceToIn(const G4VSolid* solid);
void CheckDistanceToOut(const G4VSolid* solid);
void CheckInOut(const G4VSolid* solid);

void CheckGetSurfaceArea(G4VSolid* solid);
void CheckGetCubicVolume(G4VSolid* solid);
void CheckGetPointOnSurface(G4VSolid* solid);

///////////////////////////////////////////////////////////////////////////////
//
// Unit test for Box
//
G4bool TestSolid(G4VSolid* solid)
{
  G4double dx_exp = 50, dy_exp = 60, dz_exp = 70;  // expected dimensions

  // Dump info
  //
  G4cout << G4endl;
  solid->DumpInfo();
  G4cout << G4endl;

  // Check getters
  //
  G4cout << "=== Check GetXHalfLength(), GetYHalfLength(), GetZHalfLength()" << G4endl;
  if (solid->GetEntityType() == "G4Box")
  {
    assert(((G4Box*)solid)->GetXHalfLength() == dx_exp);
    assert(((G4Box*)solid)->GetYHalfLength() == dy_exp);
    assert(((G4Box*)solid)->GetZHalfLength() == dz_exp);
  }
  else if (solid->GetEntityType() == "G4Para")
  {
    assert(((G4Para*)solid)->GetXHalfLength() == dx_exp);
    assert(((G4Para*)solid)->GetYHalfLength() == dy_exp);
    assert(((G4Para*)solid)->GetZHalfLength() == dz_exp);
  }
  else if (solid->GetEntityType() == "G4Trd")
  {
    assert(((G4Trd*)solid)->GetXHalfLength1() == dx_exp);
    assert(((G4Trd*)solid)->GetXHalfLength2() == dx_exp);
    assert(((G4Trd*)solid)->GetYHalfLength1() == dy_exp);
    assert(((G4Trd*)solid)->GetYHalfLength2() == dy_exp);
    assert(((G4Trd*)solid)->GetZHalfLength() == dz_exp);
  }
  else if (solid->GetEntityType() == "G4Trap")
  {
    assert(((G4Trap*)solid)->GetXHalfLength1() == dx_exp);
    assert(((G4Trap*)solid)->GetXHalfLength2() == dx_exp);
    assert(((G4Trap*)solid)->GetXHalfLength3() == dx_exp);
    assert(((G4Trap*)solid)->GetXHalfLength4() == dx_exp);
    assert(((G4Trap*)solid)->GetYHalfLength1() == dy_exp);
    assert(((G4Trap*)solid)->GetYHalfLength2() == dy_exp);
    assert(((G4Trap*)solid)->GetZHalfLength() == dz_exp);
  }

  // Check setters
  //
  G4cout << "=== Check SetXHalfLength(), SetYHalfLength(), SetZHalfLength()" << G4endl;
  G4double dx = 150;
  G4double dy = 160;
  G4double dz = 170;
  if (solid->GetEntityType() == "G4Box")
  {
    ((G4Box*)solid)->SetXHalfLength(dx);
    ((G4Box*)solid)->SetYHalfLength(dy);
    ((G4Box*)solid)->SetZHalfLength(dz);
    assert(((G4Box*)solid)->GetXHalfLength() == dx);
    assert(((G4Box*)solid)->GetYHalfLength() == dy);
    assert(((G4Box*)solid)->GetZHalfLength() == dz);
  }
  else if (solid->GetEntityType() == "G4Para")
  {
    ((G4Para*)solid)->SetXHalfLength(dx);
    ((G4Para*)solid)->SetYHalfLength(dy);
    ((G4Para*)solid)->SetZHalfLength(dz);
    assert(((G4Para*)solid)->GetXHalfLength() == dx);
    assert(((G4Para*)solid)->GetYHalfLength() == dy);
    assert(((G4Para*)solid)->GetZHalfLength() == dz);
  }
  else if (solid->GetEntityType() == "G4Trd")
  {
    ((G4Trd*)solid)->SetXHalfLength2(dx);
    ((G4Trd*)solid)->SetYHalfLength2(dy);
    ((G4Trd*)solid)->SetZHalfLength(dz);
    assert(((G4Trd*)solid)->GetXHalfLength2() == dx);
    assert(((G4Trd*)solid)->GetYHalfLength2() == dy);
    assert(((G4Trd*)solid)->GetZHalfLength() == dz);
  }
  else if (solid->GetEntityType() == "G4Trap")
  {
    ((G4Trap*)solid)->SetAllParameters(dz, 0, 0, dy, dx, dx, 0, dy, dx, dx, 0);
    assert(((G4Trap*)solid)->GetXHalfLength1() == dx);
    assert(((G4Trap*)solid)->GetXHalfLength2() == dx);
    assert(((G4Trap*)solid)->GetXHalfLength3() == dx);
    assert(((G4Trap*)solid)->GetXHalfLength4() == dx);
    assert(((G4Trap*)solid)->GetYHalfLength1() == dy);
    assert(((G4Trap*)solid)->GetYHalfLength2() == dy);
    assert(((G4Trap*)solid)->GetZHalfLength() == dz);
  }

  // Set original dimensions
  //
  if (solid->GetEntityType() == "G4Box")
  {
    ((G4Box*)solid)->SetXHalfLength(dx_exp);
    ((G4Box*)solid)->SetYHalfLength(dy_exp);
    ((G4Box*)solid)->SetZHalfLength(dz_exp);
  }
  else if (solid->GetEntityType() == "G4Para")
  {
    ((G4Para*)solid)->SetAllParameters(dx_exp, dy_exp, dz_exp, 0, 0, 0);
  }
  else if (solid->GetEntityType() == "G4Trd")
  {
    ((G4Trd*)solid)->SetAllParameters(dx_exp, dx_exp, dy_exp, dy_exp, dz_exp);
  }
  else if (solid->GetEntityType() == "G4Trap")
  {
    ((G4Trap*)solid)
      ->SetAllParameters(dz_exp, 0, 0, dy_exp, dx_exp, dx_exp, 0, dy_exp, dx_exp, dx_exp, 0);
  }

  // Check main methods
  //
  G4cout << "=== Check BoundingLimits()" << G4endl;
  CheckBoundingLimits(solid);

  G4cout << "=== Check Inside(p)" << G4endl;
  CheckInside(solid);

  G4cout << "=== Check SurfaceNormal(p)" << G4endl;
  CheckSurfaceNormal(solid);

  G4cout << "=== Check SafetyToIn(p), SafetyToOut(p)" << G4endl;
  CheckSafety(solid);

  G4cout << "=== Check DistanceToIn(p,v)" << G4endl;
  CheckDistanceToIn(solid);

  G4cout << "=== Check DistanceToOut(p,v)" << G4endl;
  CheckDistanceToOut(solid);

  G4cout << "=== Check DistanceToIn(p,v) vs DistanceToOut(p,v)" << G4endl;
  CheckInOut(solid);

  G4cout << "=== Check GetSurfaceArea()" << G4endl;
  CheckGetSurfaceArea(solid);

  G4cout << "=== Check GetCubicVolume()" << G4endl;
  CheckGetCubicVolume(solid);

  G4cout << "=== Check GetPointOnSurface()" << G4endl;
  CheckGetPointOnSurface(solid);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//
// Check BoudingLimits()
//
void CheckBoundingLimits(const G4VSolid* solid)
{
  G4ThreeVector pmin, pmax;
  solid->BoundingLimits(pmin, pmax);

  G4ThreeVector minCheck(kInfinity, kInfinity, kInfinity);
  G4ThreeVector maxCheck(-kInfinity, -kInfinity, -kInfinity);
  G4double x, y, z;
  for (int i = 0; i < 1000000; ++i)
  {
    G4ThreeVector p = solid->GetPointOnSurface();
    x = std::min(p.x(), minCheck.x());
    y = std::min(p.y(), minCheck.y());
    z = std::min(p.z(), minCheck.z());
    minCheck.set(x, y, z);
    x = std::max(p.x(), maxCheck.x());
    y = std::max(p.y(), maxCheck.y());
    z = std::max(p.z(), maxCheck.z());
    maxCheck.set(x, y, z);
  }
  G4double tolerance;
  assert(ApproxEqual(pmin, minCheck, tolerance = 0));
  assert(ApproxEqual(pmax, maxCheck, tolerance = 0));
}

///////////////////////////////////////////////////////////////////////////////
//
// Check Inside(p)
//
void CheckInside(const G4VSolid* solid)
{
  G4double tolerance;
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {4, 5, 6, 7}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  // check that gravity center is inside
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 8; ++i)
  {
    pcenter += vv[i] / 8.;
  }
  assert(solid->Inside(pcenter) == kInside);
  assert(ApproxEqual(pcenter, G4ThreeVector(0, 0, 0), tolerance = 1.e-15));

  // check positions of generated points
  auto Ns = 5, Nt = 5, k = 0;
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    for (auto is = k; is < Ns - k; ++is)
    {
      G4double s = is * sdel;
      for (auto it = k; it < Nt - k; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector shift = (p - pcenter).unit();
        assert(solid->Inside(p) == kSurface);
        assert(solid->Inside(p - 1 * shift) == kInside);
        assert(solid->Inside(p + 1 * shift) == kOutside);
        assert(solid->Inside(p + 10 * shift) == kOutside);
        assert(solid->Inside(p + 100 * shift) == kOutside);
        assert(solid->Inside(p + 1000 * shift) == kOutside);
        assert(solid->Inside(p + 10000 * shift) == kOutside);
        assert(solid->Inside(p + 0.25 * kTolerance * shift) == kSurface);
        assert(solid->Inside(p - 0.25 * kTolerance * shift) == kSurface);
        assert(solid->Inside(p + 0.499 * kTolerance * shift) == kSurface);
        assert(solid->Inside(p - 0.499 * kTolerance * shift) == kSurface);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetSurfaceNormal(p)
//
void CheckSurfaceNormal(const G4VSolid* solid)
{
  G4double tolerance;
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  G4int iedge[12][4] = {{0, 1, 0, 1}, {1, 2, 0, 2}, {2, 3, 0, 3}, {3, 0, 0, 4},
                        {0, 4, 4, 1}, {1, 5, 1, 2}, {2, 6, 2, 3}, {3, 7, 3, 4},
                        {4, 5, 5, 1}, {5, 6, 5, 2}, {6, 7, 5, 3}, {7, 4, 5, 4}};
  G4int inode[8][4] = {{0, 0, 4, 1}, {1, 0, 1, 2}, {2, 0, 2, 3}, {3, 0, 3, 4},
                       {4, 5, 4, 1}, {5, 5, 1, 2}, {6, 5, 2, 3}, {7, 5, 3, 4}};

  // set 3d vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  // check that gravity center is inside and at the origin
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 8; ++i)
  {
    pcenter += vv[i] / 8.;
  }
  assert(solid->Inside(pcenter) == kInside);
  assert(ApproxEqual(pcenter, G4ThreeVector(0, 0, 0), tolerance = 1.e-15));

  // check normals on faces
  auto Ns = 5, Nt = 5;
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    for (auto is = 1; is < Ns - 1; ++is)
    {
      G4double s = is * sdel;
      for (auto it = 1; it < Nt - 1; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector n = getNormal(p0, p1, p2, p3, s, t);
        G4ThreeVector shift = (p - pcenter).unit();

        G4ThreeVector pshift = p + 0. * shift;
        ;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        pshift = p + 0.25 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        pshift = p - 0.25 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        pshift = p + 0.499 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        pshift = p - 0.499 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        // check ApproxSurfaceNormal(), point is not on the surface
        pshift = p - 1000. * kTolerance * shift;
        assert(solid->Inside(pshift) != kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));

        pshift = p + 1000. * kTolerance * shift;
        assert(solid->Inside(pshift) != kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), tolerance = 1.e-15));
      }
    }
  }
  // check normals at edges
  for (auto i = 0; i < 12; ++i)
  {
    G4ThreeVector a0 = vv[iface[iedge[i][2]][0]];
    G4ThreeVector a1 = vv[iface[iedge[i][2]][1]];
    G4ThreeVector a2 = vv[iface[iedge[i][2]][2]];
    G4ThreeVector a3 = vv[iface[iedge[i][2]][3]];
    G4ThreeVector n1 = getNormal(a0, a1, a2, a3, 0.5, 0.5);

    G4ThreeVector b0 = vv[iface[iedge[i][3]][0]];
    G4ThreeVector b1 = vv[iface[iedge[i][3]][1]];
    G4ThreeVector b2 = vv[iface[iedge[i][3]][2]];
    G4ThreeVector b3 = vv[iface[iedge[i][3]][3]];
    G4ThreeVector n2 = getNormal(b0, b1, b2, b3, 0.5, 0.5);

    G4ThreeVector p = 0.5 * (vv[iedge[i][0]] + vv[iedge[i][1]]);
    G4ThreeVector n = (n1 + n2).unit();
    assert(solid->Inside(p) == kSurface);
    assert(ApproxEqual(n, solid->SurfaceNormal(p), tolerance = 0));
  }
  // check normals at vertices
  for (auto i = 0; i < 8; ++i)
  {
    G4ThreeVector a0 = vv[iface[inode[i][1]][0]];
    G4ThreeVector a1 = vv[iface[inode[i][1]][1]];
    G4ThreeVector a2 = vv[iface[inode[i][1]][2]];
    G4ThreeVector a3 = vv[iface[inode[i][1]][3]];
    G4ThreeVector n1 = getNormal(a0, a1, a2, a3, 0.5, 0.5);

    G4ThreeVector b0 = vv[iface[inode[i][2]][0]];
    G4ThreeVector b1 = vv[iface[inode[i][2]][1]];
    G4ThreeVector b2 = vv[iface[inode[i][2]][2]];
    G4ThreeVector b3 = vv[iface[inode[i][2]][3]];
    G4ThreeVector n2 = getNormal(b0, b1, b2, b3, 0.5, 0.5);

    G4ThreeVector c0 = vv[iface[inode[i][3]][0]];
    G4ThreeVector c1 = vv[iface[inode[i][3]][1]];
    G4ThreeVector c2 = vv[iface[inode[i][3]][2]];
    G4ThreeVector c3 = vv[iface[inode[i][3]][3]];
    G4ThreeVector n3 = getNormal(c0, c1, c2, c3, 0.5, 0.5);

    G4ThreeVector p = vv[inode[i][0]];
    G4ThreeVector n = (n1 + n2 + n3).unit();
    assert(solid->Inside(p) == kSurface);
    assert(ApproxEqual(n, solid->SurfaceNormal(p), tolerance = 1.e-15));
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check safety (DistanceToIn(p)/DistanceToOut(p))
//
void CheckSafety(const G4VSolid* solid)
{
  G4double tolerance;
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  // set points on the surface
  auto Ns = 5, Nt = 5, nsurf = 6 * Ns * Nt;
  std::vector<G4ThreeVector> psurf(nsurf);
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  auto ipoint = 0;
  for (auto ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    for (auto is = 0; is < Ns; ++is)
    {
      G4double s = is * sdel;
      for (auto it = 0; it < Nt; ++it)
      {
        G4double t = it * tdel;
        psurf.at(ipoint) = getPoint(p0, p1, p2, p3, s, t);
        assert(solid->Inside(psurf[ipoint]) == kSurface);
        ipoint++;
      }
    }
  }

  // Check safety distance for points on surface
  G4double shift = 0.25 * kTolerance;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToIn(p) >= 0 && solid->DistanceToIn(p) < 1.e-15);
    assert(solid->DistanceToOut(p) >= 0 && solid->DistanceToOut(p) < 1.e-15);

    p = psurf[i] - norm * shift;  // shift slightly inside
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToIn(p) == 0);
    assert(solid->DistanceToOut(p) >= 0 && solid->DistanceToOut(p) < shift * (1. + 1.e-15));

    p = psurf[i] + norm * shift;  // shift slightly outside
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToIn(p) >= 0 && solid->DistanceToIn(p) < shift * (1. + 1.e-15));
    assert(solid->DistanceToOut(p) == 0);
  }

  // Check safety distance for points shifted by 1 mm to inside
  shift = 1. * mm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] - norm * shift;
    assert(solid->Inside(p) == kInside);
    assert(solid->DistanceToOut(p) > 0);
    assert(solid->DistanceToOut(p) < shift * (1. + 1.e-15));
  }

  // Check safety distance for points shifted by 10 mm to inside
  shift = 10. * mm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] - norm * shift;
    assert(solid->Inside(p) == kInside);
    assert(solid->DistanceToOut(p) > 0);
    assert(solid->DistanceToOut(p) < shift * (1. + 1.e-15));
  }

  // Check safety distance for points shifted by 1, 10, 100 cm to outside
  shift = 1. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift * (1. + 1.e-15));
  }

  shift = 10. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift * (1. + 1.e-15));
  }

  shift = 100. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift * (1. + 1.e-15));
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToIn(p,v)
//
void CheckDistanceToIn(const G4VSolid* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  // generate points on the surface and compute normals
  auto Ns = 5, Nt = 5, npoints = 6 * Ns * Nt;
  std::vector<G4ThreeVector> pp(npoints);
  std::vector<G4ThreeVector> nn(npoints);
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  auto ipoint = 0;
  for (auto ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    for (auto is = 0; is < Ns; ++is)
    {
      G4double s = is * sdel;
      for (auto it = 0; it < Nt; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        assert(solid->Inside(p) == kSurface);
        pp.at(ipoint) = p;
        nn.at(ipoint) = solid->SurfaceNormal(p);
        ipoint++;
      }
    }
  }

  // check points on surface
  for (auto i = 0; i < npoints; ++i)
  {
    G4ThreeVector n = nn[i];
    G4ThreeVector p = pp[i];
    assert(solid->DistanceToIn(p, n) == kInfinity);
    assert(solid->DistanceToIn(p, -n) == 0);

    assert(solid->DistanceToIn(p + 0.25 * kTolerance * n, n) == kInfinity);
    assert(solid->DistanceToIn(p + 0.25 * kTolerance * n, -n) == 0);

    assert(solid->DistanceToIn(p + 0.499 * kTolerance * n, n) == kInfinity);
    assert(solid->DistanceToIn(p + 0.499 * kTolerance * n, -n) == 0);

    assert(solid->DistanceToIn(p - 0.25 * kTolerance * n, n) == kInfinity);
    assert(solid->DistanceToIn(p - 0.25 * kTolerance * n, -n) == 0);

    assert(solid->DistanceToIn(p - 0.499 * kTolerance * n, n) == kInfinity);
    assert(solid->DistanceToIn(p - 0.499 * kTolerance * n, -n) == 0);
  }

  // check points shifted to outside
  for (auto i = 0; i < npoints; ++i)
  {
    G4ThreeVector n = nn[i];
    G4ThreeVector p = pp[i] + 1 * m * n;
    assert(solid->DistanceToIn(p, n) == kInfinity);
    assert(ApproxEqual(solid->DistanceToIn(p, -n), 1 * m));

    p = pp[i] + 100 * m * n;
    assert(solid->DistanceToIn(p, n) == kInfinity);
    assert(ApproxEqual(solid->DistanceToIn(p, -n), 100 * m));
  }

  // check scratching/missing
  for (auto ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    assert(solid->Inside(pc) == kSurface);
    G4ThreeVector n = solid->SurfaceNormal(pc);  // get normal

    G4ThreeVector pp[4] = {p0, p1, p2, p3};
    G4ThreeVector p, v;
    for (auto i = 0; i < 4; ++i)
    {
      auto shift = 100. * mm;
      p = pc + shift * n;  // move the point 10 cm away from the face
      v = (p0 - p).unit();  // set direction towards the vertex
      CheckScratching(1, solid->DistanceToIn(p, v));

      p = pc + (shift - 0.25 * kTolerance) * n;  // move the point slightly towards the face
      CheckScratching(2, solid->DistanceToIn(p, v));  // scratching, expected that dist = kInfinity

      p = pc + (shift + 10 * mm) * n;  // move the point a bit more away from the face
      assert(solid->DistanceToIn(p, v)
             == kInfinity);  // tracks do not hit the solid, dist = kInfinity
    }
  }

  // check directions along the faces
  for (auto ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    G4ThreeVector n = solid->SurfaceNormal(pc);  // get normal

    G4ThreeVector pp[4] = {p0, p1, p2, p3};
    G4ThreeVector p, v;
    for (auto i = 0; i < 4; ++i)
    {
      p = pc;  // put the point in center of the face
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      if (solid->DistanceToIn(p, v) != kInfinity)
      {
        G4cout << " DistanceToIn failed, case #1:   dist = " << solid->DistanceToIn(p, v)
               << ", expected kInfinity   (distToOut = " << solid->DistanceToOut(p, v) << ")"
               << G4endl;
      }

      p = pc - 0.25 * kTolerance * n;  // move the point slightly to inside
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      if (solid->DistanceToIn(p, v) != kInfinity)
      {
        G4cout << " DistanceToIn failed, case #2:   dist = " << solid->DistanceToIn(p, v)
               << ", expected kInfinity   (distToOut = " << solid->DistanceToOut(p, v) << ")"
               << G4endl;
      }

      p = pc + 0.25 * kTolerance * n;  // move the point slightly to outside
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      CheckScratching(3, solid->DistanceToIn(p, v));
      v = ((pp[i] + pc) / 2 - p).unit();  // set the direction to the face interior
      ApproxEqual(solid->DistanceToIn(p, v), ((pp[i] + pc) / 2 - p).mag());

      p = pc + 10 * (pp[i] - pc);  // move the point outside of the solid in the face plane
      assert(solid->Inside(p) == kOutside);
      v = (pc - pp[i]).unit();
      if (solid->DistanceToIn(p, v) != kInfinity)
      {
        G4cout << " DistanceToIn failed, case #3:   dist = " << solid->DistanceToIn(p, v)
               << ", expected kInfinity" << G4endl;
      }

      p = pc + 10 * (pp[i] - pc)
          + 0.25 * kTolerance * n;  // put the point slightly above the face plane
      v = (pp[i] - p).unit();
      ApproxEqual(solid->DistanceToIn(p, v), (pp[i] - p).mag());
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToOut(p,v)
//
void CheckDistanceToOut(const G4VSolid* solid)
{
  G4ThreeVector p, n, v;
  G4double tolerance;
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);
  G4ThreeVector origin = (vv[0] + vv[1] + vv[2] + vv[3] + vv[4] + vv[5] + vv[6] + vv[7]) / 8;
  ApproxEqual(origin, G4ThreeVector(0, 0, 0), tolerance = 1.e-15);

  // generate points on the surface and compute normals
  G4int Ns = 5, Nt = 5, npoints = 6 * Ns * Nt;
  std::vector<G4ThreeVector> pp(npoints);
  std::vector<G4ThreeVector> nn(npoints);
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  for (G4int ipoint = 0, ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    for (auto is = 0; is < Ns; ++is)
    {
      for (auto it = 0; it < Nt; ++it)
      {
        G4double s = is * sdel;
        G4double t = it * tdel;
        p = getPoint(p0, p1, p2, p3, s, t);
        assert(solid->Inside(p) == kSurface);
        pp.at(ipoint) = p;
        nn.at(ipoint) = solid->SurfaceNormal(p);
        ipoint++;
      }
    }
  }

  // check calculation of return normal
  for (G4int ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    for (auto is = 1; is < Ns - 1; ++is)
    {
      for (auto it = 1; it < Nt - 1; ++it)
      {
        G4double s = is * sdel;
        G4double t = it * tdel;
        p = getPoint(p0, p1, p2, p3, s, t);
        n = getNormal(p0, p1, p2, p3, s, t);
        v = (p - origin).unit();
        G4bool validNorm;
        G4ThreeVector norm;
        ApproxEqual(solid->DistanceToOut(origin, v, true, &validNorm, &norm), (p - origin).mag());
        ApproxEqual(norm, n, tolerance = 1.e-15);
        assert(validNorm);
      }
    }
  }

  // check points on surface
  for (auto i = 0; i < npoints; ++i)
  {
    p = pp[i];
    n = nn[i];
    assert(solid->DistanceToOut(p, n) == 0);
    assert(solid->DistanceToOut(p, -n) > 0);

    assert(solid->DistanceToOut(p + 0.25 * kTolerance * n, n) == 0);
    assert(solid->DistanceToOut(p + 0.25 * kTolerance * n, -n) > 0);

    assert(solid->DistanceToOut(p + 0.499 * kTolerance * n, n) == 0);
    assert(solid->DistanceToOut(p + 0.499 * kTolerance * n, -n) > 0);

    assert(solid->DistanceToOut(p - 0.25 * kTolerance * n, n) == 0);
    assert(solid->DistanceToOut(p - 0.25 * kTolerance * n, -n) > 0);

    assert(solid->DistanceToOut(p - 0.499 * kTolerance * n, n) == 0);
    assert(solid->DistanceToOut(p - 0.499 * kTolerance * n, -n) > 0);
  }

  // check inside points
  for (G4int i = 0; i < npoints; ++i)
  {
    G4double shift = 1 * mm;
    n = nn[i];
    p = pp[i] - shift * n;
    ApproxEqual(solid->DistanceToOut(p, n), shift);

    p = pp[i];
    v = (p - origin).unit();
    ApproxEqual(solid->DistanceToOut(origin, v), (p - origin).mag());
  }

  // check directions along the faces
  for (auto ifa = 0; ifa < 6; ++ifa)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    G4ThreeVector n = solid->SurfaceNormal(pc);  // get normal

    G4ThreeVector pp[4] = {p0, p1, p2, p3};
    G4ThreeVector p, v;
    for (auto i = 0; i < 4; ++i)
    {
      p = pc;  // put the point in center of the face
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      ApproxEqual(solid->DistanceToOut(p, v), (pp[i] - p).mag());

      p = pc - 0.25 * kTolerance * n;  // move the point slightly to inside
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      if (solid->DistanceToOut(p, v) != 0)  // point is leaving the solid
      {
        G4cout << " DistanceToOut failed, case #1:   dist = " << solid->DistanceToOut(p, v)
               << ", expected 0   (distToIn = " << solid->DistanceToIn(p, v) << ")" << G4endl;
      }

      p = pc + 0.25 * kTolerance * n;  // move the point slightly to outside
      if (solid->DistanceToOut(p, v) != 0)  // point is leaving the solid
      {
        G4cout << " DistanceToOut failed, case #2:   dist = " << solid->DistanceToOut(p, v)
               << ", expected 0   (distToIn = " << solid->DistanceToIn(p, v) << ")" << G4endl;
      }

      p = pc + 0.499 * kTolerance * n;  // move the point slightly to outside
      v = (pp[i] - p).unit();  // set the direction towards the vertex
      ApproxEqual(solid->DistanceToOut(p, v), (pp[i] - p).mag());

      v = ((pp[i] + pc) / 2 - p).unit();  // set the direction to the face interior
      assert(solid->DistanceToOut(p, v) > ((pp[i] + pc) / 2 - p).mag());
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToIn(p.v) vs DistanceToOut(p,v)
// Checks directions through points on the surface,
// distIn and distOut should never be equal to zero at the same time
//
void CheckInOut(const G4VSolid* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  G4int Ns = 5, Nt = 5, npoints = 6 * Ns * Nt;
  std::vector<G4ThreeVector> pp(npoints);
  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);

  G4ThreeVector p, v;
  G4double distIn, distOut;

  // loop along faces
  //
  for (G4int ifa = 0; ifa < 6; ++ifa)
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    G4ThreeVector n = solid->SurfaceNormal(pc);  // face normal

    // generate points on the surface of the face
    for (G4int ipoint = 0, is = 0; is < Ns; ++is)
    {
      for (auto it = 0; it < Nt; ++it)
      {
        G4double s = is * sdel;
        G4double t = it * tdel;
        pp.at(ipoint++) = getPoint(p0, p1, p2, p3, s, t);
      }
    }

    // set directions
    G4ThreeVector node[4] = {p0, p1, p2, p3};
    std::vector<G4ThreeVector> vplane(8);  // directions along the face
    std::vector<G4ThreeVector> vup(8);  // directions going slightly to outside
    std::vector<G4ThreeVector> vdown(8);  // directions going slightly to inside
    for (G4int ivect = 0, i = 0; i < 4; ++i)
    {
      p = node[i];  // vertex
      vplane.at(ivect) = (p - pc).unit();
      vup.at(ivect) = (p - pc + 0.25 * kTolerance * n).unit();
      vdown.at(ivect) = (p - pc - 0.25 * kTolerance * n).unit();
      ivect++;

      p = (node[i] + node[(i + 1) % 4]) / 2;  // middle of the edge
      vplane.at(ivect) = (p - pc).unit();
      vup.at(ivect) = (p - pc + 0.25 * kTolerance * n).unit();
      vdown.at(ivect) = (p - pc - 0.25 * kTolerance * n).unit();
      ivect++;
    }

    for (G4int ipoint = 0; ipoint < npoints; ++ipoint)
    {
      p = pp[ipoint];
      for (G4int i = 0; i < 8; ++i)
      {
        v = vplane[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vup[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vdown[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);
      }

      p = pp[ipoint] + 0.25 * kTolerance * n;
      for (G4int i = 0; i < 8; ++i)
      {
        v = vplane[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vup[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vdown[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);
      }

      p = pp[ipoint] - 0.25 * kTolerance * n;
      for (G4int i = 0; i < 8; ++i)
      {
        v = vplane[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vup[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);

        v = vdown[i];
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn >= 0 && distOut >= 0 && (distIn + distOut) > 0);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetSurfaceArea()
//
void CheckGetSurfaceArea(G4VSolid* solid)
{
  G4double area = solid->GetSurfaceArea();
  G4double areaCheck = solid->G4VSolid::GetSurfaceArea();
  G4double diff = (area - areaCheck) / area;
  G4cout << "  " << solid->GetName() << " area = " << area << ", mc_estimated = " << areaCheck
         << ", diff = " << 100 * diff << "% ";
  G4cout << ((std::abs(diff) < 0.01) ? "(OK, diff < 1% as expected)"
                                     : "(failed, expected diff < 1%)")
         << G4endl;
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetCubicVolume()
//
void CheckGetCubicVolume(G4VSolid* solid)
{
  G4double volume = solid->GetCubicVolume();
  G4double volumeCheck = solid->G4VSolid::GetCubicVolume();
  G4double diff = (volume - volumeCheck) / volume;
  G4cout << "  " << solid->GetName() << " volume = " << volume << ", mc_estimated = " << volumeCheck
         << ", diff = " << 100 * diff / volume << "% ";
  G4cout << ((std::abs(diff) < 0.01) ? "(OK, diff < 1% as expected)"
                                     : "(failed, expected diff < 1%)")
         << G4endl;
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetPointOnSurface()
//
void CheckGetPointOnSurface(G4VSolid* solid)
{
  clock_t t;
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // get vertices
  std::vector<G4ThreeVector> vv(8);
  getVertices(solid, vv);

  // get face normals and areas
  std::vector<G4ThreeVector> nn(6);
  G4double ff[6];
  for (G4int i = 0; i < 6; ++i)
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    G4ThreeVector p = (p0 + p1 + p2 + p3) / 4.;
    assert(solid->Inside(p) == kSurface);
    nn[i] = solid->SurfaceNormal(p);
    ff[i] = G4GeomTools::QuadAreaNormal(p0, p1, p2, p3).mag();
  }
  G4double TotalArea = solid->GetSurfaceArea();

  // pick up random points
  G4int np[6] = {0}, nfactor = 100, ntot = TotalArea * nfactor;
  t = clock();
  for (G4int i = 0; i < ntot; ++i)
  {
    G4ThreeVector p = solid->GetPointOnSurface();
    assert(solid->Inside(p) == kSurface);
    G4ThreeVector n = solid->SurfaceNormal(p);
    np[0] += (n == nn[0]);
    np[1] += (n == nn[1]);
    np[2] += (n == nn[2]);
    np[3] += (n == nn[3]);
    np[4] += (n == nn[4]);
    np[5] += (n == nn[5]);
  }
  t = clock() - t;
  G4int other = ntot - (np[0] + np[1] + np[2] + np[3] + np[4] + np[5]);

  // print out info
  G4cout << solid->GetName() << " time : " << 1000 * (double)t / CLOCKS_PER_SEC << " ms;  "
         << ntot / 1000000. << " million points;"
         << " time per million points: " << 1000 * ((double)t / CLOCKS_PER_SEC) * 1000000. / ntot
         << " ms" << G4endl;

  G4cout << "  Area : bottom, face1, face2, face3, face4, top = " << ff[0] << ",   " << ff[1]
         << ",   " << ff[2] << ",   " << ff[3] << ",   " << ff[4] << ",   " << ff[5]
         << "   Total = " << TotalArea << G4endl;

  G4cout << "  N.pnt: bottom, face1, face2, face3, face4, top = " << np[0] << ", " << np[1] << ", "
         << np[2] << ", " << np[3] << ", " << np[4] << ", " << np[5] << " Total = " << ntot
         << " (others = " << other << ")" << G4endl;
}

///////////////////////////////////////////////////////////////////////////////
//
// Main
//
int main(int argc, char* argv[])
{
  G4double dx = 50, dy = 60, dz = 70;

  // Uncomment/comment the lines below
  // to include/exclude required solid(s) to/from testing
  //
  assert(TestSolid(new G4Box("TestBox", dx, dy, dz)));
  // assert(TestSolid(new G4Para("TestPara", dx, dy, dz, 0, 0, 0)));
  // assert(TestSolid(new G4Trd("TestTrd", dx, dx, dy, dy, dz)));
  // assert(TestSolid(new G4Trap("TestTrap", dx, dx, dy, dy, dz)));

  G4cout << "\n    All tests passed!\n" << G4endl;
  return 0;
}
