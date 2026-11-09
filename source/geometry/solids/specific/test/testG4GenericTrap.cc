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
// Unit test for G4GenericTrap
//
// Author: E.Tcherniaev (evgueni.tcherniaev@cern.ch)
//

#include "G4GenericTrap.hh"
#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Timer.hh"
#include "globals.hh"

#include "geomdefs.hh"
#include <assert.h>

#include <cmath>
#include <iomanip>

G4bool ApproxEqual(G4double a, const G4double b)
{
  return std::abs(a - b) < 1e-12;
}

G4bool ApproxEqual(const G4ThreeVector& a, const G4ThreeVector& b, G4double tol)
{
  return (a - b).mag() < tol;
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

G4bool isPlane(const G4ThreeVector& p0, const G4ThreeVector& p1, const G4ThreeVector& p2,
               const G4ThreeVector& p3)
{
  G4double volume = std::abs((p1 - p0).cross(p2 - p0).dot(p3 - p0));
  return (volume < 1e-12);
}

void CheckBoundingLimits(const G4GenericTrap* solid);
void CheckInside(const G4GenericTrap* solid);
void CheckSurfaceNormal(const G4GenericTrap* solid);
void CheckSafety(const G4GenericTrap* solid);
void CheckDistanceToIn(const G4GenericTrap* solid);
void CheckDistanceToIn_special(const G4GenericTrap* solid);
void CheckDistanceToOut(const G4GenericTrap* solid);
void CheckInOut(const G4GenericTrap* solid);

void CheckGetSurfaceArea(G4GenericTrap* solid);
void CheckGetCubicVolume(G4GenericTrap* solid);
void CheckGetPointOnSurface(G4GenericTrap* solid);

///////////////////////////////////////////////////////////////////////////////
//
// Unit test for Generic trap
//
G4bool TestGenericTrap()
{
  G4GenericTrap* prism = nullptr;
  G4GenericTrap* pyram = nullptr;
  G4GenericTrap* tetra = nullptr;
  G4GenericTrap* solid = nullptr;
  G4GenericTrap* _emec = nullptr;
  G4double dz_prism, dz_pyram, dz_tetra, dz_solid, dz_emec;
  std::vector<G4TwoVector> vv(8);

  // Construct test solids
  //
  vv[0].set(0, -50);  // triangular prism
  vv[1].set(-50, 0);
  vv[2].set(-50, 0);
  vv[3].set(50, 50);
  vv[4].set(0, -30);
  vv[5].set(-30, 0);
  vv[6].set(-30, 0);
  vv[7].set(30, 30);
  prism = new G4GenericTrap("TestPrism", dz_prism = 40, vv);

  vv[0].set(20, 10);  // pyramid
  vv[1].set(20, 10);
  vv[2].set(20, 10);
  vv[3].set(20, 10);
  vv[4].set(-60, -50);
  vv[5].set(-40, 50);
  vv[6].set(60, 50);
  vv[7].set(40, -50);
  pyram = new G4GenericTrap("TestPyramid", dz_pyram = 50, vv);

  vv[0].set(-60, 0);  // tetrahedron
  vv[1].set(50, 0);
  vv[2].set(50, 0);
  vv[3].set(-60, 0);
  vv[4].set(0, 60);
  vv[5].set(0, 60);
  vv[6].set(0, -40);
  vv[7].set(0, -40);
  tetra = new G4GenericTrap("TestTetrahedron", dz_tetra = 50, vv);

  vv[0].set(-70, -60);  // generic trap
  vv[1].set(-60, 60);
  vv[2].set(60, 60);
  vv[3].set(70, -60);
  vv[4].set(0, -60);
  vv[5].set(-60, 0);
  vv[6].set(0, 50);
  vv[7].set(60, 0);
  solid = new G4GenericTrap("TestGenTrap", dz_solid = 50, vv);

  vv[0].set(-5.535915669346718, 616.0522317995544);  // EMEC
  vv[1].set(-8.069729892664068, 616.0522317995544);
  vv[2].set(-24.39831285725495, 2008.457171565301);
  vv[3].set(-20.137911534789, 2008.457171565301);
  vv[4].set(8.12464736753118, 620.7352691121155);
  vv[5].set(5.585026161217877, 620.7352691121155);
  vv[6].set(20.29841347210058, 2023.762181401904);
  vv[7].set(24.57779306143527, 2023.762181401904);
  _emec = new G4GenericTrap("Test_EMEC", dz_emec = 14.16666666666667, vv);

  // Check accessors
  //
  std::cout << "=== Check GetZHalfLength()" << std::endl;
  assert(prism->GetZHalfLength() == dz_prism);
  assert(pyram->GetZHalfLength() == dz_pyram);
  assert(tetra->GetZHalfLength() == dz_tetra);
  assert(solid->GetZHalfLength() == dz_solid);
  assert(_emec->GetZHalfLength() == dz_emec);

  std::cout << "=== Check GetNofVertices()" << std::endl;
  assert(prism->GetNofVertices() == 8);
  assert(pyram->GetNofVertices() == 8);
  assert(tetra->GetNofVertices() == 8);
  assert(solid->GetNofVertices() == 8);
  assert(_emec->GetNofVertices() == 8);

  std::cout << "=== Check GetVertex(), GetVertices()" << std::endl;
  for (auto i = 0; i < 8; ++i)
  {
    assert(prism->GetVertices().at(i) == prism->GetVertex(i));
    assert(pyram->GetVertices().at(i) == pyram->GetVertex(i));
    assert(tetra->GetVertices().at(i) == tetra->GetVertex(i));
    assert(solid->GetVertices().at(i) == solid->GetVertex(i));
    assert(_emec->GetVertices().at(i) == _emec->GetVertex(i));
  }

  std::cout << "=== Check IsTwisted(), GetTwistAngle()" << std::endl;
  assert(!prism->IsTwisted());
  assert(!pyram->IsTwisted());
  assert(!tetra->IsTwisted());
  assert(solid->IsTwisted());
  G4int ntwist = 0;
  for (auto i = 0; i < 4; ++i)  // prism, pyram, tetra - no twisted faces

  {
    ntwist += (prism->GetTwistAngle(i) != 0);
    ntwist += (pyram->GetTwistAngle(i) != 0);
    ntwist += (tetra->GetTwistAngle(i) != 0);
  }
  assert(ntwist == 0);
  for (auto i = 0; i < 4; ++i)  // gentrap - four twisted faces
  {
    assert(solid->GetTwistAngle(i) > 0);
    assert(solid->GetTwistAngle(i) < halfpi);
    ntwist += (solid->GetTwistAngle(i) != 0);
  }
  assert(ntwist == 4);
  ntwist = 0;
  for (auto i = 0; i < 4; ++i)  // emec - two twisted faces
  {
    assert(_emec->GetTwistAngle(i) <= 0);  // twsit angle for plane faces is 0
    assert(_emec->GetTwistAngle(i) > -halfpi);
    ntwist += (_emec->GetTwistAngle(i) != 0);
  }
  assert(ntwist == 2);

  // Check main methods
  //
  std::cout << "=== Check BoundingLimits()" << std::endl;
  CheckBoundingLimits(prism);
  CheckBoundingLimits(pyram);
  CheckBoundingLimits(tetra);
  CheckBoundingLimits(solid);
  CheckBoundingLimits(_emec);

  std::cout << "=== Check Inside(p)" << std::endl;
  CheckInside(prism);
  CheckInside(pyram);
  CheckInside(tetra);
  CheckInside(solid);
  CheckInside(_emec);

  std::cout << "=== Check SurfaceNormal(p)" << std::endl;
  CheckSurfaceNormal(prism);
  CheckSurfaceNormal(pyram);
  CheckSurfaceNormal(tetra);
  CheckSurfaceNormal(solid);
  CheckSurfaceNormal(_emec);

  std::cout << "=== Check SafetyToIn(p), SafetyToOut(p)" << std::endl;
  CheckSafety(prism);
  CheckSafety(pyram);
  CheckSafety(tetra);
  CheckSafety(solid);
  CheckSafety(_emec);

  std::cout << "=== Check DistanceToIn(p,v)" << std::endl;
  CheckDistanceToIn(prism);
  CheckDistanceToIn(pyram);
  CheckDistanceToIn(tetra);
  CheckDistanceToIn(solid);
  CheckDistanceToIn(_emec);

  std::cout << "=== Check DistanceToIn_special(p,v)" << std::endl;
  CheckDistanceToIn_special(prism);
  CheckDistanceToIn_special(pyram);
  CheckDistanceToIn_special(tetra);
  CheckDistanceToIn_special(solid);
  CheckDistanceToIn_special(_emec);

  std::cout << "=== Check DistanceToOut(p,v)" << std::endl;
  CheckDistanceToOut(prism);
  CheckDistanceToOut(pyram);
  CheckDistanceToOut(tetra);
  CheckDistanceToOut(solid);
  CheckDistanceToOut(_emec);

  std::cout << "=== Check DistanceToIn(p,v) vs DistanceToOut(p,v)" << std::endl;
  CheckInOut(prism);
  CheckInOut(pyram);
  CheckInOut(tetra);
  CheckInOut(solid);
  CheckInOut(_emec);

  std::cout << "=== Check GetSurfaceArea()" << std::endl;
  CheckGetSurfaceArea(prism);
  CheckGetSurfaceArea(pyram);
  CheckGetSurfaceArea(tetra);
  CheckGetSurfaceArea(solid);
  CheckGetSurfaceArea(_emec);

  std::cout << "=== Check GetCubicVolume()" << std::endl;
  CheckGetCubicVolume(prism);
  CheckGetCubicVolume(pyram);
  CheckGetCubicVolume(tetra);
  CheckGetCubicVolume(solid);
  CheckGetCubicVolume(_emec);

  std::cout << "=== Check GetPointOnSurface()" << std::endl;
  CheckGetPointOnSurface(prism);
  CheckGetPointOnSurface(pyram);
  CheckGetPointOnSurface(tetra);
  CheckGetPointOnSurface(solid);
  CheckGetPointOnSurface(_emec);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
//
// Check BoudingLimits()
//
void CheckBoundingLimits(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kHalfTolerance = 0.5 * kTolerance;

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
  assert(pmin.x() - kHalfTolerance < minCheck.x());
  assert(pmin.y() - kHalfTolerance < minCheck.y());
  assert(pmin.z() - kHalfTolerance < minCheck.z());
  assert(pmax.x() + kHalfTolerance > maxCheck.x());
  assert(pmax.y() + kHalfTolerance > maxCheck.y());
  assert(pmax.z() + kHalfTolerance > maxCheck.z());

  G4double magnitude = (pmax - pmin).mag();
  assert((minCheck.x() - pmin.x()) / magnitude < 0.001);
  assert((minCheck.y() - pmin.y()) / magnitude < 0.001);
  assert((minCheck.z() - pmin.z()) / magnitude < 0.001);
  assert((pmax.x() - maxCheck.x()) / magnitude < 0.001);
  assert((pmax.y() - maxCheck.y()) / magnitude < 0.001);
  assert((pmax.z() - maxCheck.z()) / magnitude < 0.001);
}

///////////////////////////////////////////////////////////////////////////////
//
// Check Inside(p)
//
void CheckInside(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {4, 5, 6, 7}};

  // set vertices in 3D
  G4double dz = solid->GetZHalfLength();
  std::vector<G4ThreeVector> vv(8);
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
    pcenter += (vv[i] + vv[k]) / 8;
  }
  // check that gravity center is inside
  assert(solid->Inside(pcenter) == kInside);
  // check positions of generated points
  auto Ns = 5, Nt = 5, k = 0;
  if (solid->GetName() == "Test_EMEC")
  {
    Ns = Nt = 4;
    k = 1;
  }
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
        if (solid->GetName() != "Test_EMEC")
        {
          assert(solid->Inside(p - 10 * shift) == kInside);
        }
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
void CheckSurfaceNormal(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // set vertices in 3D
  G4double dz = solid->GetZHalfLength();
  std::vector<G4ThreeVector> vv(8);
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
    pcenter += (vv[i] + vv[k]) / 8;
  }
  // check that gravity center is inside
  assert(solid->Inside(pcenter) == kInside);
  // check normals within surfaces
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
        if (((p2 - p0).cross(p1 - p3)).mag() == 0.) continue;  // zero area
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector n = getNormal(p0, p1, p2, p3, s, t);
        G4ThreeVector shift = (p - pcenter).unit();

        G4ThreeVector pshift = p + 0. * shift;
        ;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 0.1 * kTolerance));

        pshift = p + 0.25 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 0.1 * kTolerance));

        pshift = p - 0.25 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 0.1 * kTolerance));

        pshift = p + 0.499 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 0.1 * kTolerance));

        pshift = p - 0.499 * kTolerance * shift;
        assert(solid->Inside(pshift) == kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 0.1 * kTolerance));

        // check ApproxSurfaceNormal(), point is not on the surface
        pshift = p - 1000. * kTolerance * shift;
        assert(solid->Inside(pshift) != kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 100. * kTolerance));

        pshift = p + 1000. * kTolerance * shift;
        assert(solid->Inside(pshift) != kSurface);
        assert(ApproxEqual(n, solid->SurfaceNormal(pshift), 100. * kTolerance));
      }
    }
  }
  // check normals at lateral edges
  for (auto i = 0; i < 4; ++i)
  {
    G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                         {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};
    auto ifa1 = i + 1;
    G4ThreeVector a0 = vv[iface[ifa1][0]];
    G4ThreeVector a1 = vv[iface[ifa1][1]];
    G4ThreeVector a2 = vv[iface[ifa1][2]];
    G4ThreeVector a3 = vv[iface[ifa1][3]];
    if (((a2 - a0).cross(a1 - a3)).mag() == 0.) continue;  // zero area
    G4ThreeVector n1 = getNormal(a0, a1, a2, a3, 0.5, 1.);

    auto ifa2 = (i + 1) % 4 + 1;
    G4ThreeVector b0 = vv[iface[ifa2][0]];
    G4ThreeVector b1 = vv[iface[ifa2][1]];
    G4ThreeVector b2 = vv[iface[ifa2][2]];
    G4ThreeVector b3 = vv[iface[ifa2][3]];
    if (((b2 - b0).cross(b1 - b3)).mag() == 0.) continue;  // zero area
    G4ThreeVector n2 = getNormal(b0, b1, b2, b3, 0.5, 0.);

    G4ThreeVector pedge = 0.5 * (b0 + b1);
    G4ThreeVector nedge = (n1 + n2).unit();
    assert(solid->Inside(pedge) == kSurface);
    assert(ApproxEqual(nedge, solid->SurfaceNormal(pedge), 0.1 * kTolerance));
    ;
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check safety (DistanceToIn(p)/DistanceToOut(p))
//
void CheckSafety(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // set vertices in 3D
  G4double dz = solid->GetZHalfLength();
  std::vector<G4ThreeVector> vv(8);
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
  }
  // set points on the surface
  auto Ns = 100, Nt = 100, nsurf = 6 * Ns * Nt;
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
  G4double epsilon = kTolerance / 4.;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToIn(p) >= 0);
    assert(solid->DistanceToIn(p) < 0.01 * kTolerance);
    assert(solid->DistanceToOut(p) >= 0);
    assert(solid->DistanceToOut(p) < 0.01 * kTolerance);

    p = psurf[i] - norm * epsilon;
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToIn(p) >= 0);  // point may be located on outer side of the surface
    assert(solid->DistanceToIn(p) < epsilon);
    assert(solid->DistanceToOut(p) >= 0);
    assert(solid->DistanceToOut(p) < epsilon + 0.01 * kTolerance);

    p = psurf[i] + norm * epsilon;
    assert(solid->Inside(p) == kSurface);
    assert(solid->DistanceToOut(p) == 0);
    assert(solid->DistanceToIn(p) >= 0);
    assert(solid->DistanceToIn(p) < epsilon + 0.01 * kTolerance);
  }

  // Check safety distance for points shifted by 1 mm to inside
  G4double shift = 1. * mm;
  auto kinside = 0;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] - norm * shift;
    if (solid->Inside(p) == kInside)
    {
      kinside++;
      assert(solid->DistanceToOut(p) > 0);
      assert(solid->DistanceToOut(p) < shift + 0.01 * kTolerance);
    }
    else
    {
      assert(solid->DistanceToIn(p) > 0);
      assert(solid->DistanceToIn(p) < shift + 0.01 * kTolerance);
    }
  }
  if (solid->GetName() == "Test_EMEC")
  {
    assert(kinside > 0.80 * nsurf);  // check that at least 80% of points are inside
  }
  else
  {
    assert(kinside > 0.90 * nsurf);  // check that at least 90% of points are inside
  }

  // Check safety distance for points shifted by 10 mm to inside
  shift = 10. * mm;
  kinside = 0;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] - norm * shift;
    if (solid->Inside(p) == kInside)
    {
      kinside++;
      assert(solid->DistanceToOut(p) > 0);
      assert(solid->DistanceToOut(p) < shift + 0.01 * kTolerance);
    }
    else
    {
      assert(solid->DistanceToIn(p) > 0);
      assert(solid->DistanceToIn(p) < shift + 0.01 * kTolerance);
    }
  }
  if (solid->GetName() == "Test_EMEC")
  {
    assert(kinside != 0);  // check that at least some points are inside
  }
  else
  {
    assert(kinside > 0.80 * nsurf);  // check that at least 80% of points are inside
  }

  // Check safety distance for points shifted by 1 cm to outside
  shift = 1. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift + 0.01 * kTolerance);
  }

  // Check safety distance for points shifted by 10 cm to outside
  shift = 10. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift + 0.01 * kTolerance);
  }

  // Check safety distance for points shifted by 100 cm to outside
  shift = 100. * cm;
  for (auto i = 0; i < nsurf; ++i)
  {
    G4ThreeVector p = psurf[i];
    G4ThreeVector norm = solid->SurfaceNormal(p);
    p = psurf[i] + norm * shift;
    assert(solid->Inside(p) == kOutside);
    assert(solid->DistanceToIn(p) > 0);
    assert(solid->DistanceToIn(p) < shift + 0.01 * kTolerance);
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToIn(p,v)
// Checks directions through points on the surface
//
void CheckDistanceToIn(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double halfTolerance = 0.5 * kTolerance;
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};
  G4double shift = 0, dist = 0;

  // set 3D vertices
  G4double dz = solid->GetZHalfLength();
  G4ThreeVector vv[8];
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
  }

  // Check directions along normals.
  // Points on a face (except its edge), are moved along the normal.
  //
  auto Ns = 11, Nt = 11, ks = 1, kt = 1;

  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    // next two loops generate points on the surface
    for (auto is = ks; is < Ns - ks; ++is)
    {
      G4double s = is * sdel;
      for (auto it = kt; it < Nt - kt; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector vector = solid->SurfaceNormal(p);

        // check points on surface
        G4ThreeVector pshift = p;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        pshift = p + 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        pshift = p - 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        // shift point to outside by 1 mm
        shift = 1. * mm;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        // shift point to outside by 10 cm
        shift = 10. * cm;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        // shift point to outside by 1 m
        shift = 1. * m;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);

        // shift point to outside by 10 m
        shift = 10. * m;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        dist = solid->DistanceToIn(pshift, vector);
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);
      }
    }
  }

  // Check directions to the center of the solid
  //
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 8; ++i)
  {
    pcenter += vv[i] / 8.;
  };

  // check that gravity center is inside
  assert(solid->Inside(pcenter) == kInside);
  dist = (solid->DistanceToIn(pcenter, G4ThreeVector(1, 0, 1).unit()));
  if (dist != 0)
  {
    std::cout << "\n"
              << solid->GetName() << ": DistanceToIn(pcenter) = " << dist << ", expected 0" << "\n"
              << std::endl;
  }
  assert(solid->DistanceToIn(pcenter, G4ThreeVector(1, 0, 1).unit()) == 0);

  // set points on the surface
  Ns = 11, Nt = 11, ks = 0, kt = 0;
  if (solid->GetName() == "Test_EMEC")
  {
    Ns = 12, Nt = 12, ks = 3, kt = 3;
  }

  sdel = 1. / G4double(Ns - 1);
  tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    // next two loops generate points on the surface
    for (auto is = ks; is < Ns - ks; ++is)
    {
      G4double s = is * sdel;
      for (auto it = kt; it < Nt - kt; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector vector = (p - pcenter).unit();
        if (solid->Inside(p - vector) != kInside)
        {
          std::cout << std::endl;
          if (solid->Inside(p - vector) == kSurface)
            std::cout << solid->GetName() << ": Point is on Surface, expected Inside" << std::endl;
          if (solid->Inside(p - vector) == kOutside)
            std::cout << solid->GetName() << ": Point is Outside, expected Inside" << std::endl;
          std::cout << std::endl;
        }
        assert(solid->Inside(p - vector) == kInside);

        // check points on surface
        G4ThreeVector pshift = p;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        pshift = p + 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        pshift = p - 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToIn(pshift, -vector);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);

        // shift point to outside by 1 mm
        shift = 1. * mm;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        if (solid->GetName() == "Test_EMEC") continue;  // skip the rest for Test_EMEC

        // shift point to outside by 10 cm
        shift = 10. * cm;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        // shift point to outside by 1 m
        shift = 1. * m;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);

        // shift point to outside by 10 m
        shift = 10. * m;
        pshift = p + shift * vector;
        dist = solid->DistanceToIn(pshift, -vector);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToIn = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          if (solid->Inside(pshift) == kInside)
            std::cout << "!!! Point is Inside, expected Outside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Outside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToIn_special(p,v)
// Checks directions along the faces
//
void CheckDistanceToIn_special(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double halfTolerance = 0.5 * kTolerance;
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};
  G4double shift = 0, dist = 0;

  // set 3D vertices
  G4double dz = solid->GetZHalfLength();
  G4ThreeVector vv[8];
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
    assert(solid->Inside(vv[i]) == kSurface);
    assert(solid->Inside(vv[k]) == kSurface);
  }

  // compute central point
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 8; ++i)
  {
    pcenter += vv[i] / 8.;
  };

  // set lateral face middle points
  G4ThreeVector ff[4];
  for (auto i = 0; i < 4; ++i)
  {
    ff[i].set(0, 0, 0);
    for (auto k = 0; k < 4; ++k)
    {
      ff[i] += vv[iface[i + 1][k]] / 4.;
    }
    assert(solid->Inside(ff[i]) == kSurface);
  }

  // Check directions from the vertices to the neighbour face middle points and vice versa
  //
  for (auto i = 0; i < 8; ++i)
  {
    for (auto k = 0; k < 4; ++k)
    {
      G4bool icheck = false;
      if (i == iface[k + 1][0]) icheck = true;
      if (i == iface[k + 1][1]) icheck = true;
      if (i == iface[k + 1][2]) icheck = true;
      if (i == iface[k + 1][3]) icheck = true;
      if (!icheck) continue;
      G4ThreeVector direction = (ff[k] - vv[i]).unit();
      G4double leng = (ff[k] - vv[i]).mag();
      auto position = solid->Inside((ff[k] + vv[i]) / 2.);

      // point flies from a vertex to a middle point
      dist = solid->DistanceToIn(vv[i], direction);
      if (position == kInside)
      {
        if (dist != 0.)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0.);
      }
      if (position == kOutside)
      {
        if (std::abs(dist - leng) >= halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << leng
                    << ", discrepancy = " << dist - leng << "\n"
                    << std::endl;
        }
        assert(std::abs(dist - leng) < halfTolerance);
      }

      // point is shited by 100 m
      G4double shift = 100. * m;
      dist = solid->DistanceToIn(vv[i] - direction * shift, direction);
      if (position == kInside)
      {
        if (std::abs(dist - shift) >= halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << shift
                    << ", discrepancy = " << dist - shift << "\n"
                    << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);
      }
      G4double shift_leng = shift + leng;
      if (position == kOutside)
      {
        if (std::abs(dist - shift_leng) >= halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << shift_leng
                    << ", discrepancy = " << dist - shift_leng << "\n"
                    << std::endl;
        }
        assert(std::abs(dist - shift_leng) < halfTolerance);
      }

      // point flies from a middle point to a vertex
      dist = solid->DistanceToIn(ff[k], -direction);
      if (position == kInside)
      {
        if (dist != 0.)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0.);
      }
      if (position == kOutside)
      {
        if (dist != kInfinity)
        {
          std::cout << "\n"
                    << solid->GetName() << " vertex#" << i << " center_of_face#" << (k + 1)
                    << ": DistanceToIn = " << dist << ", expected " << kInfinity << "\n"
                    << std::endl;
        }
        assert(dist == kInfinity);
      }
    }
  }

  // Check directions between the bottom vertices and the top vertices
  //
  for (auto i = 0; i < 4; ++i)
  {
    for (auto k = 4; k < 8; ++k)
    {
      if (std::abs((k - 4) - i) == 2) continue;  // skip directions along big diagonals
      G4ThreeVector direction = (vv[k] - vv[i]).unit();
      auto position = solid->Inside((vv[k] + vv[i]) / 2.);
      G4double dist1 = solid->DistanceToIn(vv[i], direction);
      G4double dist2 = solid->DistanceToIn(vv[k], -direction);
      if (position == kSurface)
      {
        assert((dist1 + dist2) == kInfinity || (dist1 + dist2) == 2 * kInfinity);
      }
      else
      {
        assert(dist1 == dist2);
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToOut(p,v)
// Checks directions to outside through points on the surface
//
void CheckDistanceToOut(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double halfTolerance = 0.5 * kTolerance;
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};
  G4double shift = 0, dist = 0;
  G4bool validNorm = false;
  G4ThreeVector norm(0, 0, 0);

  // set 3D vertices
  G4double dz = solid->GetZHalfLength();
  G4ThreeVector vv[8];
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
  }

  // Check directions along normals.
  // Points on a face (except its edge), are moved along the normal.
  //
  auto Ns = 11, Nt = 11, ks = 1, kt = 1;
  if (solid->GetName() == "Test_EMEC")
  {
    Ns = 11, Nt = 11, ks = 2, kt = 2;
  }

  G4double sdel = 1. / G4double(Ns - 1);
  G4double tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    // next two loops generate points on the surface
    for (auto is = ks; is < Ns - ks; ++is)
    {
      G4double s = is * sdel;
      for (auto it = kt; it < Nt - kt; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        // G4ThreeVector vector = solid->SurfaceNormal(p);
        G4ThreeVector vector = getNormal(p0, p1, p2, p3, s, t);

        // check points on surface
        G4ThreeVector pshift = p;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToOut(pshift, vector, true, &validNorm, &norm);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToOut = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);
        assert(validNorm == isPlane(p0, p1, p2, p3));
        if (validNorm)
        {
          for (auto k = 0; k < 8; ++k)
          {
            assert(norm.dot(vv[k] - p) < halfTolerance);
          }
        }
        else
        {
          assert(ApproxEqual(norm, vector, 0.1 * halfTolerance));
        }

        pshift = p + 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToOut(pshift, vector, true, &validNorm, &norm);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToOut = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);
        assert(validNorm == isPlane(p0, p1, p2, p3));
        if (validNorm)
        {
          for (auto k = 0; k < 8; ++k)
          {
            assert(norm.dot(vv[k] - p) < halfTolerance);
          }
        }
        else
        {
          assert(ApproxEqual(norm, vector, 0.1 * halfTolerance));
        }

        pshift = p - 0.25 * kTolerance * vector;
        assert(solid->Inside(pshift) == kSurface);
        dist = solid->DistanceToOut(pshift, vector, true, &validNorm, &norm);
        if (dist != 0)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToOut = " << dist << ", expected " << 0.0 << "\n"
                    << std::endl;
        }
        assert(dist == 0);
        assert(validNorm == isPlane(p0, p1, p2, p3));
        if (validNorm)
        {
          for (auto k = 0; k < 8; ++k)
          {
            assert(norm.dot(vv[k] - p) < halfTolerance);
          }
        }
        else
        {
          assert(ApproxEqual(norm, vector, 0.1 * halfTolerance));
        }

        // shift point to inside by 0.5 mm
        shift = 0.5 * mm;
        pshift = p - shift * vector;
        dist = solid->DistanceToOut(pshift, vector, true, &validNorm, &norm);
        if (std::abs(dist - shift) > halfTolerance)
        {
          std::cout << "\n"
                    << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToOut = " << dist << ", expected = " << shift
                    << ", discrepancy = " << dist - shift << std::endl;
          std::cout << "Safety to out = " << solid->DistanceToOut(pshift) << std::endl;
          if (solid->Inside(pshift) == kOutside)
            std::cout << "!!! Point is Outside, expected Inside !!!" << std::endl;
          if (solid->Inside(pshift) == kSurface)
            std::cout << "!!! Point is on Surface, expected Inside !!!" << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - shift) < halfTolerance);
        assert(validNorm == isPlane(p0, p1, p2, p3));
        if (validNorm)
        {
          for (auto k = 0; k < 8; ++k)
          {
            assert(norm.dot(vv[k] - p) < halfTolerance);
          }
        }
        else
        {
          assert(ApproxEqual(norm, vector, 0.1 * halfTolerance));
        }
      }
    }
  }

  // Check directions from the center of the solid
  //
  G4ThreeVector pcenter(0, 0, 0);
  for (auto i = 0; i < 8; ++i)
  {
    pcenter += vv[i] / 8.;
  };
  assert(solid->Inside(pcenter) == kInside);

  // set points on the surface
  Ns = 11, Nt = 11, ks = 1, kt = 1;
  if (solid->GetName() == "Test_EMEC")
  {
    Ns = 12, Nt = 12, ks = 3, kt = 3;
  }

  sdel = 1. / G4double(Ns - 1);
  tdel = 1. / G4double(Nt - 1);
  for (auto i = 0; i < 6; ++i)  // loop along faces
  {
    G4ThreeVector p0 = vv[iface[i][0]];
    G4ThreeVector p1 = vv[iface[i][1]];
    G4ThreeVector p2 = vv[iface[i][2]];
    G4ThreeVector p3 = vv[iface[i][3]];
    // next two loops generate points on the surface
    for (auto is = ks; is < Ns - ks; ++is)
    {
      G4double s = is * sdel;
      for (auto it = kt; it < Nt - kt; ++it)
      {
        G4double t = it * tdel;
        G4ThreeVector p = getPoint(p0, p1, p2, p3, s, t);
        G4ThreeVector vector = (p - pcenter).unit();
        if (solid->Inside(p - vector) != kInside)
        {
          std::cout << std::endl;
          if (solid->Inside(p - vector) == kSurface)
            std::cout << solid->GetName() << ": Point is on Surface, expected Inside" << std::endl;
          if (solid->Inside(p - vector) == kOutside)
            std::cout << solid->GetName() << ": Point is Outside, expected Inside" << std::endl;
          std::cout << std::endl;
        }
        assert(solid->Inside(p - vector) == kInside);

        // check direction from center to points
        G4double leng = (p - pcenter).mag();
        dist = solid->DistanceToOut(pcenter, vector, true, &validNorm, &norm);
        if (std::abs(dist - leng) >= halfTolerance)
        {
          std::cout << std::endl;
          std::cout << solid->GetName() << " face#" << i << " (is=" << is << ", it=" << it << ")"
                    << ": DistanceToOut = " << dist << ", expected " << leng << std::endl;
          std::cout << std::endl;
        }
        assert(std::abs(dist - leng) < halfTolerance);
        assert(validNorm == isPlane(p0, p1, p2, p3));
        if (validNorm)
        {
          for (auto k = 0; k < 8; ++k)
          {
            assert(norm.dot(vv[k] - p) < halfTolerance);
          }
        }
        else
        {
          assert(ApproxEqual(norm, getNormal(p0, p1, p2, p3, s, t), 0.1 * halfTolerance));
        }
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check DistanceToIn(p.v) vs DistanceToOut(p,v)
//
void CheckInOut(const G4GenericTrap* solid)
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double halfTolerance = 0.5 * kTolerance;
  G4int iface[6][4] = {{0, 1, 2, 3}, {0, 4, 5, 1}, {1, 5, 6, 2},
                       {2, 6, 7, 3}, {3, 7, 4, 0}, {7, 6, 5, 4}};

  // set vertices in 3D
  G4double dz = solid->GetZHalfLength();
  std::vector<G4ThreeVector> vv(8);
  for (auto i = 0; i < 4; ++i)
  {
    auto k = i + 4;
    G4TwoVector a = solid->GetVertex(i);
    G4TwoVector b = solid->GetVertex(k);
    vv[i].set(a.x(), a.y(), -dz);
    vv[k].set(b.x(), b.y(), dz);
  }

  // check directions from center of face to vertices
  for (auto ifa = 0; ifa < 6; ++ifa)
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    G4ThreeVector norm = solid->SurfaceNormal(pc);
    for (auto iv = 0; iv < 4; ++iv)
    {
      G4ThreeVector vec = vv[iface[ifa][iv]] - pc;
      G4double leng = vec.mag();
      if (leng < halfTolerance) continue;
      vec = vec.unit();
      // no shift
      G4ThreeVector p = pc;
      G4double distIn = solid->DistanceToIn(p, vec);
      G4double distOut = solid->DistanceToOut(p, vec);
      assert(distIn == 0 || distIn == kInfinity);
      if (distIn == kInfinity) assert(std::abs(distOut - leng) < halfTolerance || distOut == 0.);
      if (distIn == 0.) assert(std::abs(distOut - leng) < halfTolerance);
      for (auto k = 1; k < 7; ++k)
      {
        // tiny shift to outside
        p = pc + (k / 10.) * halfTolerance * norm;
        assert(solid->Inside(p) == kSurface);
        distIn = solid->DistanceToIn(p, vec);
        distOut = solid->DistanceToOut(p, vec);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);

        // tiny shift to inside
        p = pc - (k / 10.) * halfTolerance * norm;
        assert(solid->Inside(p) == kSurface);
        distIn = solid->DistanceToIn(p, vec);
        distOut = solid->DistanceToOut(p, vec);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);
      }
    }
  }

  // check directions from center of face to edge midpoints
  for (auto ifa = 0; ifa < 6; ++ifa)
  {
    G4ThreeVector p0 = vv[iface[ifa][0]];
    G4ThreeVector p1 = vv[iface[ifa][1]];
    G4ThreeVector p2 = vv[iface[ifa][2]];
    G4ThreeVector p3 = vv[iface[ifa][3]];
    G4ThreeVector pc = (p0 + p1 + p2 + p3) / 4;  // center of the face
    G4ThreeVector norm = solid->SurfaceNormal(pc);
    for (auto iv = 0; iv < 4; ++iv)
    {
      G4ThreeVector midpoint = (vv[iface[ifa][iv]] + vv[iface[ifa][(iv + 1) % 4]]) / 2.;
      G4ThreeVector vec = midpoint - pc;
      G4double leng = vec.mag();
      if (leng < halfTolerance) continue;
      vec = vec.unit();
      // no shift
      G4ThreeVector p = pc;
      G4double distIn = solid->DistanceToIn(p, vec);
      G4double distOut = solid->DistanceToOut(p, vec);
      assert(distIn == 0 || distIn == kInfinity);
      if (distIn == kInfinity) assert(std::abs(distOut - leng) < halfTolerance || distOut == 0.);
      if (distIn == 0.) assert(std::abs(distOut - leng) < halfTolerance);
      for (auto k = 1; k < 6; ++k)
      {
        // tiny shift to outside
        p = pc + (k / 10.) * halfTolerance * norm;
        assert(solid->Inside(p) == kSurface);
        distIn = solid->DistanceToIn(p, vec);
        distOut = solid->DistanceToOut(p, vec);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);

        G4ThreeVector v = (midpoint - p).unit();
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);

        // tiny shift to inside
        p = pc - (k / 10.) * halfTolerance * norm;
        assert(solid->Inside(p) == kSurface);
        distIn = solid->DistanceToIn(p, vec);
        distOut = solid->DistanceToOut(p, vec);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);

        v = (midpoint - p).unit();
        distIn = solid->DistanceToIn(p, v);
        distOut = solid->DistanceToOut(p, v);
        assert(distIn != 0. || distOut != 0.);
        if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
        assert(solid->Inside(p + distOut * vec) == kSurface);
      }
    }
  }

  // check random points on the surface
  for (auto i = 0; i < 1000000; ++i)
  {
    G4ThreeVector pnt = solid->GetPointOnSurface();
    G4ThreeVector norm = solid->SurfaceNormal(pnt);
    G4ThreeVector vec = (norm.orthogonal()).unit();
    G4double distIn = solid->DistanceToIn(pnt, vec);
    G4double distOut = solid->DistanceToOut(pnt, vec);
    if (distIn == 0. && distOut == 0.)
    {
      std::cout << std::setprecision(18) << solid->GetName() << ":  p = " << pnt << " v = " << vec
                << std::endl;
    }
    assert(distIn != 0. || distOut != 0.);
    for (auto k = 1; k < 6; ++k)
    {
      // tiny shift to outside
      G4ThreeVector p = pnt + (k / 10.) * halfTolerance * norm;
      assert(solid->Inside(p) == kSurface);
      distIn = solid->DistanceToIn(p, vec);
      distOut = solid->DistanceToOut(p, vec);
      if (distIn == 0. && distOut == 0.)
      {
        std::cout << std::setprecision(18) << solid->GetName() << ": +k = " << k << " p = " << p
                  << " v = " << vec << std::endl;
      }
      assert(distIn != 0. || distOut != 0.);
      if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
      assert(solid->Inside(p + distOut * vec) == kSurface);

      // tiny shift to inside
      p = pnt - (k / 10.) * halfTolerance * norm;
      assert(solid->Inside(p) == kSurface);
      distIn = solid->DistanceToIn(p, vec);
      distOut = solid->DistanceToOut(p, vec);
      if (distIn == 0. && distOut == 0.)
      {
        std::cout << std::setprecision(18) << solid->GetName() << ": -k = " << k << " p = " << p
                  << " v = " << vec << std::endl;
      }
      assert(distIn != 0. || distOut != 0.);
      if (distIn != kInfinity) assert(solid->Inside(p + distIn * vec) == kSurface);
      assert(solid->Inside(p + distOut * vec) == kSurface);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetSurfaceArea()
//
void CheckGetSurfaceArea(G4GenericTrap* solid)
{
  G4double area = solid->GetSurfaceArea();
  G4double areaCheck = solid->G4VSolid::GetSurfaceArea();
  std::cout << solid->GetName() << "\tarea = " << area << ", mc_estimated = " << areaCheck
            << ", diff = " << 100 * (area - areaCheck) / area << " %" << std::endl;
  G4double diff = std::abs(area - areaCheck) / area;
  assert(std::abs(area - areaCheck) / area < 0.03);  // expected < 3%
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetCubicVolume()
//
void CheckGetCubicVolume(G4GenericTrap* solid)
{
  G4double volume = solid->GetCubicVolume();
  G4double volumeCheck = solid->G4VSolid::GetCubicVolume();
  std::cout << solid->GetName() << "\tvolume = " << volume << ", mc_estimated = " << volumeCheck
            << ", diff = " << 100 * (volume - volumeCheck) / volume << " %" << std::endl;
  G4double diff = std::abs(volume - volumeCheck) / volume;
  assert(std::abs(volume - volumeCheck) / volume < 0.01);  // expected < 1%
}

///////////////////////////////////////////////////////////////////////////////
//
// Check GetPointOnSurface()
//
void CheckGetPointOnSurface(G4GenericTrap* solid)
{
  G4double dz = solid->GetZHalfLength();
  const std::vector<G4TwoVector>& vv = solid->GetVertices();
  clock_t t;

  // compite areas
  G4double TotalArea = solid->GetSurfaceArea();
  G4double BotArea = std::abs(G4GeomTools::QuadArea(vv[0], vv[1], vv[2], vv[3]));
  G4double TopArea = std::abs(G4GeomTools::QuadArea(vv[4], vv[5], vv[6], vv[7]));
  G4double LatArea = TotalArea - BotArea - TopArea;

  // pick up random points
  G4int nbot = 0, ntop = 0, nlat = 0, nfactor = 100, ntot = TotalArea * nfactor;
  t = clock();
  for (int i = 0; i < ntot; i++)
  {
    G4ThreeVector p = solid->GetPointOnSurface();
    assert(solid->Inside(p) == kSurface);
    nbot += (p.z() == -dz);
    ntop += (p.z() == dz);
  }
  nlat = ntot - nbot - ntop;
  t = clock() - t;

  // print out info
  std::cout << "\n"
            << solid->GetName() << "\t"
            << "time : " << 1000 * (double)t / CLOCKS_PER_SEC << " ms   \t" << ntot / 1000000.
            << " million points / "
            << "time per million points : " << 1000 * ((double)t / CLOCKS_PER_SEC) * 1000000. / ntot
            << " ms" << std::endl;

  std::cout << "  Area  bottom, top, lateral, total = " << BotArea << ",\t" << TopArea << ",\t"
            << LatArea << ",\t" << TotalArea << std::endl;
  std::cout << "  N.pnt bottom, top, lateral, total = " << nbot << ",\t" << ntop << ",\t" << nlat
            << ",\t" << ntot << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
//
// Main
//
int main(int argc, char* argv[])
{
  assert(TestGenericTrap());
  std::cout << "\n    All tests for G4GenericTrap passed!\n" << std::endl;

  return 0;
}
