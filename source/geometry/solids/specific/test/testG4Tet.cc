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
// Unit test for G4Tet
//
// Author: E.Tcherniaev (evgueni.tcherniaev@cern.ch)
//

// ensure asserts are compiled in
#undef NDEBUG

#include <assert.h>
#include <cmath>
#include <iomanip>

#include "geomdefs.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeomTools.hh"
#include "G4Timer.hh"

#include "G4Tet.hh"

G4bool ApproxEqual(G4double a, const G4double b) {
  return std::abs(a - b) < 1e-10;
}

G4bool ApproxEqual(const G4ThreeVector& a, const G4ThreeVector& b) {
  return (a - b).mag2() < 1e-12;
}

///////////////////////////////////////////////////////////////////////////////
//
// Unit test for Tet
//
bool TestTet()
{
  G4double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kHalfTolerance = 0.5 * kTolerance;

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check surface area, volume and other basic methods
  //
  const G4ThreeVector p0(0., 0., 0.), p1(2.1, 0., 0.), p2(0., 2.2, 0.), p3(0., 0., 2.3);
  G4Tet solid("TestTet", p0, p1, p2, p3);

  std::cout << "=== Check SetVertices()" << std::endl;
  G4ThreeVector v0, v1, v2, v3;
  solid.SetVertices(p1, p2, p3, p0);
  solid.GetVertices(v0, v1, v2, v3);
  assert(v0 == p1);
  assert(v1 == p2);
  assert(v2 == p3);
  assert(v3 == p0);
  solid.SetVertices(p0, p1, p2, p3);

  std::cout << "=== Check GetVertices()" << std::endl;

  std::vector<G4ThreeVector> pp = solid.GetVertices();
  std::cout << " number of vertices: " << pp.size() << std::endl;
  for (G4int i=0; i < pp.size(); ++i)
  {
    std::cout << "   p[" << i << "] = " << pp[i] << std::endl;
  }
  assert(pp[0] == p0);
  assert(pp[1] == p1);
  assert(pp[2] == p2);
  assert(pp[3] == p3);

  std::cout << "=== Check SreamInfo()" << std::endl;
  solid.StreamInfo(std::cout);

  // Check Surface area
  int Npoints = 1000000;
  std::cout << "=== Check GetSurfaceArea()" << std::endl;
  double area = solid.GetSurfaceArea();
  double areaMath = G4GeomTools::TriangleAreaNormal(p0, p1, p2).mag() +
                    G4GeomTools::TriangleAreaNormal(p0, p2, p3).mag() +
                    G4GeomTools::TriangleAreaNormal(p0, p3, p1).mag() +
                    G4GeomTools::TriangleAreaNormal(p1, p2, p3).mag();
  double areaCheck = solid.EstimateSurfaceArea(Npoints, -1.);
  std::cout << " area = " << area << "   exact = " << areaMath << "   mc_estimated = " << areaCheck
            << " (" << Npoints / 1000000. << " million points)" << std::endl;
  assert(area == areaMath);

  // Check Cubic volume
  std::cout << "=== Check GetCubicVolume()" << std::endl;
  G4double vol = solid.GetCubicVolume();
  G4double volMath = (p1 - p0).dot(G4GeomTools::TriangleAreaNormal(p1, p2, p3).unit()) *
                     G4GeomTools::TriangleAreaNormal(p1, p2, p3).mag() / 3.;
  G4double volCheck = solid.EstimateCubicVolume(Npoints, 0.001);
  std::cout << " volume = " << vol << "   exact = " << volMath << "   mc_estimated = " << volCheck
            << " (" << Npoints / 1000000. << " million points)" << std::endl;
  assert(vol == volMath);

  // Check Extent
  std::cout << "=== Check BoundingLimits()" << std::endl;
  G4ThreeVector minCheck(kInfinity, kInfinity, kInfinity);
  G4ThreeVector maxCheck(-kInfinity, -kInfinity, -kInfinity);
  for (int i = 0; i < Npoints; ++i) {
    G4ThreeVector p = solid.GetPointOnSurface();
    minCheck.set(std::min(p.x(), minCheck.x()), std::min(p.y(), minCheck.y()), std::min(p.z(), minCheck.z()));
    maxCheck.set(std::max(p.x(), maxCheck.x()), std::max(p.y(), maxCheck.y()), std::max(p.z(), maxCheck.z()));
  }
  G4ThreeVector minExtent, maxExtent;
  solid.BoundingLimits(minExtent, maxExtent);
  std::cout << " calculated:    min = " << minExtent << " max = " << maxExtent << std::endl;
  std::cout << " mc_estimated:  min = " << minCheck << " max = " << maxCheck
            << " (" << Npoints / 1000000. << " million points)" << std::endl;

  assert(std::abs(minExtent.x() - minCheck.x()) < 0.01);
  assert(std::abs(minExtent.y() - minCheck.y()) < 0.01);
  assert(std::abs(minExtent.z() - minCheck.z()) < 0.01);
  assert(std::abs(maxExtent.x() - maxCheck.x()) < 0.01);
  assert(std::abs(maxExtent.y() - maxCheck.y()) < 0.01);
  assert(std::abs(maxExtent.z() - maxCheck.z()) < 0.01);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check Inside()
  //
  std::cout << "=== Check Inside()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);

  G4double kin  = 0.999;
  G4double kout = 1.001;

  G4ThreeVector pC = (p0 + p1 + p2 + p3) / 4.; // center of tetrahedron

  G4ThreeVector pf0 = (p0 + p1 + p2) / 3.; // centers of faces
  G4ThreeVector pf1 = (p0 + p2 + p3) / 3.;
  G4ThreeVector pf2 = (p0 + p3 + p1) / 3.;
  G4ThreeVector pf3 = (p1 + p2 + p3) / 3.;

  G4ThreeVector pe01 = (p0 + p1) / 2.; // centers of edges
  G4ThreeVector pe02 = (p0 + p2) / 2.;
  G4ThreeVector pe03 = (p0 + p3) / 2.;
  G4ThreeVector pe12 = (p1 + p2) / 2.;
  G4ThreeVector pe13 = (p1 + p3) / 2.;
  G4ThreeVector pe23 = (p2 + p3) / 2.;

  assert(solid.Inside(p0) == kSurface);
  assert(solid.Inside(p1) == kSurface);
  assert(solid.Inside(p2) == kSurface);
  assert(solid.Inside(p3) == kSurface);

  assert(solid.Inside(pf0) == kSurface);
  assert(solid.Inside(pf1) == kSurface);
  assert(solid.Inside(pf2) == kSurface);
  assert(solid.Inside(pf3) == kSurface);

  assert(solid.Inside(pe01) == kSurface);
  assert(solid.Inside(pe02) == kSurface);
  assert(solid.Inside(pe03) == kSurface);
  assert(solid.Inside(pe12) == kSurface);
  assert(solid.Inside(pe13) == kSurface);
  assert(solid.Inside(pe23) == kSurface);

  assert(solid.Inside(pC) == kInside);
  assert(solid.Inside(pC + (p0 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (p0 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (p1 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (p2 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (p3 - pC) * kin) == kInside);

  assert(solid.Inside(pC + (pf0 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pf1 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pf2 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pf3 - pC) * kin) == kInside);

  assert(solid.Inside(pC + (pe01 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pe02 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pe03 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pe12 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pe13 - pC) * kin) == kInside);
  assert(solid.Inside(pC + (pe23 - pC) * kin) == kInside);

  assert(solid.Inside(pC + (p0 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (p0 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (p1 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (p2 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (p3 - pC) * kout) == kOutside);

  assert(solid.Inside(pC + (pf0 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pf1 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pf2 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pf3 - pC) * kout) == kOutside);

  assert(solid.Inside(pC + (pe01 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pe02 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pe03 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pe12 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pe13 - pC) * kout) == kOutside);
  assert(solid.Inside(pC + (pe23 - pC) * kout) == kOutside);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check SurfaceNormal()
  //
  std::cout << "=== Check SurfaceNormal()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);
  G4ThreeVector norm;

  // check normal at faces
  G4ThreeVector nf0 = G4GeomTools::TriangleAreaNormal(p0, p1, p2).unit();
  G4ThreeVector nf1 = G4GeomTools::TriangleAreaNormal(p0, p2, p3).unit();
  G4ThreeVector nf2 = G4GeomTools::TriangleAreaNormal(p0, p3, p1).unit();
  G4ThreeVector nf3 = G4GeomTools::TriangleAreaNormal(p1, p2, p3).unit();
  if (nf0.dot(pf0 - pC) < 0.) nf0 = -nf0;
  if (nf1.dot(pf1 - pC) < 0.) nf1 = -nf1;
  if (nf2.dot(pf2 - pC) < 0.) nf2 = -nf2;
  if (nf3.dot(pf3 - pC) < 0.) nf3 = -nf3;
  assert(solid.SurfaceNormal(pf0) == nf0);
  assert(solid.SurfaceNormal(pf1) == nf1);
  assert(solid.SurfaceNormal(pf2) == nf2);
  assert(solid.SurfaceNormal(pf3) == nf3);

  // check normal at edges
  G4ThreeVector ne01 = (nf0 + nf2).unit();
  G4ThreeVector ne02 = (nf0 + nf1).unit();
  G4ThreeVector ne03 = (nf1 + nf2).unit();
  G4ThreeVector ne12 = (nf0 + nf3).unit();
  G4ThreeVector ne13 = (nf2 + nf3).unit();
  G4ThreeVector ne23 = (nf1 + nf3).unit();
  assert(solid.SurfaceNormal(pe01) == ne01);
  assert(solid.SurfaceNormal(pe02) == ne02);
  assert(solid.SurfaceNormal(pe03) == ne03);
  assert(solid.SurfaceNormal(pe12) == ne12);
  assert(solid.SurfaceNormal(pe13) == ne13);
  assert(solid.SurfaceNormal(pe23) == ne23);

  // check normal at vertices
  G4ThreeVector n0 = (nf0 + nf1 + nf2).unit();
  G4ThreeVector n1 = (nf0 + nf2 + nf3).unit();
  G4ThreeVector n2 = (nf0 + nf1 + nf3).unit();
  G4ThreeVector n3 = (nf1 + nf2 + nf3).unit();
  assert(solid.SurfaceNormal(p0) == n0);
  assert(solid.SurfaceNormal(p1) == n1);
  assert(solid.SurfaceNormal(p2) == n2);
  assert(solid.SurfaceNormal(p3) == n3);

  // check point outside surface, ApproxSurfaceNormal()
  assert(solid.SurfaceNormal(pf0 + 0.1 * nf0) == nf0);
  assert(solid.SurfaceNormal(pf0 - 0.1 * nf0) == nf0);
  assert(solid.SurfaceNormal(pf1 + 0.1 * nf1) == nf1);
  assert(solid.SurfaceNormal(pf1 - 0.1 * nf1) == nf1);
  assert(solid.SurfaceNormal(pf2 + 0.1 * nf2) == nf2);
  assert(solid.SurfaceNormal(pf2 - 0.1 * nf2) == nf2);
  assert(solid.SurfaceNormal(pf3 + 0.1 * nf3) == nf3);
  assert(solid.SurfaceNormal(pf3 - 0.1 * nf3) == nf3);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check SafetyToIn()
  //
  std::cout << "=== Check SafetyToIn()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);

  assert(ApproxEqual(solid.DistanceToIn(pf0 + 0.5 * nf0), 0.5));
  assert(ApproxEqual(solid.DistanceToIn(pf1 + 0.5 * nf1), 0.5));
  assert(ApproxEqual(solid.DistanceToIn(pf2 + 0.5 * nf2), 0.5));
  assert(ApproxEqual(solid.DistanceToIn(pf3 + 0.5 * nf3), 0.5));

  assert(std::abs(solid.DistanceToIn(pf0)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToIn(pf1)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToIn(pf2)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToIn(pf3)) < kHalfTolerance);

  assert(solid.DistanceToIn(pf0 - 0.5 * nf0) == 0.);
  assert(solid.DistanceToIn(pf1 - 0.5 * nf1) == 0.);
  assert(solid.DistanceToIn(pf2 - 0.5 * nf2) == 0.);
  assert(solid.DistanceToIn(pf3 - 0.5 * nf3) == 0.);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check SafetyToOut()
  //
  std::cout << "=== Check SafetyToOut()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);

  assert(ApproxEqual(solid.DistanceToOut(pf0 - 0.2 * nf0), 0.2));
  assert(ApproxEqual(solid.DistanceToOut(pf1 - 0.2 * nf1), 0.2));
  assert(ApproxEqual(solid.DistanceToOut(pf2 - 0.2 * nf2), 0.2));
  assert(ApproxEqual(solid.DistanceToOut(pf3 - 0.2 * nf3), 0.2));

  assert(std::abs(solid.DistanceToOut(pf0)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToOut(pf1)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToOut(pf2)) < kHalfTolerance);
  assert(std::abs(solid.DistanceToOut(pf3)) < kHalfTolerance);

  assert(solid.DistanceToOut(pf0 + 0.5 * nf0) == 0.);
  assert(solid.DistanceToOut(pf1 + 0.5 * nf1) == 0.);
  assert(solid.DistanceToOut(pf2 + 0.5 * nf2) == 0.);
  assert(solid.DistanceToOut(pf3 + 0.5 * nf3) == 0.);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check DistanceToIn()
  //
  std::cout << "=== Check DistanceToIn()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);

  // set coordinates for points in grid
  double del   = kTolerance/3.;
  double xmin = minExtent.x();
  double ymin = minExtent.y();
  double zmin = minExtent.z();
  double xmax = maxExtent.x();
  double ymax = maxExtent.y();
  double zmax = maxExtent.z();

  static const int np = 9;
  double xxx[np] = { xmin-1, xmin-del, xmin, xmin+del, 0.5, xmax-del, xmax, xmax+del, xmax+1 };
  double yyy[np] = { ymin-1, ymin-del, ymin, ymin+del, 0.5, ymax-del, ymax, ymax+del, ymax+1 };
  double zzz[np] = { zmin-1, zmin-del, zmin, zmin+del, 0.5, zmax-del, zmax, zmax+del, zmax+1 };

  // check direction parallel to +Z
  for (int ix = 0; ix < np; ++ix)
  {
    for (int iy = 0; iy < np; ++iy)
    {
      for (int iz = 0; iz < np; ++iz)
      {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0., 0., 1.);
        double dist = solid.DistanceToIn(p, v);
        // Check inside points ("wrong" side)
        if (solid.Inside(p) == kInside)
        {
          assert(dist == 0);
          continue;
        }
        // Check relative to bounding box
        if (p.x() < xmin + kHalfTolerance || p.x() > xmax - kHalfTolerance ||
            p.y() < ymin + kHalfTolerance || p.y() > ymax - kHalfTolerance)
        {
          assert(dist == kInfinity);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface)
        {
          assert(dist == 0.);
          continue;
        }
        // Check outside points
        assert(dist == ((p.z() < zmin) ? zmin - p.z() : kInfinity));
      }
    }
  }

  // check direction parallel to -Z
  for (int ix = 0; ix < np; ++ix)
  {
    for (int iy = 0; iy < np; ++iy)
    {
      for (int iz = 0; iz < np; ++iz)
      {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0., 0.,-1.);
        double dist = solid.DistanceToIn(p, v);
        // Check inside points ("wrong" side)
        if (solid.Inside(p) == kInside)
        {
          assert(dist == 0);
          continue;
        }
        // Check relative to bounding box
        if (p.x() < xmin + kHalfTolerance || p.x() > xmax - kHalfTolerance ||
            p.y() < ymin + kHalfTolerance || p.y() > ymax - kHalfTolerance)
        {
          assert(dist == kInfinity);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface) {
          assert(dist == kInfinity);
          continue;
        }
        G4ThreeVector n = ((p2 - p1).cross(p3 - p1)).unit();
        double exp = p.z() - (n.dot(p1) - n.x()*p.x() - n.y()*p.y()) / n.z();
        assert(ApproxEqual(dist, ((p.z() < zmin) ? kInfinity : exp)));
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check DistanceToOut()
  //
  std::cout << "=== Check DistanceToOut()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);

  // check direction parallel to -Z
  for (int ix = 0; ix < np; ++ix)
  {
    for (int iy = 0; iy < np; ++iy)
    {
      for (int iz = 0; iz < np; ++iz)
      {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0., 0.,-1.);
        double dist = solid.DistanceToOut(p, v);
        assert(dist == ((p.z() <= zmin + kHalfTolerance) ? 0. : p.z() - zmin));
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check GetPointOnSurface()
  //
  std::cout << "=== Check GetPointOnSurface()" << std::endl;
  solid.SetVertices(p0, p1, p2, p3);
  G4Timer time;
  clock_t t;

  double sx =  G4GeomTools::TriangleAreaNormal(p0, p2, p3).mag();
  double sy =  G4GeomTools::TriangleAreaNormal(p0, p1, p3).mag();
  double sz =  G4GeomTools::TriangleAreaNormal(p0, p1, p2).mag();
  double sxyz = G4GeomTools::TriangleAreaNormal(p1, p2, p3).mag();

  int nx = 0, ny = 0, nz = 0, nxyz = 0, nfactor = 100000, ntot = area * nfactor;

  t = clock();
  for (int i = 0; i < ntot; i++) {
    G4ThreeVector rndPoint = solid.GetPointOnSurface();
    assert(solid.Inside(rndPoint) == kSurface);
    if (rndPoint.x() == 0.)
      ++nx;
    else if (rndPoint.y() == 0.)
      ++ny;
    else if (rndPoint.z() == 0.)
      ++nz;
    else
      ++nxyz;
  }
  t = clock() - t;

  std::cout << "sx,sy,sz,sxyz = " << sx << ",\t" << sy << ",\t" << sz << ",\t" << sxyz << std::endl;
  std::cout << "nx,ny,nz,nxyz = " << nx << ",\t" << ny << ",\t" << nz << ",\t" << nxyz << std::endl;
  assert(std::abs(nx - sx * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(ny - sy * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nz - sz * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nxyz - sxyz * nfactor) < 2. * std::sqrt(ntot));

  std::cout << "Time : " << (double)t/CLOCKS_PER_SEC << " sec   " << ntot / 1000000. << " million points" << std::endl;
  std::cout << "Time per million points : " << ((double)t/CLOCKS_PER_SEC) * 1000000. / ntot << " sec" << std::endl;

  return true;
}

int main(int argc, char *argv[])
{
  assert(TestTet());
  std::cout << "\n    All tests for G4Tet passed!\n" << std::endl;

  return 0;
}
