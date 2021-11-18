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
// Unit test for G4Ellipsoid
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

#include "G4Ellipsoid.hh"
#include "G4GeomTools.hh"

#include "G4Timer.hh"

G4bool ApproxEqual(double a, const double b) {
  return std::abs(a - b) < 1e-10;
}

G4bool ApproxEqual(const G4ThreeVector& a, const G4ThreeVector& b) {
  return (a - b).mag2() < 1e-12;
}

///////////////////////////////////////////////////////////////////////////////
//
// Estimate normal to the surface at given z and phi
//
G4ThreeVector EstimateNormal(double a, double b, double c, double z, double phi)
{
  double delta = 0.001;
  double z1    = z - delta;
  double z2    = z + delta;
  double phi1  = phi - delta;
  double phi2  = phi + delta;
  double rho1  = std::sqrt((1. + z1 / c) * (1. - z1 / c));
  double rho2  = std::sqrt((1. + z2 / c) * (1. - z2 / c));

  G4ThreeVector p1(rho1 * std::cos(phi1) * a, rho1 * std::sin(phi1) * b, z1);
  G4ThreeVector p2(rho1 * std::cos(phi2) * a, rho1 * std::sin(phi2) * b, z1);
  G4ThreeVector p3(rho2 * std::cos(phi1) * a, rho2 * std::sin(phi1) * b, z2);
  G4ThreeVector p4(rho2 * std::cos(phi2) * a, rho2 * std::sin(phi2) * b, z2);

  return ((p4 - p1).cross(p3 - p2)).unit();
}

///////////////////////////////////////////////////////////////////////////////
//
// Unit test for Ellipsoid
//
bool TestEllipsoid()
{
  double kTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  double kHalfTolerance = 0.5 * kTolerance;

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check surface area, volume and other basic methods
  //
  std::cout << "=== Check Get & Set methods" << std::endl;

  double a, b, c, zbottom, ztop;

  G4Ellipsoid solid("Test_Ellipsoid", a = 3., b = 4., c = 5., zbottom = -4.5, ztop = 3.5);
  std::cout << "Ellipsoid (" << a << ", " << b << ", " << c << ", " << zbottom << ", " << ztop << ")" << std::endl;

  assert(solid.GetDx() == a);
  assert(solid.GetDy() == b);
  assert(solid.GetDz() == c);
  assert(solid.GetSemiAxisMax(0) == a);
  assert(solid.GetSemiAxisMax(1) == b);
  assert(solid.GetSemiAxisMax(2) == c);
  assert(solid.GetZBottomCut() == zbottom);
  assert(solid.GetZTopCut() == ztop);

  solid.SetZCuts(-c - 1., c + 1);
  assert(solid.GetZBottomCut() == -c);
  assert(solid.GetZTopCut() == c);

  solid.SetZCuts(-4., 2);
  assert(solid.GetZBottomCut() == -4.);
  assert(solid.GetZTopCut() == 2.);

  solid.SetZCuts(0., 0.);
  assert(solid.GetZBottomCut() == -c);
  assert(solid.GetZTopCut() == c);

  solid.SetZCuts(zbottom, ztop);
  assert(solid.GetZBottomCut() == zbottom);
  assert(solid.GetZTopCut() == ztop);

  solid.SetSemiAxis(2., 3., 4.);
  assert(solid.GetDx() == 2.);
  assert(solid.GetDy() == 3.);
  assert(solid.GetDz() == 4.);
  assert(solid.GetZBottomCut() == -4.);
  assert(solid.GetZTopCut() == 3.5);

  solid.SetSemiAxis(a, b, c);
  solid.SetZCuts(zbottom, ztop);
  assert(solid.GetDx() == a);
  assert(solid.GetDy() == b);
  assert(solid.GetDz() == c);
  assert(solid.GetZBottomCut() == zbottom);
  assert(solid.GetZTopCut() == ztop);

  std::cout << "=== Check SreamInfo()" << std::endl;
  solid.StreamInfo(std::cout);

  // Check Surface area
  int Npoints = 1000000;
  std::cout << "=== Check GetSurfaceArea()" << std::endl;

  // check sphere
  solid.SetSemiAxis(5., 5., 5.);
  solid.SetZCuts(0., 0.);
  double area      = solid.GetSurfaceArea();
  double areaMath  = 4. * CLHEP::pi * 25.;
  double areaCheck = solid.EstimateSurfaceArea(Npoints, -1.);
  std::cout << " sphere(5) = " << area << "   exact = " << areaMath << "   mc_estimated = " << areaCheck << " ("
            << Npoints / 1000000. << " million points)" << std::endl;
  assert(std::abs(area - areaMath) < 0.01 * area);

  // check prolate spheroid
  solid.SetSemiAxis(3., 3., 5.);
  solid.SetZCuts(0., 0.);
  area      = solid.GetSurfaceArea();
  double e  = 4. / 5.;
  areaMath  = CLHEP::twopi * 3. * (3. + 5. * std::asin(e) / e);
  areaCheck = solid.EstimateSurfaceArea(Npoints, -1.);
  std::cout << " spheroid(3,3,5) = " << area << "   exact = " << areaMath << "   mc_estimated = " << areaCheck << " ("
            << Npoints / 1000000. << " million points)" << std::endl;
  assert(std::abs(area - areaMath) < 0.01 * area);

  // check oblate spheroid
  solid.SetSemiAxis(5., 5., 3.);
  solid.SetZCuts(0., 0.);
  area      = solid.GetSurfaceArea();
  areaMath  = CLHEP::twopi * 25. + CLHEP::pi * 9. * std::log((1. + e) / (1 - e)) / e;
  areaCheck = solid.EstimateSurfaceArea(Npoints, -1.);
  std::cout << " spheroid(5,5,3) = " << area << "   exact = " << areaMath << "   mc_estimated = " << areaCheck << " ("
            << Npoints / 1000000. << " million points)" << std::endl;
  assert(std::abs(area - areaMath) < 0.01 * area);

  // check ellipsoid under test
  solid.SetSemiAxis(a, b, c);
  solid.SetZCuts(zbottom, ztop);
  area      = solid.GetSurfaceArea();
  areaCheck = solid.EstimateSurfaceArea(Npoints, -1.);
  std::cout << " ellipsoid(3,4,5, -4.5,3.5) = " << area << "   mc_estimated = " << areaCheck << " ("
            << Npoints / 1000000. << " million points)" << std::endl;
  assert(std::abs(area - areaCheck) < 0.01 * area);

  // Check Cubic volume
  std::cout << "=== Check GetCubicVolume()" << std::endl;
  solid.SetSemiAxis(a, b, c);
  solid.SetZCuts(zbottom, ztop);
  double vol      = solid.GetCubicVolume();
  double volCheck = solid.EstimateCubicVolume(Npoints, 0.001);
  std::cout << " volume = " << vol << "   mc_estimated = " << volCheck << " (" << Npoints / 1000000.
            << " million points)" << std::endl;
  assert(std::abs(vol - volCheck) < 0.01 * vol);

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
  std::cout << " mc_estimated:  min = " << minCheck << " max = " << maxCheck << " (" << Npoints / 1000000.
            << " million points)" << std::endl;

  assert(std::abs(minExtent.x() - minCheck.x()) < 0.001 * std::abs(minExtent.x()));
  assert(std::abs(minExtent.y() - minCheck.y()) < 0.001 * std::abs(minExtent.y()));
  assert(minExtent.z() == minCheck.z());
  assert(std::abs(maxExtent.x() - maxCheck.x()) < 0.001 * std::abs(maxExtent.x()));
  assert(std::abs(maxExtent.y() - maxCheck.y()) < 0.001 * std::abs(maxExtent.y()));
  assert(maxExtent.z() == maxCheck.z());

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check Inside()
  //
  std::cout << "=== Check Inside(p)" << std::endl;
  solid.SetSemiAxis(a, b, c);
  solid.SetZCuts(zbottom, ztop);
  int NZ      = 30;
  int NPHI    = 30;
  double DZ   = 2. * c / NZ;
  double DPHI = CLHEP::twopi / NPHI;
  for (int iz = 0; iz < NZ; ++iz) {
    double z = -c + iz * DZ;
    for (int iphi = 0; iphi < NPHI; ++iphi) {
      double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
      double phi = iphi * DPHI;
      double rho = std::sqrt((1. + z / c) * (1. - z / c));
      double px  = rho * std::cos(phi) * a + eps;
      double py  = rho * std::sin(phi) * b + eps;
      double pz  = z + eps;
      if (z < zbottom) pz = zbottom;
      if (z > ztop) pz = ztop;
      G4ThreeVector p(px, py, pz);
      assert(solid.Inside(p) == kSurface);
      assert(solid.Inside(p * 0.999) == kInside);
      assert(solid.Inside(p * 1.001) == kOutside);
    }
  }
  assert(solid.Inside(G4ThreeVector(0., 0., zbottom)) == kSurface);
  assert(solid.Inside(G4ThreeVector(0., 0., zbottom - kTolerance)) == kOutside);
  assert(solid.Inside(G4ThreeVector(0., 0., zbottom + kTolerance)) == kInside);
  assert(solid.Inside(G4ThreeVector(0., 0., ztop)) == kSurface);
  assert(solid.Inside(G4ThreeVector(0., 0., ztop - kTolerance)) == kInside);
  assert(solid.Inside(G4ThreeVector(0., 0., ztop + kTolerance)) == kOutside);
  assert(solid.Inside(G4ThreeVector(0., 0., 0.)) == kInside);
  assert(solid.Inside(G4ThreeVector(0., 0., -c)) == kOutside);
  assert(solid.Inside(G4ThreeVector(0., 0., +c)) == kOutside);
  assert(solid.Inside(G4ThreeVector(-a, 0., 0.)) == kSurface);
  assert(solid.Inside(G4ThreeVector(-a - 2.*kTolerance, 0., 0.)) == kOutside);
  assert(solid.Inside(G4ThreeVector(-a + 2.*kTolerance, 0., 0.)) == kInside);
  assert(solid.Inside(G4ThreeVector(0., b, 0.)) == kSurface);
  assert(solid.Inside(G4ThreeVector(0., b - 2.*kTolerance, 0.)) == kInside);
  assert(solid.Inside(G4ThreeVector(0., b + 2.*kTolerance, 0.)) == kOutside);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check Normal()
  //
  std::cout << "=== Check SurfaceNormal(p)" << std::endl;
  solid.SetSemiAxis(a, b, c);
  solid.SetZCuts(zbottom, ztop);
  G4ThreeVector normal(0.);

  // Check normals on lateral surface
  NZ   = 30;
  NPHI = 30;
  DZ   = (ztop - zbottom) / NZ;
  DPHI = CLHEP::twopi / NPHI;
  for (int iz = 1; iz < NZ - 1; ++iz) {
    double z = zbottom + iz * DZ;
    for (int iphi = 0; iphi < NPHI; ++iphi) {
      double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
      double phi = iphi * DPHI;
      double rho = std::sqrt((1. + z / c) * (1. - z / c));
      double px  = rho * std::cos(phi) * a + eps;
      double py  = rho * std::sin(phi) * b + eps;
      double pz  = z + eps;
      normal     = solid.SurfaceNormal(G4ThreeVector(px, py, pz));
      assert(ApproxEqual(normal, EstimateNormal(a, b, c, pz, phi)));
    }
  }

  // Check normals at zbottom edge
  for (int iphi = 0; iphi < NPHI; ++iphi) {
    double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
    double phi = iphi * DPHI;
    double rho = std::sqrt((1. + zbottom / c) * (1. - zbottom / c));
    double px  = rho * std::cos(phi) * a + eps;
    double py  = rho * std::sin(phi) * b + eps;
    double pz  = zbottom + eps;
    normal     = solid.SurfaceNormal(G4ThreeVector(px, py, pz));
    assert(ApproxEqual(normal, (EstimateNormal(a, b, c, pz, phi) + G4ThreeVector(0., 0., -1.)).unit()));
  }

  // Check normals at ztop edge
  for (int iphi = 0; iphi < NPHI; ++iphi) {
    double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
    double phi = iphi * DPHI;
    double rho = std::sqrt((1. + ztop / c) * (1. - ztop / c));
    double px  = rho * std::cos(phi) * a + eps;
    double py  = rho * std::sin(phi) * b + eps;
    double pz  = ztop + eps;
    normal     = solid.SurfaceNormal(G4ThreeVector(px, py, pz));
    assert(ApproxEqual(normal, (EstimateNormal(a, b, c, pz, phi) + G4ThreeVector(0., 0., 1.)).unit()));
  }

  // Check normals on zbottom cut
  assert(solid.SurfaceNormal(G4ThreeVector( 0.0,  0.0,  zbottom)) == G4ThreeVector(0., 0., -1.));
  assert(solid.SurfaceNormal(G4ThreeVector( 0.5,  0.0, zbottom)) == G4ThreeVector(0., 0., -1.));
  assert(solid.SurfaceNormal(G4ThreeVector( 0.0, -0.5, zbottom)) == G4ThreeVector(0., 0., -1.));
  assert(solid.SurfaceNormal(G4ThreeVector(-0.5,  0.5, zbottom)) == G4ThreeVector(0., 0., -1.));
  assert(solid.SurfaceNormal(G4ThreeVector(-0.6, -0.6, zbottom)) == G4ThreeVector(0., 0., -1.));

  // Check normals on ztop cut
  assert(solid.SurfaceNormal(G4ThreeVector( 0.0,  0.0, ztop)) == G4ThreeVector(0., 0., 1.));
  assert(solid.SurfaceNormal(G4ThreeVector( 0.5,  0.0, ztop)) == G4ThreeVector(0., 0., 1.));
  assert(solid.SurfaceNormal(G4ThreeVector( 0.0, -0.5, ztop)) == G4ThreeVector(0., 0., 1.));
  assert(solid.SurfaceNormal(G4ThreeVector(-0.5,  0.5, ztop)) == G4ThreeVector(0., 0., 1.));
  assert(solid.SurfaceNormal(G4ThreeVector(-0.6, -0.6, ztop)) == G4ThreeVector(0., 0., 1.));

  // Check points not on surface
  normal = solid.SurfaceNormal(G4ThreeVector(0., 0., 0.));
  assert(normal.mag() == 1.);

  assert(solid.SurfaceNormal(G4ThreeVector(0.0, 0.0, zbottom + 1.)) == G4ThreeVector(0., 0., -1.));
  assert(solid.SurfaceNormal(G4ThreeVector(0.1, 0.1, zbottom + 1.)) == G4ThreeVector(0., 0., -1.));

  assert(solid.SurfaceNormal(G4ThreeVector( 0.0,  0.0, ztop - 1.)) == G4ThreeVector(0., 0., 1.));
  assert(solid.SurfaceNormal(G4ThreeVector(-0.2, -0.3, ztop - 1.)) == G4ThreeVector(0., 0., 1.));

  for (int iz = 1; iz < NZ - 1; ++iz) {
    double z = zbottom + iz * DZ;
    for (int iphi = 0; iphi < NPHI; ++iphi) {
      double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
      double phi = iphi * DPHI;
      double rho = std::sqrt((1. + z / c) * (1. - z / c));
      double px  = rho * std::cos(phi) * a + eps;
      double py  = rho * std::sin(phi) * b + eps;
      double pz  = z + eps;
      G4ThreeVector p(px, py, pz);
      assert(solid.Inside(p) == kSurface);
      assert(solid.Inside(p * 1.1) == kOutside);
      assert(solid.Inside(p * 0.8) == kInside);
      normal = solid.SurfaceNormal(p * 1.1);
      assert(ApproxEqual(normal, EstimateNormal(a, b, c, pz, phi)));
      normal = solid.SurfaceNormal(p * 0.8);
      assert(ApproxEqual(normal, EstimateNormal(a, b, c, pz, phi)));
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check SafetyToIn()
  //
  std::cout << "=== Check SafetyToIn(p)" << std::endl;
  solid.SetSemiAxis(a = 3, b = 4, c = 5);
  solid.SetZCuts(zbottom = -4.5, ztop = 3.5);

  // Check consistence with the convention
  NZ   = 30;
  NPHI = 30;
  DZ   = 2. * c / NZ;
  DPHI = CLHEP::twopi / NPHI;
  for (int iz = 0; iz < NZ; ++iz) {
    double z = -c + iz * DZ;
    for (int iphi = 0; iphi < NPHI; ++iphi) {
      double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
      double phi = iphi * DPHI;
      double rho = std::sqrt((1. + z / c) * (1. - z / c));
      double px  = rho * std::cos(phi) * a + eps;
      double py  = rho * std::sin(phi) * b + eps;
      double pz  = z + eps;
      if (z < zbottom) pz = zbottom;
      if (z > ztop) pz = ztop;
      G4ThreeVector p(px, py, pz);
      // point on surface
      assert(solid.Inside(p) == kSurface);
      double safety = solid.DistanceToIn(p);
      assert(safety >= 0. &&  safety < kHalfTolerance);
      // point is outside
      p = G4ThreeVector(px, py, pz) * 2.;
      assert(solid.Inside(p) == kOutside);
      assert(solid.DistanceToIn(p) > 0.);
      // point is inside
      p = G4ThreeVector(px, py, pz) * 0.5;
      assert(solid.Inside(p) == kInside);
      assert(solid.DistanceToIn(p) == 0.);
    }
  }

  // Check particular points to verify that the algorithm works as expected
  assert(solid.DistanceToIn(G4ThreeVector(+10, 0, 0)) == 10. - a);
  assert(solid.DistanceToIn(G4ThreeVector(-10, 0, 0)) == 10. - a);
  assert(solid.DistanceToIn(G4ThreeVector(0, +10, 0)) == 10. - b);
  assert(solid.DistanceToIn(G4ThreeVector(0, -10, 0)) == 10. - b);
  assert(solid.DistanceToIn(G4ThreeVector(0, 0, +10)) == 10. - ztop);
  assert(solid.DistanceToIn(G4ThreeVector(0, 0, -10)) == 10. + zbottom);
  assert(solid.DistanceToIn(G4ThreeVector(a, b, ztop)) > 0.);
  assert(ApproxEqual(solid.DistanceToIn(G4ThreeVector(a, b, c)), (std::sqrt(3.) - 1.) * a));

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check SafetyToOut()
  //
  std::cout << "=== Check SafetyToOut(p)" << std::endl;
  solid.SetSemiAxis(a = 3, b = 4, c = 5);
  solid.SetZCuts(zbottom = -4.5, ztop = 3.5);

  // Check consistence with the convention
  NZ   = 30;
  NPHI = 30;
  DZ   = 2. * c / NZ;
  DPHI = CLHEP::twopi / NPHI;
  for (int iz = 0; iz < NZ; ++iz) {
    double z = -c + iz * DZ;
    for (int iphi = 0; iphi < NPHI; ++iphi) {
      double eps = 0.5 * kHalfTolerance * (2. * G4UniformRand() - 1.);
      double phi = iphi * DPHI;
      double rho = std::sqrt((1. + z / c) * (1. - z / c));
      double px  = rho * std::cos(phi) * a + eps;
      double py  = rho * std::sin(phi) * b + eps;
      double pz  = z + eps;
      if (z < zbottom) pz = zbottom;
      if (z > ztop) pz = ztop;
      // point on surface
      G4ThreeVector p(px, py, pz);
      assert(solid.Inside(p) == kSurface);
      double safety = solid.DistanceToOut(p);
      //std::cout << " safety = " << safety << std::endl;
      assert(safety >= 0. &&  safety < kHalfTolerance);
      // point is outside
      p = G4ThreeVector(px, py, pz) * 2.;
      assert(solid.Inside(p) == kOutside);
      assert(solid.DistanceToOut(p) == 0.);
      // point is inside
      p = G4ThreeVector(px, py, pz) * 0.5;
      assert(solid.Inside(p) == kInside);
      assert(solid.DistanceToOut(p) > 0.);
    }
  }

  // Check particular points to verify that the algorithm works as expected
  assert(solid.DistanceToOut(G4ThreeVector(0, 0, 0)) == a);
  assert(solid.DistanceToOut(G4ThreeVector(0, 0, 2.5)) == ztop - 2.5);

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check DistanceToIn()
  //
  std::cout << "=== Check DistanceToIn(p,v)" << std::endl;
  solid.SetSemiAxis(a = 3, b = 4, c = 5);
  solid.SetZCuts(zbottom = -4.5, ztop = 3.5);
  double sctop = std::sqrt((c - ztop) * (2. - c + ztop));
  double scbot = std::sqrt((c + zbottom) * (2. - c - zbottom));
  double del   = kTolerance / 3.;

  // set coordinates for points in grid
  static const int np = 11;

  double xxx[np] = {-a - 1, -a - del, -a, -a + del, -1, 0, 1, a - del, a, a + del, a + 1};
  double yyy[np] = {-b - 1, -b - del, -b, -b + del, -3.8, 0, 3.8, b - del, b, b + del, b + 1};
  double zzz[np] = {-4.5 - 1, -4.5 - del, -4.5, -4.5 + del, -1, 0, 1, 3.5 - del, 3.5, 3.5 + del, 3.5 + 1};

  // check directions parallel to +Z
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        double dist = solid.DistanceToIn(p, G4ThreeVector(0, 0, 1));
        // Check inside points ("wrong" side)
        if (solid.Inside(p) == kInside) {
          assert(dist == 0);
          continue;
        }
        // Check relative to bounding box
        if (std::abs(p.x()) > a - kHalfTolerance || std::abs(p.y()) > b - kHalfTolerance) {
          assert(dist == kInfinity);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface) {
          if (p.z() > zbottom + 1) {
            assert(dist == kInfinity);
          } else {
            assert(dist == 0.);
          }
          continue;
        }
        // Check outside points
        if (p.z() >= 0. || (p.x() * p.x() / (a * a) + p.y() * p.y() / (b * b)) >= 1.) {
          assert(dist == kInfinity);
        } else {
          if (p.x() * p.x() / (a * a * scbot * scbot) + p.y() * p.y() / (b * b * scbot * scbot) <= 1.) {
            assert(ApproxEqual(dist, zbottom - p.z()));
          } else {
            double z = c * std::sqrt(1. - p.x() * p.x() / (a * a) - p.y() * p.y() / (b * b));
            assert(ApproxEqual(dist, -z - p.z()));
          }
        }
      }
    }
  }

  // check directions parallel to -Z
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        double dist = solid.DistanceToIn(p, G4ThreeVector(0, 0, -1));
        // Check inside points ("wrong" side)
        if (solid.Inside(p) == kInside) {
          assert(dist == 0);
          continue;
        }
        // Check relative to bounding box
        if (std::abs(p.x()) > a - kHalfTolerance || std::abs(p.y()) > b - kHalfTolerance) {
          assert(dist == kInfinity);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface) {
          if (p.z() < ztop - 1) {
            assert(dist == kInfinity);
          } else {
            assert(dist == 0.);
          }
          continue;
        }
        // Check outside points
        if (p.z() <= 0 || (p.x() * p.x() / (a * a) + p.y() * p.y() / (b * b)) >= 1.) {
          assert(dist == kInfinity);
        } else {
          if (p.x() * p.x() / (a * a * sctop * sctop) + p.y() * p.y() / (b * b * sctop * sctop) <= 1.) {
            assert(ApproxEqual(dist, p.z() - ztop));
          } else {
            double z = c * std::sqrt(1. - p.x() * p.x() / (a * a) - p.y() * p.y() / (b * b));
            assert(ApproxEqual(dist, p.z() - z));
          }
        }
      }
    }
  }

  // check directions parallel to +Y
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        double dist = solid.DistanceToIn(p, G4ThreeVector(0, 1, 0));
        // Check relative to bounding box
        if (std::abs(p.x()) > a - kHalfTolerance ||
            p.z() < zbottom + kHalfTolerance || p.z() > ztop - kHalfTolerance)
        {
          assert(dist == kInfinity);
          continue;
        }
        // Check inside points
        if (solid.Inside(p) == kInside)
        {
          assert(dist == 0);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface) {
          assert (dist == (p.y() < 0) ? 0 : kInfinity);
          continue;
        }
        // Check outside points
        if (p.y() >= 0 || (p.x() * p.x() / (a * a) + p.z() * p.z() / (c * c)) >= 1.)
        {
          assert(dist == kInfinity);
        } else {
          double y = b * std::sqrt(1. - p.x() * p.x() / (a * a) - p.z() * p.z() / (c * c));
          assert(ApproxEqual(dist, -y - p.y()));
        }
      }
    }
  }

  // check directions parallel to -X
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        double dist = solid.DistanceToIn(p, G4ThreeVector(-1, 0, 0));
        // Check relative to bounding box
        if (std::abs(p.y()) > b - kHalfTolerance ||
            p.z() < zbottom + kHalfTolerance || p.z() > ztop - kHalfTolerance)
        {
          assert(dist == kInfinity);
          continue;
        }
        // Check inside points
        if (solid.Inside(p) == kInside)
        {
          assert(dist == 0);
          continue;
        }
        // Check points on surface
        if (solid.Inside(p) == kSurface) {
          assert (dist == (p.x() > 0) ? 0 : kInfinity);
          continue;
        }
        // Check outside points
        if (p.x() <= 0 || (p.y() * p.y() / (b * b) + p.z() * p.z() / (c * c)) >= 1.) {
          assert(dist == kInfinity);
        } else {
          double x = a * std::sqrt(1. - p.y() * p.y() / (b * b) - p.z() * p.z() / (c * c));
          assert(ApproxEqual(dist, p.x() - x));
        }
      }
    }
  }

  // check far points
  double Kfar = 1.e+5;
  int nz      = 40;
  int nphi    = 36;
  for (int iz = 0; iz < nz + 1; ++iz) {
    for (int iphi = 0; iphi < nphi; ++iphi) {
      double z   = -1. + iz * (2. / nz);
      double phi = iphi * (CLHEP::twopi / nphi);
      double rho = std::sqrt(1. - z * z);
      double x   = rho * std::cos(phi);
      double y   = rho * std::sin(phi);
      if (z < zbottom / c) {
        x = a * x * zbottom / (c * z);
        y = b * y * zbottom / (c * z);
        z = zbottom;
      } else if (z > ztop / c) {
        x = a * x * ztop / (c * z);
        y = b * y * ztop / (c * z);
        z = ztop;
      } else {
        x *= a;
        y *= b;
        z *= c;
      }
      G4ThreeVector p(x, y, z);
      G4ThreeVector v = -p.unit();
      double dist = solid.DistanceToIn(Kfar * p, v);
      assert(std::abs(dist - (Kfar - 1.) * p.mag()) < kHalfTolerance);
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check DistanceToOut()
  //
  std::cout << "=== Check DistanceToOut()" << std::endl;
  solid.SetSemiAxis(a = 3, b = 4, c = 5);
  solid.SetZCuts(zbottom = -4.5, ztop = 3.5);

  // check directions parallel to +Z
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0., 0., 1.);
        bool validNorm = false;
        G4ThreeVector norm(0, 0, 0);
        double dist = solid.DistanceToOut(p, v, true, &validNorm, &norm);
        // Check if point is outside ("wrong" side)
        if (solid.Inside(p) == kOutside) {
          assert(dist == 0);
        } else {
          if (p.z() > ztop - kHalfTolerance || (p.x()*p.x()/(a*a) + p.y()*p.y()/(b * b)) >= 1.) {
            assert(dist == 0.);
            assert(validNorm && norm == solid.SurfaceNormal(p));
          } else if (solid.Inside(G4ThreeVector(p.x(), p.y(), ztop)) == kSurface) {
            assert(ApproxEqual(dist, ztop - p.z()));
            assert(validNorm && norm == G4ThreeVector(0, 0, 1));
          } else {
            double z = c * std::sqrt(1. - p.x()*p.x()/(a*a) - p.y()*p.y()/(b*b));
            assert(ApproxEqual(dist, z - p.z()));
            assert(validNorm && ApproxEqual(norm, solid.SurfaceNormal(p + dist*v)));
          }
        }
      }
    }
  }

  // check directions parallel to -Z
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0., 0.,-1.);
        bool validNorm = false;
        G4ThreeVector norm(0, 0, 0);
        double dist = solid.DistanceToOut(p, v, true, &validNorm, &norm);
        // Check if point is outside ("wrong" side)
        if (solid.Inside(p) == kOutside) {
          assert(dist == 0);
        } else {
          if (p.z() < zbottom + kHalfTolerance || (p.x()*p.x()/(a*a) + p.y()*p.y()/(b*b)) >= 1.) {
            assert(dist == 0.);
            assert(validNorm && norm == solid.SurfaceNormal(p));
          } else if (solid.Inside(G4ThreeVector(p.x(), p.y(), zbottom)) == kSurface) {
            assert(ApproxEqual(dist, p.z() - zbottom));
            assert(validNorm && norm == G4ThreeVector(0, 0,-1));
          } else {
            double z = c * std::sqrt(1. - p.x()*p.x()/(a*a) - p.y()*p.y()/(b*b));
            assert(ApproxEqual(dist, p.z() + z));
            assert(validNorm && ApproxEqual(norm, solid.SurfaceNormal(p + dist*v)));
          }
        }
      }
    }
  }

  // check directions parallel to -Y
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(0.,-1., 0.);
        bool validNorm = false;
        G4ThreeVector norm(0, 0, 0);
        double dist = solid.DistanceToOut(p, v, true, &validNorm, &norm);
        // Check if point is outside ("wrong" side)
        if (solid.Inside(p) == kOutside) {
          assert(dist == 0);
        } else {
          if ((p.x()*p.x()/(a*a) + p.z()*p.z()/(c * c)) >= 1.) {
            assert(dist == 0.);
            assert(validNorm && norm == solid.SurfaceNormal(p));
          } else {
            if (p.y() < -b + kHalfTolerance) { // check relative to  bounding box
              assert(dist == 0.);
              assert(validNorm && norm == solid.SurfaceNormal(p));
            } else {
              double y = b * std::sqrt(1. - p.x()*p.x()/(a*a) - p.z()*p.z()/(c*c));
              assert(ApproxEqual(dist, p.y() + y));
              G4ThreeVector pnew = p + dist*v;
              if (pnew.z() > zbottom + kHalfTolerance && pnew.z() < ztop - kHalfTolerance) {
                assert(validNorm && ApproxEqual(norm, solid.SurfaceNormal(pnew)));
              }
            }
          }
        }
      }
    }
  }

  // check directions parallel to +X
  for (int ix = 0; ix < np; ++ix) {
    for (int iy = 0; iy < np; ++iy) {
      for (int iz = 0; iz < np; ++iz) {
        G4ThreeVector p(xxx[ix], yyy[iy], zzz[iz]);
        G4ThreeVector v(1., 0., 0.);
        bool validNorm = false;
        G4ThreeVector norm(0, 0, 0);
        double dist = solid.DistanceToOut(p, v, true, &validNorm, &norm);
        // Check if point is outside ("wrong" side)
        if (solid.Inside(p) == kOutside) {
          assert(dist == 0);
        } else {
          if ((p.y()*p.y()/(b*b) + p.z()*p.z()/(c*c)) >= 1.) {
            assert(dist == 0.);
            assert(validNorm && norm == solid.SurfaceNormal(p));
          } else {
            if (p.x() > a - kHalfTolerance) { // check relative to bounding box
              assert(dist == 0.);
              assert(validNorm && norm == solid.SurfaceNormal(p));
            } else {
              double x = a * std::sqrt(1. - p.y()*p.y()/(b*b) - p.z()*p.z()/(c*c));
              assert(ApproxEqual(dist, x - p.x()));
              G4ThreeVector pnew = p + dist*v;
              if (pnew.z() > zbottom + kHalfTolerance && pnew.z() < ztop - kHalfTolerance) {
                assert(validNorm && ApproxEqual(norm, solid.SurfaceNormal(pnew)));
              }
            }
          }
        }
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////////
  //
  // Check PointOnSurface()
  //
  std::cout << "=== Check GetPointOnSurface()" << std::endl;
  solid.SetSemiAxis(a = 3, b = 4, c = 5);
  solid.SetZCuts(zbottom = -4.5, ztop = 3.5);
  area         = solid.GetSurfaceArea();
  double hbot  = 1. + zbottom / c;
  double htop  = 1. - ztop / c;
  double szneg = CLHEP::pi * a * b * hbot * (2. - hbot);
  double szpos = CLHEP::pi * a * b * htop * (2. - htop);
  double sside = area - szneg - szpos;

  G4Timer time;
  clock_t t;

  t = clock();
  int nzneg = 0, nzpos = 0, nside = 0, nfactor = 100000, ntot = area * nfactor;
  for (int i = 0; i < ntot; i++) {
    G4ThreeVector rndPoint = solid.GetPointOnSurface();
    assert(solid.Inside(rndPoint) == kSurface);
    if (rndPoint.z() == zbottom)
      ++nzneg;
    else if (rndPoint.z() == ztop)
      ++nzpos;
    else
      ++nside;
  }
  t = clock() - t;

  std::cout << "szneg,sside,szpos = " << szneg << ", \t" << sside << ", \t" << szpos << std::endl;
  std::cout << "nzneg,nside,nzpos = " << nzneg << ", \t" << nside << ", \t" << nzpos << std::endl;

  assert(std::abs(nzneg - szneg * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nside - sside * nfactor) < 2. * std::sqrt(ntot));
  assert(std::abs(nzpos - szpos * nfactor) < 2. * std::sqrt(ntot));

  std::cout << "Time : " << (double)t/CLOCKS_PER_SEC << " sec   " << ntot / 1000000. << " million points" << std::endl;
  std::cout << "Time per million points : " << ((double)t/CLOCKS_PER_SEC) * 1000000. / ntot << " sec" << std::endl;

  return true;
}

int main(int argc, char *argv[])
{
  assert(TestEllipsoid());
  std::cout << "\n    All tests for G4Ellipsoid passed!\n" << std::endl;

  return 0;
}
